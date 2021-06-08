#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <mpi.h>
#include <string.h>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>      // std::ostringstream
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#include "../libs/Eigen/Core"
#include <fstream>
#define EIGEN_MPL2_ONLY

//~ #include "libs/LBFG/LBFGS.h"

#include "../DataHandling/Star.h"
#include "../DataHandling/DataLoading.h"
#include "../Optimizer/DescentFunctor.h"
#include "../GenericFunctions/FileHandler.h"
#include "../Likelihood/LogLikelihood.h"
#include "../Likelihood/LogLikelihoodPrior.h"
#include "../GenericFunctions/timeCodes.h"
#include "GlobalVariables.h"
#include "../Optimizer/ManualOptimizer.h"
#include "../DataHandling/CommandArguments.h"


#include "../Likelihood/ProbabilityFunctions.h"
using Eigen::VectorXd;


//store MPI data as global variables
//bad practice in general, but given MPI structure, seems reasonable to break that rule
int ProcessRank;
int JobSize;


//DATA STORAGE + Command Line Options
std::vector<std::vector<Star>> Data;
CommandArgs Args;
int TotalStars;
int MaxStarsInCore;


//RootProcess is the main action loop of the 0-ranked core. 
//It initiates the LBFGS algorithm, and controls the workflow of the other cores
VectorXd RootProcess()
{
	GlobalLog(1,
		std::cout << "\nRoot Process is intialising gradient descent framework. "; printTime();
		std::cout << "\tAttempting to minimise " << totalRawParams << " parameters (mapped to " << totalTransformedParams << " in transform space)" << std::endl;
	);
	
	//tell the workers to resize their parameter vectors + prepare for start
	int nParameters = totalRawParams;
	int nParametersForWorkers = totalTransformedParams; 
	MPI_Bcast(&nParametersForWorkers, 1, MPI_INT, RootID, MPI_COMM_WORLD);
	VectorXd x = initialisedVector(nParameters,Args.LoadInStartVector,Args.StartVectorLocation);

	//initialise the functor & the solver
	DescentFunctor fun = DescentFunctor(ProcessRank,Data,totalTransformedParams,Args.OutputDirectory,TotalStars);
	Optimizer<DescentFunctor> op = Optimizer<DescentFunctor>(nParameters,fun);
	
	//set up the criteria for termination
	op.Condition.gConvergence = Args.GradLim;
	op.Condition.MaxSteps = Args.MaxSteps;
	op.Condition.fConvergence = 5e-7;
	op.Condition.xConvergence = 0.02;
	op.Condition.SaveSteps = SaveSteps;

	op.Progress.ProgressDir = Args.OutputDirectory + "/";
	// GO GO GO GO!
	op.Minimize(x,N_SGD_Batches,Args.FreezeSteps);
		
	GlobalLog(0,
		std::cout << "\nSOLVER ENDED: " << op.Converged << std::endl;
		std::cout << "\nSolver condition:\n" << op.GetStatus() << std::endl;
	);
	
	//broadcast to workers that the minimization procedure has finished
	int circuitBreaker = -1;
	MPI_Bcast(&circuitBreaker, 1, MPI_INT, RootID, MPI_COMM_WORLD);	
	
	
	return x;
}

//this is the main action loop of all core with rank > 0 
//after initialising key objects, enters into a loop of calculting the liklihood based on a position broadcast from root
//loop continues until a circuit breaker signal is sent from root. 
void WorkerProcess()
{
	//recieve initial broadcast (this serves as a basic check of MPI functionality, rather than actually needing this data....)
	int dimensionality;
	MPI_Bcast(&dimensionality, 1, MPI_INT, RootID, MPI_COMM_WORLD);
	
	std::vector<double> pos = std::vector<double>(dimensionality,0.0);
	
	LogLikelihood L = LogLikelihood(Data,ProcessRank);
	
	//empty vectors for broadcasting reasons (do I really need these?!)
	std::vector<double> emptyVec(dimensionality,0.0);
	double emptyS=0;
	int emptyS2=0;
	int targetBatch;
	int effectiveBatches;
	bool hasFinished = false;
	while (hasFinished == false)
	{
		//check for circuit breaker signal (= -1 when circuit is broken) if no signal, calculate next iteration of L
		MPI_Bcast(&targetBatch, 1, MPI_INT, RootID, MPI_COMM_WORLD);

		if (targetBatch >= 0)
		{	

			MPI_Bcast(&effectiveBatches, 1, MPI_INT, RootID, MPI_COMM_WORLD);
			//recive new position data, copy it into position vector, then calculate likelihood contribution
			MPI_Bcast(&pos[0], dimensionality, MPI_DOUBLE, RootID, MPI_COMM_WORLD);
			
			L.Calculate(pos,targetBatch,effectiveBatches,N_SGD_Batches);
			const double l = L.Value; //for some reason, have to copy into a temporary value here - MPI fails otherwise(?)
			int nS = L.StarsUsed;
			
			//broadcast results back to root 
			MPI_Reduce(&nS, &emptyS2, 1,MPI_INT, MPI_SUM, RootID,MPI_COMM_WORLD);
			MPI_Reduce(&l,&emptyS,1,MPI_DOUBLE,MPI_SUM,RootID,MPI_COMM_WORLD);
			MPI_Reduce(&L.Gradient[0], &emptyVec[0], dimensionality,MPI_DOUBLE, MPI_SUM, RootID,MPI_COMM_WORLD);
			
		}
		else
		{
			hasFinished = true;
			GlobalLog(2,
				std::cout << "\tWorker " << ProcessRank << " recieved the signal to end the calculation " << std::endl;
			);
		}
	}
}


void PostTransform(VectorXd & z, std::vector<double> * TransformedPosition, std::vector<double> * CleanedPosition, LogLikelihoodPrior * L)
{
	//convert converged position into converged-transformed coordinates. 
	std::vector<double> mut_gaps = std::vector<double>(Nt,0);
	std::string gapFile = "../../ModelInputs/gaps_prior.dat";
	double timeFactor = (double)TotalScanningTime / Nt;
	int it = 0;
	bool inGap = false;
	int borderWidth = 0;
	int modifiedBorderWidth = borderWidth * timeFactor;
	bool inBorder= false;
	int trueTime = 0;
	int lastEnd = -9999;
	forLineVectorInFile(gapFile,' ',
					
		int gapStart = std::stoi(FILE_LINE_VECTOR[0]);
		int gapEnd = std::stoi(FILE_LINE_VECTOR[1]);
		
		trueTime = floor(it * timeFactor);
		while (trueTime < gapEnd)
		{
			int leftDistance = std::min(abs(trueTime - gapStart),abs(trueTime - lastEnd));
			int rightDistance = abs(trueTime - gapEnd);
			
			bool inGap = (trueTime >= gapStart) && (trueTime <= gapEnd);
			
			bool nearGapEdge = (leftDistance < modifiedBorderWidth) || (rightDistance < modifiedBorderWidth);
			double insertValue = mut_normal;
			if (inGap)
			{
				insertValue = mut_gap;
				//~ freezeOuts[it] = true;

			}
			if (nearGapEdge)
			{
				insertValue = mut_border;
			}
		
			mut_gaps[it] = insertValue;
			//~ std::cout << "\t " <<it << "  " << trueTime << "   " << insertValue << "   " << leftDistance << "   " << rightDistance << std::endl;
			
			++it;
			trueTime = floor((double)it * timeFactor);
			
		}
		//~ std::cout << "Gap finished at it = " << it << " t = " << trueTime << std::endl;
		lastEnd = gapEnd;
	);
	
	while (it<Nt)
	{
		mut_gaps[it] = mut_normal;
		++it;
	}
	L->MakeCovarianceMatrix();
	
	double u = exp(-1.0/lt);
	double ua = sqrt(1.0-u*u);
	double previous = z[Nt-1]; // First case is trivial
	TransformedPosition[0][Nt-1] = mut_gaps[Nt-1] + sigmat * previous;
	for (int i = Nt - 2; i >= 0; i--) 
	{
    	previous = ua * z[i] + u * previous;
    	TransformedPosition[0][i] = mut_gaps[i] + sigmat * previous;
	}
	
	std::vector<int> needlet_u;
	std::vector<int> needlet_v;
    std::vector<double> needlet_w;
	std::string needlet_file = "../../ModelInputs/needlets_"+std::to_string(healpix_order)+"_"+std::to_string(needlet_order)+".csv";
	int i = 0;
    forLineVectorInFile(needlet_file,',',
 
		if (i > 0)
		{
	        needlet_u.push_back(std::stoi(FILE_LINE_VECTOR[0]));
	        needlet_v.push_back(std::stoi(FILE_LINE_VECTOR[1]));
	        needlet_w.push_back(std::stod(FILE_LINE_VECTOR[2]));
		}
        ++i;
    );    
    int needletN = needlet_u.size();
	
	std::vector<double> bVector = std::vector<double>(Nm*Ns,0);
	for (int s = 0; s < Ns; ++s)
	{
		for (int i = 0; i < L->choleskyN; ++i)
		{
			bVector[s*Nm+L->cholesky_u[i]] += L->cholesky_w[i] * z[Nt+s*Nm+L->cholesky_v[i]];
		}
	}

	// yml
	for (int i = 0; i < needletN; ++i)
	{
		for (int m = 0; m < Nm; ++m)
		{
			TransformedPosition[0][Nt+needlet_u[i]*Nm+m] += needlet_w[i]*bVector[needlet_v[i]*Nm+m];
		}
	}
	
	for (int i = 0; i < totalTransformedParams; ++i)
	{
		CleanedPosition[0][i] = TransformedPosition[0][i];
		if (i < Nt && mut_gaps[i] == mut_gap)
		{
			CleanedPosition[0][i] = mut_gap;
		}
		
	}
	
}

void PostProcess(VectorXd x)
{
	//clean the memory
	Data.resize(1);
	Data[0].resize(1);
	MPI_Bcast(&x[0],totalRawParams,MPI_DOUBLE,RootID,MPI_COMM_WORLD);
	
	std::vector<double> TransformedPosition(totalTransformedParams,mum_prior);
	std::vector<double> CleanedPosition(totalTransformedParams,mum_prior);
	LogLikelihoodPrior L = LogLikelihoodPrior(Data,ProcessRank);
	
	PostTransform(x,&TransformedPosition,&CleanedPosition, &L);
	
	
	std::string trueDataSource = "../../Data/MainData/";
	
	std::vector<File> files = GetAssignments(ProcessRank,trueDataSource);
	std::string base = Args.OutputDirectory + "/PostProcessing/";
	if (ProcessRank == 0)
	{
		mkdirSafely(base);
	}
	std::string gapFile = "../../ModelInputs/gaps_prior.dat";
	std::vector<int> gapStarts;
	std::vector<int> gapEnds;
	forLineVectorInFile(gapFile,' ',
		gapStarts.push_back(stoi(FILE_LINE_VECTOR[0]));
		gapEnds.push_back(stoi(FILE_LINE_VECTOR[1]));
	);
	
	for (int i = 0; i < files.size(); ++i)
	{
		std::cout << ProcessRank << " is saving to " << files[i].Name << std::endl;
		std::fstream outfile;
		outfile.open(base + std::to_string(files[i].Bin + magOffset) + ".dat",std::ios::out);
		int n = 0;
		outfile << "StarID, OriginalContribution, FlattenedGap\n"; 
		forLineVectorInFile(files[i].Name,',',

			Data[0][0] = Star(FILE_LINE_VECTOR,files[i].Bin,gapStarts,gapEnds);
			
			L.Calculate(TransformedPosition,0,1,1);
			double v1 = L.Value;
			L.Calculate(CleanedPosition,0,1,1);
			double v2 = L.Value;
			
			outfile << n <<", ";
			outfile << std::setprecision(10) << v1 << ", ";
			outfile << std::setprecision(10) << v2 << "\n";
			++n;
		);
		outfile.close();
		
	}
	
}

void Welcome()
{
	GlobalLog(0,
		if (ProcessRank == RootID)
		{
			std::cout << "\n\n----------------------------------------------\n";
			std::cout << "\n~~ Gaia Selection Function Optimization ~~\n\n";
			std::cout << "Root process online. " << JobSize - 1 << " workers connected.\n";
			printTime();
			std::cout << "\n----------------------------------------------\n\n";
			std::cout << std::endl;
		}
	);
	MPI_Barrier(MPI_COMM_WORLD);
	
	PrintStatus(Args.OutputDirectory);
	srand(Args.RandomSeed);
}


void GradTest()
{
	VectorXd xV = initialisedVector(totalTransformedParams,false,"");
	std::vector<double> x;
	for (int i = 0; i < totalTransformedParams; ++i)
	{
		x.push_back(xV[i]);
	}
	std::cout << std::endl;
	//~ Star s = Data[0][0];
	int nStar = 1;
	Data.resize(1);
	Star s = Data[0][1];
	Data[0].resize(2);
	//~ Data[0][0] = s;
	
	std::cout << Data[0][0].gBin << std::endl;
	
	LogLikelihoodPrior L = LogLikelihoodPrior(Data,ProcessRank);
	
	double dx = 1e-6;
	L.Calculate(x,0,1,1);
	double trueVal = L.Value;
	std::vector<double> analyticalGrad = L.Gradient;
	int g = 0;
	int space = 0;
	for (int i = 0; i < std::min(300,totalTransformedParams); ++i)
	{
		std::vector<double> xCopy = x;
		xCopy[i] += dx;
		L.Calculate(xCopy,0,1,1);
		
		double numericalGrad = (L.Value - trueVal)/dx;
		
		if (i < Nt)
		{
			std::cout << "Time " << i;
		}
		else
		{
			std::cout << "Mag " << g << " Space " << space;
			++g;
			if (g == Nm)
			{
				++space;
				g = 0;
			}
		}
		
		std::cout << "\t"  << analyticalGrad[i] << "\t" << numericalGrad;
		 
		std::cout << std::endl;
		
	}

}

int main(int argc, char *argv[])
{
	auto start = std::chrono::system_clock::now();
	
	//MPI initialization commands
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcessRank);
	MPI_Comm_size(MPI_COMM_WORLD, &JobSize);	
	Args.ReadArguments(argc,argv,ProcessRank);
	
	Welcome();
	
	LoadData(ProcessRank,JobSize,Data,TotalStars,Args.DataSource);
	VectorXd x;
	
	//~ GradTest();
	if (ProcessRank == RootID) 
	{
		x = RootProcess();
	}
	else
	{
		x = VectorXd::Zero(totalRawParams);
		WorkerProcess();	
	}

	MPI_Barrier(MPI_COMM_WORLD);
	//~ PostProcess(x);


	//exit gracefully
	MPI_Barrier(MPI_COMM_WORLD);
	auto end = std::chrono::system_clock::now();
	GlobalLog(0,
		if (ProcessRank == RootID)
		{
			std::cout << "All workers reached end of line. Duration was: " << formatDuration(start,end) << "\n";
		}
	);
	
	MPI_Finalize();
	return 0;
}
