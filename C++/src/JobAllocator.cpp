#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <mpi.h>
#include <string.h>
#include <chrono>
#include <ctime>
#include <iomanip>

#define EIGEN_STACK_ALLOCATION_LIMIT 0
#include "libs/Eigen/Core"
#include <fstream>
#define EIGEN_MPL2_ONLY

//~ #include "libs/LBFG/LBFGS.h"

#include "libs/cppoptlib/meta.h"
#include "libs/cppoptlib/problem.h"
#include "libs/cppoptlib/solver/gradientdescentsolver.h"
#include "libs/cppoptlib/solver/conjugatedgradientdescentsolver.h"
#include "libs/cppoptlib/solver/lbfgssolver.h"
#include "libs/cppoptlib/solver/neldermeadsolver.h"

#include "libs/cppoptlib/solver/newtondescentsolver.h"
#include "libs/cppoptlib/solver/cmaessolver.h"
#include "Star.h"
#include "DescentFunctor.h"
#include "FileHandler.h"
#include "LogLikelihood.h"
#include "LogLikelihoodPrior.h"
#include "timeCodes.h"
#include "GlobalVariables.h"

using Eigen::VectorXd;



//store MPI data as global variables
//bad practice in general, but given MPI structure, seems reasonable to break that rule
int ProcessRank;
int JobSize;
int RootID = 0; //<- declare that process 0 is always Root.

std::string OutputDirectory = "Output";
std::vector<Star> Data;
std::vector<int> Bins;
std::vector<std::string> Files;
	




VectorXd RootMinimiser(VectorXd &x, int steps, double lim,DescentFunctor &fun)
{


	DescentFunctor::TCriteria realCriteria = DescentFunctor::TCriteria::defaults();
    cppoptlib::LbfgsSolver<DescentFunctor> solver;
	realCriteria.iterations = steps;
//~ //	realCriteria.xDelta = 0;
	realCriteria.gradNorm = lim;
	solver.setStopCriteria(realCriteria);
	
	solver.minimize(fun,x);

	return x;
}

//RootProcess is the main action loop of the 0-ranked core. 
//It initiates the LBFGS algorithm, and controls the workflow of the other cores
void RootProcess()
{
	std::cout << "\nRoot Process is intialising gradient descent framework. "; printTime();
	std::cout << "\tAttempting to minimise " << totalRawParams << " (mapped to " << totalTransformedParams << " in transform space)" << std::endl;
	//tell the workers to resize their parameter vectors
	
	int nParameters = totalRawParams;
	int nParametersForWorkers = totalTransformedParams; 
	MPI_Bcast(&nParametersForWorkers, 1, MPI_INT, RootID, MPI_COMM_WORLD);
	
	
	
	int nLoops = 1;
	

	VectorXd x = initialisedVector(nParameters);
	DescentFunctor fun(ProcessRank,Data,Bins,totalTransformedParams,OutputDirectory);
	
	int logStopper = -3;
	double condition = pow(10,logStopper);
	for (int i = 0; i < nLoops;++i)
	{
		
		x = RootMinimiser(x,200,condition,fun);
		logStopper -=2;
		if (i < nLoops - 1)
		{
			condition = pow(10,logStopper);
		}
		else
		{
			condition = 0;
		}
		
	}

	fun.SavePosition(true);
	//broadcast to workers that the minimization procedure has finished
	int circuitBreaker = -1;
	MPI_Bcast(&circuitBreaker, 1, MPI_INT, RootID, MPI_COMM_WORLD);
	
}

//this is the main action loop of all core with rank > 0 
//after initialising key objects, enters into a loop of calculting the liklihood based on a position broadcast from root
//loop continues until a circuit breaker signal is sent from root. 
void WorkerProcess()
{
	//recieve initial broadcast (this serves as a basic check of MPI functionality, rather than actually needing this data....)
	int dimensionality;
	MPI_Bcast(&dimensionality, 1, MPI_INT, RootID, MPI_COMM_WORLD);
	
	
	
	//initialise the liklihood object and position vector which will be reused 
	LogLikelihood L = LogLikelihood(Data,Bins,dimensionality,ProcessRank);
	VectorXd pos = VectorXd::Zero(dimensionality);
	
	
	//empty vectors for broadcasting reasons (do I really need these?!)
	std::vector<double> emptyVec(dimensionality,0.0);
	double emptyS=0;
	

	int loopChecker;
	bool hasFinished = false;
	while (hasFinished == false)
	{
		//check for circuit breaker signal (= -1 when circuit is broken)
		MPI_Bcast(&loopChecker, 1, MPI_INT, RootID, MPI_COMM_WORLD);

		//if no signal, calculate next iteration of L
		if (loopChecker >= 0)
		{
			
			//recive new position data, copy it into position vector
			MPI_Bcast(&pos[0], dimensionality, MPI_DOUBLE, RootID, MPI_COMM_WORLD);
			
			//key bit! Likelihood updat
			L.Calculate(pos);
			const double l = L.Value; //for some reason, have to copy into a temporary value here - MPI fails otherwise(?)
			
			//broadcast results back to root - they're summed up pointwise across all workers, giving the total L and gradient contributions
			MPI_Reduce(&l,&emptyS,1,MPI_DOUBLE,MPI_SUM,RootID,MPI_COMM_WORLD);
			MPI_Reduce(&L.Gradient[0], &emptyVec[0], dimensionality,MPI_DOUBLE, MPI_SUM, RootID,MPI_COMM_WORLD);
		}
		else
		{
			hasFinished = true;
			std::cout << "\tWorker " << ProcessRank << " recieved the signal to end the calculation " << std::endl;
		}
	}
}


void GetAssignments(int id)
{
	std::string fileRoot = "../../TestSets/zeros";
	std::string assignmentFile = "coreAssignments.dat";
	
	forLineVectorInFile(assignmentFile,',',
		
		int core = stoi(FILE_LINE_VECTOR[0]);
		if (core == id)
		{
			for (int i = 1; i < FILE_LINE_VECTOR.size(); i+=2)
			{
				Files.push_back(fileRoot + FILE_LINE_VECTOR[i]);
				Bins.push_back(stoi(FILE_LINE_VECTOR[i+1]));
			}
		}
	);
}

void LoadData(int id)
{
	std::cout << "\tProcess " << ProcessRank << " beginning data readin" << std::endl;
	
	auto start = std::chrono::system_clock::now();
	GetAssignments(id);
	bool isReporter = (ProcessRank == JobSize - 1);
	int meaningfullyLargeNumber = 1e8;
	int readIn = 0;
	int lastCheckPoint = 0;
	
	
	//int RandFrac = 1e6;
	long int linesChecked = 0;
	for (int i = 0; i < Files.size(); ++i)
	{
		std::string file = Files[i];
		int gBin = Bins[i];
		//use a fancy macro (FileHandler.h) to read in data line by line, and split it into a std::vector<std::string> for the data container to process
		
		forLineVectorInFile(file,',',
		
			//~ int r = rand() % RandFrac;
			//~ if (r == 0)
			//~ {
			Star s = Star(FILE_LINE_VECTOR,gBin);

			Data.push_back(s);
			//}	
			++linesChecked;
		);
	}
	
	auto end = std::chrono::system_clock::now();
	std::string duration = formatDuration(start,end);
	std::cout << "\tProcess " << ProcessRank << " has loaded in " << Data.size() << " datapoints in " << duration << std::endl; 
}



void processArgs(int argc, char *argv[])
{
	bool outDirFlag = false;
	
	for (int i = 1; i < argc; ++i)
	{
		
		std::string arg = argv[i];
		
		if (outDirFlag == true)
		{
			OutputDirectory = arg;
			outDirFlag = false;
		
		}
		
		
		if (arg == "-f")
		{
			outDirFlag = true;
		}
		
		
	}
	
	mkdirReturn dirReport = mkdirSafely(OutputDirectory);
	mkdirReturn dirReport2= mkdirSafely(OutputDirectory + "/" + TempDirName);
	if ((dirReport.Successful || dirReport2.Successful) == false)
	{
		std::cout << "\n ERROR: Could not locate or create the output directory " << OutputDirectory << " catastrophic error.\n";
		exit(1);
	}
}


int main(int argc, char *argv[])
{
	//MPI initialization commands
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcessRank);
	MPI_Comm_size(MPI_COMM_WORLD, &JobSize);
	
	processArgs(argc,argv);
	srand(3);
	
	auto start = std::chrono::system_clock::now();
	
	if (ProcessRank == RootID)
	{
		std::cout << "Root process online. " << JobSize - 1 << " workers connected. ";
		printTime();
		std::cout << std::endl;
	}
	

	LoadData(ProcessRank);
	if (ProcessRank == RootID) 
	{
		RootProcess();
	}
	else
	{
		WorkerProcess();	
	}
	
	
	
	//exit gracefully
	auto end = std::chrono::system_clock::now();
	std::cout << "Process " << ProcessRank << " reports job has finished. Closing MPI and exiting gracefully \n";
	
	if (ProcessRank == RootID)
	{
		std::cout << "Duration was: " << formatDuration(start,end) << "\n";
	}

	MPI_Finalize();
	return 0;
}







