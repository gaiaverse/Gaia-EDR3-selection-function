#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <mpi.h>
#include <string.h>
#include <chrono>
#include <ctime>
#include <iomanip>
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
#include "Likelihood.h"
#include "timeCodes.h"
#include "GlobalVariables.h"

using Eigen::VectorXd;
//using namespace LBFGSpp;


//store MPI data as global variables
//bad practice in general, but given MPI structure, seems reasonable to break that rule
int ProcessRank;
int JobSize;
int RootID = 0; //<- declare that process 0 is always Root.


std::vector<Star> Data;
std::vector<int> Bins;
std::vector<std::string> Files;
	//RootProcess is the main action loop of the 0-ranked core. 
//It initiates the LBFGS algorithm, and controls the workflow of the other cores

void FinalResult(Eigen::VectorXd & finalpos)
{
	std::fstream outfile;
	outfile.open("testoutput.dat",std::ios::out);
	
	for (int i = 0; i < finalpos.size(); ++i)
	{
		outfile << finalpos[i] << "\n";
	}
	outfile.close();
	
}


VectorXd RootMinimiser(VectorXd &x, int steps, double lim)
{
	int nParameters = Nh+Ng*(Nt + 1);
	DescentFunctor fun(ProcessRank,Data,Bins,nParameters);
	
	std::cout << "\t New solver loop initialised" << std::endl;
	DescentFunctor::TCriteria realCriteria = DescentFunctor::TCriteria::defaults();
    cppoptlib::LbfgsSolver<DescentFunctor> solver;
	realCriteria.iterations = steps;
//	realCriteria.xDelta = 0;
	realCriteria.gradNorm = 0;
	solver.setStopCriteria(realCriteria);
	solver.minimize(fun,x);

	
	return x;
}


void RootProcess()
{
	std::cout << "\nRoot Process is intialising gradient descent framework. "; printTime();
	//tell the workers to resize their parameter vectors
	

	int nParameters = Nh+Ng*(Nt + 1);
	MPI_Bcast(&nParameters, 1, MPI_INT, RootID, MPI_COMM_WORLD);
	
	
	int nLoops = 1;
	VectorXd x = initialisedVector(nParameters);

	std::cout << "The size of the position vector is " << sizeof(x) << std::endl;

	int logStopper = -5;
	double condition = pow(10,logStopper);
	for (int i = 0; i < nLoops;++i)
	{
		
		x = RootMinimiser(x,1000000,condition);
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
    // position vector - load with initial guess, will be overwritten by the final estimate
  
	FinalResult(x);
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
	Likelihood L = Likelihood(Data,Bins,dimensionality,ProcessRank);
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
			std::cout << "Worker " << ProcessRank << " recieved the signal to end the calculation " << std::endl;
		}
	}
}


void GetAssignments(int id)
{
	std::string fileRoot = "../../MainData/";
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
		
		std::cout << "\t\t" << ProcessRank << " is opening " << file << std::endl;
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

void gradientCheck()
{
	std::fstream file;
	file.open("gradientTest.txt",std::ios::out);
	int dim = Nh + Ng*(Nt + 1);
	VectorXd y = initialisedVector(dim);
	std::cout << y.transpose() << std::endl;
	Likelihood L = Likelihood(Data,Bins,dim,ProcessRank);
	int w = 35;
	file <<std::left << std::setw(w) << "x0" FILEGAP "PriorMu" FILEGAP "AnalyticalGrad" FILEGAP "NumericalGrad" << "\n";
	
	double dm = 0.001;
	double ddm = 1e-3;
	int xsId = Nh + Ng;
	for (double xs = -1; xs < 1; xs+=dm)
	{
		VectorXd x = y;
		x[xsId] = xs;
		
		L.Calculate(x);
		double lTrue = L.Value;
		
		double dLdlm = L.Gradient[xsId];
		VectorXd xUp = x;
		VectorXd xDown = x;
		xUp[xsId] += ddm;
		xDown[xsId] -= ddm;
		
		L.Calculate(xUp);
		double lUp = L.Value;
		
		L.Calculate(xDown);
		double lDown = L.Value;
		
		double numGrad = (lUp- lTrue)/(ddm);
		
		file <<std::left << std::setw(w) << xs FILEGAP lTrue FILEGAP dLdlm FILEGAP numGrad << "\n";
	}
	
	file.close();
	
}


int main(int argc, char *argv[])
{
	//MPI initialization commands
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcessRank);
	MPI_Comm_size(MPI_COMM_WORLD, &JobSize);
	
	srand(3);
	
	auto start = std::chrono::system_clock::now();
	
	if (ProcessRank == RootID)
	{
		std::cout << "Root process online. " << JobSize - 1 << " workers connected. ";
		printTime();
		std::cout << std::endl;
	}
	
	//enter workers into their main action loops
	LoadData(ProcessRank);
	
	//gradientCheck();
	if (ProcessRank == RootID) 
	{
		RootProcess();
	}
	else
	{
		WorkerProcess();	
	}
	
	
	auto end = std::chrono::system_clock::now();
	
	std::cout << "Process " << ProcessRank << " reports job has finished. Closing MPI and exiting gracefully \n";
	
	if (ProcessRank == RootID)
	{
		std::cout << "Duration was: " << formatDuration(start,end) << "\n";
	}
	//exit gracefully
	MPI_Finalize();
	return 0;
}







