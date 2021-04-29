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
#include "DescentClasses/Star.h"
#include "DescentClasses/DescentFunctor.h"
#include "GenericFunctions/FileHandler.h"
#include "LikelihoodClasses/LogLikelihood.h"
#include "LikelihoodClasses/LogLikelihoodPrior.h"
#include "GenericFunctions/timeCodes.h"
#include "GlobalVariables.h"

using Eigen::VectorXd;



//store MPI data as global variables
//bad practice in general, but given MPI structure, seems reasonable to break that rule
int ProcessRank;
int JobSize;
int RootID = 0; //<- declare that process 0 is always Root.
int RandomSeed = time(NULL);
int burnInSteps = 1;
std::string OutputDirectory = "Output";
std::vector<Star> Data;
std::vector<int> Bins;
std::vector<std::string> Files;
	

//RootProcess is the main action loop of the 0-ranked core. 
//It initiates the LBFGS algorithm, and controls the workflow of the other cores


struct optimizerReturn
{
	VectorXd X;
	std::string Status;
	std::string Condition;
};

optimizerReturn launchMinimizer(DescentFunctor & fun, int maxSteps, double gradLim,VectorXd x)
{
	DescentFunctor::TCriteria realCriteria = DescentFunctor::TCriteria::defaults();
    cppoptlib::LbfgsSolver<DescentFunctor> solver;
	realCriteria.iterations = maxSteps;
	//~ realCriteria.xDelta = stepLim;
	realCriteria.gradNorm = gradLim;
	solver.setStopCriteria(realCriteria);
	
	
	
	//let it loose!
	solver.minimize(fun,x);
	
	optimizerReturn r;
	r.X = x;
	std::ostringstream q,p;
	
	q << solver.status();
	p << solver.criteria();
	
	r.Status = q.str();
	r.Condition = p.str();
	return r;
}

void RootProcess()
{
	GlobalLog(1,
		std::cout << "\nRoot Process is intialising gradient descent framework. "; printTime();
		std::cout << "\tAttempting to minimise " << totalRawParams << " parameters (mapped to " << totalTransformedParams << " in transform space)" << std::endl;
	);
	
	//tell the workers to resize their parameter vectors + prepare for start
	int nParameters = totalRawParams;
	int nParametersForWorkers = totalTransformedParams; 
	MPI_Bcast(&nParametersForWorkers, 1, MPI_INT, RootID, MPI_COMM_WORLD);
	
	
	VectorXd x = initialisedVector(nParameters);
	DescentFunctor fun(ProcessRank,Data,Bins,totalTransformedParams,OutputDirectory);
	
	//set up the criteria for termination

	int nLoops = 1;
	int maxSteps = 1000; 
	double gradLim = 1e-3;
	double stepLim = 1e-200;
	
	optimizerReturn r;
	for (int i = 0; i < burnInSteps;++i)
	{
		r = launchMinimizer(fun,1,gradLim,x);
		x = r.X;
		
		GlobalLog(0,
			std::cout << "\nBurnin ended ENDED: " << r.Status << std::endl;
			std::cout << "\nSolver condition:\n" << r.Condition << std::endl;
			std::cout << "\nBurnin directed to: " << x.transpose() << std::endl;
		);
	}
	
	r = launchMinimizer(fun,maxSteps,gradLim,x);
	
	x = r.X;
	GlobalLog(0,
		std::cout << "\nSOLVER ENDED: " << r.Status << std::endl;
		std::cout << "\nSolver condition:\n" << r.Condition << std::endl;

	);
	
	
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
			GlobalLog(2,
				std::cout << "\tWorker " << ProcessRank << " recieved the signal to end the calculation " << std::endl;
			);
		}
	}
}


void GetAssignments(int id)
{
	std::string fileRoot = "../../TestSets/ones/";
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
	GlobalLog(1,
		if (ProcessRank == RootID)
		{
			std::cout << "Data Readin begun" <<std::endl;
		}
	);
	auto start = std::chrono::system_clock::now();
	GetAssignments(id);
	bool isReporter = (ProcessRank == JobSize - 1);
	int meaningfullyLargeNumber = 1e8;
	int readIn = 0;
	int lastCheckPoint = 0;
	
	
	for (int i = 0; i < Files.size(); ++i)
	{
		std::string file = Files[i];
		int gBin = Bins[i];
		//use a fancy macro (FileHandler.h) to read in data line by line, and split it into a std::vector<std::string> for the data container to process
		
		forLineVectorInFile(file,',',
			Star s = Star(FILE_LINE_VECTOR,gBin);
			Data.push_back(s);
		);
	}
	
	
	
	GlobalLog(2,
		auto end = std::chrono::system_clock::now();
		std::string duration = formatDuration(start,end);
		std::cout << "\tProcess " << ProcessRank << " has loaded in " << Data.size() << " datapoints in " << duration << std::endl; 
	);
}

void processArgs(int argc, char *argv[])
{
	bool outDirFlag = false;
	bool seedFlag = false;
	bool burnFlag = false;
	for (int i = 1; i < argc; ++i)
	{
		
		std::string arg = argv[i];
		
		if (outDirFlag == true)
		{
			OutputDirectory = arg;
			outDirFlag = false;
		
		}
		if (seedFlag == true)
		{
			RandomSeed = std::stoi(arg);
			if (ProcessRank == RootID)
			{
				GlobalLog(2,
					std::cout << "Root reports random seed set to " << RandomSeed << "\n";
				);
			}
			seedFlag = false;
		}
		if (burnFlag == true)
		{
			burnInSteps = std::stoi(arg);
			if (ProcessRank == RootID)
			{
				GlobalLog(2,
					std::cout << "Root reports burnin tme set to " << burnInSteps << "\n";
				);
			}
			burnFlag == false;
		}
		
		
		if (arg == "-f")
		{
			outDirFlag = true;
		}
		if (arg == "-s")
		{
			seedFlag = true;
		}
		if (arg == "-b")
		{
			burnFlag = true;
		}
		
		
	}
	
	mkdirReturn dirReport = mkdirSafely(OutputDirectory);
	mkdirReturn dirReport2= mkdirSafely(OutputDirectory + "/" + TempDirName);
	if ((dirReport.Successful || dirReport2.Successful) == false)
	{
		ERROR(1,"Could not locate or create the output directory " + OutputDirectory + " or subdirectories therein.");
	}
}

int main(int argc, char *argv[])
{
	
	//MPI initialization commands
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcessRank);
	MPI_Comm_size(MPI_COMM_WORLD, &JobSize);

	std::cout << "hello!" << std::endl;

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

	processArgs(argc,argv);
	PrintStatus(OutputDirectory);
	srand(RandomSeed);
	
	auto start = std::chrono::system_clock::now();
	

	//pTestSuite();
	LoadData(ProcessRank);
	//mapper();
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
	
	GlobalLog(1,
		std::cout << "Process " << ProcessRank << " reports job has finished. Waiting for rest. \n";
	);
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	GlobalLog(0,
		if (ProcessRank == RootID)
		{
			std::cout << "All workers reached end of line. Duration was: " << formatDuration(start,end) << "\n";
		}
	);
	
	MPI_Finalize();
	return 0;
}







