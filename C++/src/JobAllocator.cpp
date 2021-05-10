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
#include "DataOperators/Star.h"
#include "DataOperators/DataLoading.h"
#include "DescentClasses/DescentFunctor.h"
#include "GenericFunctions/FileHandler.h"
#include "LikelihoodClasses/LogLikelihood.h"
#include "LikelihoodClasses/LogLikelihoodPrior.h"
#include "GenericFunctions/timeCodes.h"
#include "GlobalVariables.h"
#include "DescentClasses/ManualOptimizer.h"
#include "DataOperators/CommandArguments.h"
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
	VectorXd x = initialisedVector(nParameters,Args.LoadInStartVector,Args.StartVectorLocation);

	//initialise the functor & the solver
	DescentFunctor fun = DescentFunctor(ProcessRank,Data,totalTransformedParams,Args.OutputDirectory,TotalStars);
	Optimizer<DescentFunctor> op = Optimizer<DescentFunctor>(nParameters,fun);
	
	//set up the criteria for termination
	op.Condition.gConvergence = Args.GradLim;
	op.Condition.MaxSteps = Args.MaxSteps;
	op.Progress.ProgressDir = Args.OutputDirectory + "/";
	// GO GO GO GO!
	op.Minimize(x,N_SGD_Batches);
		
	GlobalLog(0,
		std::cout << "\nSOLVER ENDED: " << op.Converged << std::endl;
		std::cout << "\nSolver condition:\n" << op.GetStatus() << std::endl;
	);
	
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
			L.Calculate(pos,targetBatch,effectiveBatches);
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
	if (ProcessRank == RootID) 
	{
		RootProcess();
	}
	else
	{
		WorkerProcess();	
	}

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
