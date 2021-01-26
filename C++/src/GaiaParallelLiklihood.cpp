#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <mpi.h>
#include <string.h>
#include <chrono>
#include <ctime>
#include "libs/Eigen/Core"
#define EIGEN_MPL2_ONLY

#include "libs/LBFG/LBFGS.h"

#include "customClasses.h"
#include "DescentFunctor.h"
#include "FileHandler.h"
#include "Liklihood.h"
#include "timeCodes.h"
using Eigen::VectorXd;
using namespace LBFGSpp;


//store MPI data as global variables
//bad practice in general, but given MPI structure, seems reasonable to break that rule
int ProcessRank;
int JobSize;
int RootID = 0; //<- declare that process 0 is always Root.


std::vector<Star> Data;


//RootProcess is the main action loop of the 0-ranked core. 
//It initiates the LBFGS algorithm, and controls the workflow of the other cores

void RootProcess()
{
	std::cout << "\nRoot Process is intialising gradient descent framework. "; printTime();
	//tell the workers to resize their parameter vectors
	int nParameters = 4;
	MPI_Bcast(&nParameters, 1, MPI_INT, RootID, MPI_COMM_WORLD);
	
	
	//initialise the LBFGS parameters, then create and initialise the solver /functor pair
	LBFGSParam<double> param;
    param.epsilon = 1e-5;
    param.max_iterations = 100;
    
    LBFGSSolver<double> solver(param);
    DescentFunctor fun(ProcessRank,Data,nParameters);
	
    // position vector - load with initial guess, will be overwritten by the final estimate
    VectorXd x = VectorXd::Zero(nParameters);
    double fx;
    
    //initialise the minimization procedure
    int niter = solver.minimize(fun, x, fx);

	//broadcast to workers that the minimization procedure has finished
	int circuitBreaker = -1;
	MPI_Bcast(&circuitBreaker, 1, MPI_INT, RootID, MPI_COMM_WORLD);
	
	std::cout << "\n\nResult: \n" << std::endl;
	std::cout << x << std::endl;
}

//this is the main action loop of all core with rank > 0 
//after initialising key objects, enters into a loop of calculting the liklihood based on a position broadcast from root
//loop continues until a circuit breaker signal is sent from root. 
void WorkerProcess()
{
	//recieve initial broadcast (this serves as a basic check of MPI functionality, rather than actually needing this data....)
	int vecSize;
	MPI_Bcast(&vecSize, 1, MPI_INT, RootID, MPI_COMM_WORLD);
	
	
	
	//initialise the liklihood object and position vector which will be reused 
	Liklihood L = Liklihood(Data,vecSize,ProcessRank);
	VectorXd pos = VectorXd::Zero(vecSize);
	
	
	//empty vectors for broadcasting reasons (do I really need these?!)
	std::vector<double> emptyVec(vecSize,0.0);
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
			MPI_Bcast(&pos[0], vecSize, MPI_DOUBLE, RootID, MPI_COMM_WORLD);
			
			//key bit! Liklihood updat
			L.Calculate(pos);
			const double l = L.Value; //for some reason, have to copy into a temporary value here - MPI fails otherwise(?)
			
			//broadcast results back to root - they're summed up pointwise across all workers, giving the total L and gradient contributions
			MPI_Reduce(&l,&emptyS,1,MPI_DOUBLE,MPI_SUM,RootID,MPI_COMM_WORLD);
			MPI_Reduce(&L.Gradient[0], &emptyVec[0], vecSize,MPI_DOUBLE, MPI_SUM, RootID,MPI_COMM_WORLD);
		}
		else
		{
			hasFinished = true;
		}
	}
}


void LoadData(int id)
{
	//read in the data assigned to this worker - for the full code will need to work out the assignment protocols
	std::string fileName = "Data/MockData_" + std::to_string(id) + ".dat";	
	
	std::cout << "\tProcess " << ProcessRank << " beginning data readin from " << fileName<< std::endl;
	
	auto start = std::chrono::system_clock::now();

	bool isReporter = (ProcessRank == JobSize - 1);
	int meaningfullyLargeNumber = 1e8;
	int readIn = 0;
	int lastCheckPoint = 0;
	//use a fancy macro (FileHandler.h) to read in data line by line, and split it into a std::vector<std::string> for the data container to process
	forLineVectorInFile(fileName,',',
		Star s = Star(FILE_LINE_VECTOR);
		Data.push_back(s);
		
		if (isReporter)
		{
			++readIn;
			if (readIn >= lastCheckPoint + meaningfullyLargeNumber)
			{
				lastCheckPoint = readIn;
				auto checkpoint = std::chrono::system_clock::now();
				std::string duration = formatDuration(start,checkpoint);
				std::cout << "\t\tProcess " << ProcessRank << " has found " << readIn << " datapoints after " << duration << std::endl; 
			}
		}
	);
	
	auto end = std::chrono::system_clock::now();
	std::string duration = formatDuration(start,end);
	std::cout << "\tProcess " << ProcessRank << " has loaded in " << Data.size() << " datapoints in " << duration << std::endl; 
}

int main(int argc, char *argv[])
{
	//MPI initialization commands
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcessRank);
	MPI_Comm_size(MPI_COMM_WORLD, &JobSize);
	
	
	
	
	auto start = std::chrono::system_clock::now();
	std::vector<double> v;
	for (int i = 0; i < 1e8; ++i)
	{
		v.push_back(exp(-i*0.54));
	}
	std::cout << v.size() << std::endl;
	
	//~ if (ProcessRank == RootID)
	//~ {
		//~ std::cout << "Root process online. " << JobSize - 1 << " workers connected. ";
		//~ printTime();
		//~ std::cout << std::endl;
	//~ }
	//~ MPI_Barrier(MPI_COMM_WORLD);
	
	//~ //enter workers into their main action loops
	//~ LoadData(ProcessRank);
	//~ if (ProcessRank == RootID) 
	//~ {
		//~ RootProcess();
	//~ }
	//~ else
	//~ {
		//~ WorkerProcess();
	//~ }
	
	auto end = std::chrono::system_clock::now();
	if (ProcessRank == RootID)
	{
		std::cout << "Root process reports job has finished. \n";
		std::cout << "Duration was: " << formatDuration(start,end) << "\n";
		std::cout << "Closing MPI and exiting gracefully" << std::endl;
	}
	//exit gracefully
	MPI_Finalize();
	return 0;
}







