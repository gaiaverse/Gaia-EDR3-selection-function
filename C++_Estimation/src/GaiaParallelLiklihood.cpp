#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>

#include "libs/Eigen/Core"
#include "libs/LBFG/LBFGS.h"

#include "customClasses.h"
#include "DescentFunctor.h"
#include "FileHandler.h"
#include "Liklihood.h"
int ProcessRank;
int JobSize;
int RootID = 0;

using Eigen::VectorXd;
using namespace LBFGSpp;


std::vector<Star> Data;

void RootProcess()
{
	int nParameters = 3;
	MPI_Bcast(&nParameters, 1, MPI_INT, RootID, MPI_COMM_WORLD);
	LBFGSParam<double> param;
    param.epsilon = 1e-15;
    param.max_iterations = 30000;
 
    // Create solver and function object
    LBFGSSolver<double> solver(param);
 
    DescentFunctor fun(ProcessRank,Data,nParameters);
	
    // Initial guess
    VectorXd x = VectorXd::Zero(nParameters);
    // x will be overwritten to be the best point found
    double fx;
    int niter = solver.minimize(fun, x, fx);
	
	
	int circuitBreaker = -1;
	MPI_Bcast(&circuitBreaker, 1, MPI_INT, RootID, MPI_COMM_WORLD);
	
	std::cout << "\n\nResult: \n" << std::endl;
	std::cout << x << std::endl;
}

void WorkerProcess()
{
	//load in data
	
	
	bool hasFinished = false;
	int vecSize;
	MPI_Bcast(&vecSize, 1, MPI_INT, RootID, MPI_COMM_WORLD);
	
	std::vector<double> pos;
	pos.resize(vecSize);
	Liklihood L = Liklihood(Data,vecSize,ProcessRank);
	
	std::vector<double> emptyVec(vecSize,0.0);
	std::vector<double> gradL(vecSize,0.0);
	double emptyS=0;
	int loopChecker;
	while (hasFinished == false)
	{
		MPI_Bcast(&loopChecker, 1, MPI_INT, RootID, MPI_COMM_WORLD);

		if (loopChecker >= 0)
		{
			
			MPI_Bcast(&pos[0], vecSize, MPI_DOUBLE, RootID, MPI_COMM_WORLD);
			
			//calculate liklihood function
			//return a scalar L, and a vector of length(vecSize) of the gradient
			
			L.Calculate(pos);
			
			const double l = L.Value;
			
			
			//MPI_Reduce(&gradL, &L, vecSize + 1,MPI_DOUBLE, MPI_SUM, RootID,MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
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
	std::string fileName = "Data/MockData_" + std::to_string(id) + ".dat";
	
	forLineVectorInFile(fileName,',',
		Star s = Star(FILE_LINE_VECTOR);
		Data.push_back(s);
	);
	std::cout << "Process " << ProcessRank << " has loaded in " << Data.size() << " datapoints" << std::endl; 
}

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcessRank);
	MPI_Comm_size(MPI_COMM_WORLD, &JobSize);
	
	//printf("Rank %d online \n", ProcessRank);
	
	
	LoadData(ProcessRank);
	if (ProcessRank == RootID) 
	{
		RootProcess();
	}
	else
	{
		WorkerProcess();
	}
	MPI_Finalize();
	

	return 0;
}







