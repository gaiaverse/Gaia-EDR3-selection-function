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

#include "libs/LBFG/LBFGS.h"

#include "Star.h"
#include "DescentFunctor.h"
#include "FileHandler.h"
#include "Likelihood.h"
#include "timeCodes.h"
#include "GlobalVariables.h"
using Eigen::VectorXd;
using namespace LBFGSpp;


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

void FinalResult(Eigen::VectorXd finalpos)
{
	//what to do here?
}

void RootProcess()
{
	std::cout << "\nRoot Process is intialising gradient descent framework. "; printTime();
	//tell the workers to resize their parameter vectors
	

	int nParameters = Nh+Ng*(Nt + 1);
	MPI_Bcast(&nParameters, 1, MPI_INT, RootID, MPI_COMM_WORLD);
	
	
	//initialise the LBFGS parameters, then create and initialise the solver /functor pair
	LBFGSParam<double> param;
    param.epsilon = 1e-5;
    param.max_iterations = 100;
    
    LBFGSSolver<double> solver(param);
    DescentFunctor fun(ProcessRank,Data,Bins,nParameters);
	
    // position vector - load with initial guess, will be overwritten by the final estimate
    VectorXd x = VectorXd::Zero(nParameters);
    
    //initialisation of hyperhyperparameters
    std::vector<double> hyperhyper = {-1,-1,0,-1,0};
    for (int i = 0; i < Nh; ++i)
    {
		x[i] = hyperhyper[i];
	}

	//initialisation of boring old hyperparameters
	x.segment(Nh,Ng	).array() += 2.0;

    double fx;
    
    //initialise the minimization procedure
    int niter = solver.minimize(fun, x, fx);

	//broadcast to workers that the minimization procedure has finished
	int circuitBreaker = -1;
	MPI_Bcast(&circuitBreaker, 1, MPI_INT, RootID, MPI_COMM_WORLD);
	
	FinalResult(x);
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
	
	
	int RandFrac = 1e6;
	long int linesChecked = 0;
	for (int i = 0; i < Files.size(); ++i)
	{
		std::string file = Files[i];
		int gBin = Bins[i];
		//use a fancy macro (FileHandler.h) to read in data line by line, and split it into a std::vector<std::string> for the data container to process
		
		std::cout << "\t\t" << ProcessRank << " is opening " << file << std::endl;
		forLineVectorInFile(file,',',
		
			int r = rand() % RandFrac;
			if (r == 0)
			{
				Star s = Star(FILE_LINE_VECTOR,gBin);
				std::cout << ProcessRank << " got star " << linesChecked << std::endl;
				//~ Data.push_back(s);
			}	
			++linesChecked;
		);
	}
	
	auto end = std::chrono::system_clock::now();
	std::string duration = formatDuration(start,end);
	std::cout << "\tProcess " << ProcessRank << " has loaded in " << Data.size() << " datapoints in " << duration << std::endl; 
}

void gradientCheck(double dx, std::fstream & output)
{
	int dimensionality = Nh + Ng*(Nt + 1);
	Likelihood L = Likelihood(Data,Bins,dimensionality,ProcessRank);
	
	Eigen::VectorXd position = VectorXd::Random(dimensionality);
	
	position[0] = 0;
	position[1] = 0;
	position[2] = 0;
	position[3] = 0;
	position[4] = 0;
	
	std::vector<double> assign = {-1.75, -1.65, -1.55, -1.45, -1.35, -1.25, -1.15, -1.05, -0.95,
       -0.85, -0.75, -0.65, -0.55, -0.45, -0.35, -0.25, -0.15, -0.05,
        0.05,  0.15,  0.25,  0.35,  0.45,  0.55,  0.65,  0.75,  0.85,
        0.95,  1.05,  1.15,  1.25,  1.35,  1.45,  1.55,  1.65};
	for (int i = 0; i < Ng; ++i)
	{
		position[Nh + i] = assign[i];
	}
	
	//true value
	//~ std::ofstream output;
	//~ output.open("gradientTest.txt",std::ios::app);
	L.Calculate(position);
	double trueVal = L.Value;
	Eigen::VectorXd trueGrad = L.Gradient;
	Eigen::VectorXd gradApprox1 = VectorXd::Zero(dimensionality);
	Eigen::VectorXd gradApprox2 = VectorXd::Zero(dimensionality);
	int w = 20;
	//output << std::left << std::setw(w) << "ParameterValue" << std::setw(w) << "Prior" << std::setw(w)  << "TrueGrad"<< std::setw(w) << "1WayApprox"<<std::setw(w) << "delta1"<< std::setw(w) << "2WayApprox" << std::setw(w) << "delta2" << std::endl;
	for (int i = 0; i < dimensionality; ++i)
	{
		Eigen::VectorXd posUp = position;
		Eigen::VectorXd posDown = position;

		posUp[i] += dx;
		posDown[i] -= dx;
		
		L.Calculate(posUp);
		long double Lup = L.Value;
		L.Calculate(posDown);
		long double Ldown = L.Value;
		
		gradApprox1[i] = (Lup - trueVal)/dx;
		gradApprox2[i] = (Lup - Ldown)/(2*dx);
		

		//output << std::left <<std::scientific<< std::setw(w) << position[i] << std::setw(w) << trueVal << std::setw(w) <<  trueGrad[i]<< std::setw(w) << gradApprox1[i] << std::setw(w) << std::setw(w) << gradApprox2[i] /trueGrad[i] - 1 << std::setw(w) << gradApprox2[i] << std::setw(w) << gradApprox2[i] /trueGrad[i] - 1 << std::endl;
	}
	
	long double diff1 = (gradApprox1 - trueGrad).norm();
	long double diff2 = (gradApprox2 - trueGrad).norm();
	
	long double OneWayDifference = diff1 / (gradApprox1.norm() + trueGrad.norm());
	long double TwoWayDifference = diff2/ (gradApprox2.norm() + trueGrad.norm());
	
	//~ std::vector<long double> vals = {gradApprox1.norm(), gradApprox2.norm(), trueGrad.norm(), diff1, diff2, OneWayDifference, TwoWayDifference};
	//~ std::vector<std::string> names = {"OneWayGrad (norm)","TwoWayGrad (norm)","TrueGrad(norm)", "diff1", "diff2", "One-Way Difference", "Two-Way Difference"};
	
	//~ for (int i = 0; i < vals.size(); ++i)
	//~ {
		//~ output << names[i] << " is:\t" << vals[i] <<std::endl;
	//~ }
	
	output << std::left << std::setw(w) << log10(SingularityPreventer) << std::setw(w) << log10(dx) << std::setw(w) << trueVal  << std::setw(w) << TwoWayDifference <<std::endl;
	//~ output.close();
}

double SingularityPreventer = 1;

int main(int argc, char *argv[])
{
	//MPI initialization commands
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcessRank);
	MPI_Comm_size(MPI_COMM_WORLD, &JobSize);
	
	srand(ProcessRank);
	
	auto start = std::chrono::system_clock::now();
	
	if (ProcessRank == RootID)
	{
		std::cout << "Root process online. " << JobSize - 1 << " workers connected. ";
		printTime();
		std::cout << std::endl;
	}
	std::fstream output;
	int w = 20;
	output.open("gradientTest.txt",std::ios::out);
	output << std::left << std::setw(w) << "Perturber" << std::setw(w) << "dx" << std::setw(w) << "PriorMu"  << std::setw(w) << "GradScore" <<std::endl;
	for (double l = - 17; l < -8; l+=1)
	{
		std::cout << l << std::endl;
		SingularityPreventer = pow(10,l);
		for (double p = -10; p < -1; p+=1)
		{
			gradientCheck(pow(10,p),output);
		}
	}
	output.close();
	//enter workers into their main action loops
	//LoadData(ProcessRank);
	
	
	//~ if (ProcessRank == RootID) 
	//~ {
		//~ RootProcess();
	//~ }
	//~ else
	//~ {
		//~ WorkerProcess();	
	//~ }
	
	
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







