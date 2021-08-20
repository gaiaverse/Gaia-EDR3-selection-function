#pragma once
#include <vector>
#include <mpi.h>
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#define EIGEN_MPL2_ONLY
#include "../libs/Eigen/Core"

//~ #include "Star.h"
#include "../Likelihood/LogLikelihoodPrior.h"
#include "../Main/GlobalConstants.h"
#include "../libs/JSL/JSL.h"
#include "../Likelihood/ProbabilityFunctions.h"

#include "../Main/EfficiencyVector.h"
using Eigen::VectorXd;

//DescentFunctor is a function-like class which acts as a wrapper for the gradient descent algorithm. 
//The overloaded operator () allows the class to be called as a function by LBFGs, but the classlike nature allows the function to 
//access data without needing to continually reload it into the function


class DescentFunctor
{
	private:

		//since the descent functor runs only once (i.e. on Root), we would still like root to use its CPU cycles to do some calculating, so we have a copy of 
		//the structures needed to do liklihood analysis stored within
		const std::vector<std::vector<Star>> &Data; 
		LogLikelihoodPrior L;

		
		std::chrono::time_point<std::chrono::system_clock> Start;

		//running values for the loglikelihood and gradient 

		std::vector<std::vector<std::vector<double>>> HyperBuffer;
		int HyperBufferLoc;
		int HyperBufferSize = 30;
		bool HyperSaved = false;
		int hyperStep = 0;
		std::string OutputDir;		
		std::vector<double> bVector;
		
		
  
		int NStars;
		int StarsInLastBatch;
		
		EfficiencyVector Efficiency;
		
		
		void SaveHyperBuffer();
		void ResetPosition();
		
		int MaxBatches;
	public:
		int LoopID;

		double Value;
		std::vector<double> Gradient;
		
		std::vector<double> mut_gaps;

	    DescentFunctor(const std::vector<std::vector<Star>> & data, std::string outdir, std::string loadPosition, int nStars, int maxBatches): Data(data), L(data),Efficiency(loadPosition)
	    {
			
			NStars = nStars;
			LoopID = 0;
			Start = std::chrono::system_clock::now();
			MaxBatches = maxBatches;		
			
			OutputDir = outdir;
			
			Value = 0;

			HyperBuffer = std::vector(NVariancePops,std::vector(hyperOrder+2,std::vector(HyperBufferSize,0.0)));
			HyperBufferLoc = 0;
		}
	    void DistributeCalculations(const VectorXd &y, int batchID, int effectiveBatches);
 
		
		void Calculate(const VectorXd &x, int batchID, int effectiveBatches);
		void Calculate(const VectorXd &x);
		void SavePosition(bool finalSave, int saveStep, bool uniqueSave, const VectorXd & x);
		
};


