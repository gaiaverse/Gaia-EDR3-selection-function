#pragma once
#include <vector>
#include <mpi.h>
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#define EIGEN_MPL2_ONLY
#include "../libs/Eigen/Core"

//~ #include "Star.h"
#include "../LikelihoodClasses/LogLikelihoodPrior.h"

#include "../GenericFunctions/timeCodes.h"


#include "../libs/cppoptlib/meta.h"
#include "../libs/cppoptlib/problem.h"
#include "../libs/cppoptlib/solver/bfgssolver.h"
using Eigen::VectorXd;

using namespace cppoptlib;


//DescentFunctor is a function-like class which acts as a wrapper for the gradient descent algorithm. 
//The overloaded operator () allows the class to be called as a function by LBFGs, but the classlike nature allows the function to 
//access data without needing to continually reload it into the function


class DescentFunctor
{
	private:
		int RunningID;

		//since the descent functor runs only once (i.e. on Root), we would still like root to use its CPU cycles to do some calculating, so we have a copy of 
		//the structures needed to do liklihood analysis stored within
		const std::vector<std::vector<Star>> &Data; 
		LogLikelihoodPrior L;

		
		std::chrono::time_point<std::chrono::system_clock> Start;

		//running values for the loglikelihood and gradient 

		
		std::string OutputDir;
		//to prevent double evaluations at the same point, prevlock saves current position against a threshold
		VectorXd PrevLock;
		const double lockLim = 1e-15;
		
		
		//Needlet stuff - has to be public
		int needletN;
		std::vector<int> needlet_u;
    	std::vector<int> needlet_v;
    	std::vector<double> needlet_w;
    	std::vector<double> forwardBVector;
		int NStars;
		int StarsInLastBatch;
		//holder for transformed values
		std::vector<double> TransformedPosition;
		std::vector<double> TransformedGradient;
		std::vector<double> bVector;

		void ForwardTransform(const VectorXd &z);
		void BackwardTransform();		
		void ResetPosition();
		

	public:
		int LoopID;

		double Value;
		std::vector<double> Gradient;
	
	    DescentFunctor(int n,const std::vector<std::vector<Star>> & data, int nParams,std::string outdir, int nStars): Data(data), L(data, n)
	    {
				NStars = nStars;
				RunningID = n;
				LoopID = 0;
				Start = std::chrono::system_clock::now();
								
				PrevLock = VectorXd::Random(totalRawParams);
				
				TransformedPosition = std::vector<double>(totalTransformedParams,0);
				TransformedGradient = std::vector<double>(totalTransformedParams,0);
				OutputDir = outdir;
				
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
			    needletN = needlet_u.size();
			    bVector = std::vector<double>(Nm*Ns,0);
			    
				Value = 0;
				Gradient = std::vector<double>(totalRawParams,0);
		}
	    void DistributeCalculations(const VectorXd &y, int batchID, int effectiveBatches);
 
		
		void Calculate(const VectorXd &x, int batchID, int effectiveBatches);
		void Calculate(const VectorXd &x);
		void SavePosition(bool finalSave);
};


