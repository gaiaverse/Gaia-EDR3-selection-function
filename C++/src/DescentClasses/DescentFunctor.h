#pragma once
#include <vector>
#include <mpi.h>
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#define EIGEN_MPL2_ONLY
#include "../libs/Eigen/Core"

#include "Star.h"
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


class DescentFunctor: public Problem<double> 
{
	private:
		int RunningID;

		//since the descent functor runs only once (i.e. on Root), we would still like root to use its CPU cycles to do some calculating, so we have a copy of 
		//the structures needed to do liklihood analysis stored within
		const std::vector<Star> &Data; 
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
		int NStars;
		//holder for transformed values
		VectorXd TransformedPosition;
		VectorXd TransformedGradient;
		

		void ForwardTransform(const VectorXd &z);
		void BackwardTransform();		
		void ResetPosition();
		
		double CurrentValue;
		VectorXd CurrentGradient;
	public:
		int LoopID;
		using typename cppoptlib::Problem<double>::Scalar;
		using typename cppoptlib::Problem<double>::TVector;
	
		double Value;
		VectorXd Gradient;
	
	    DescentFunctor(int n,const std::vector<Star> & data, int nParams,std::string outdir, int nStars): Data(data), L(data,n)
	    {
				NStars = nStars;
				RunningID = n;
				LoopID = 0;
				Start = std::chrono::system_clock::now();
				CurrentValue = 0;
				CurrentGradient = VectorXd::Zero(totalRawParams);
				
				PrevLock = VectorXd::Random(totalRawParams);
				
				TransformedPosition = VectorXd::Zero(totalTransformedParams);
				TransformedGradient = VectorXd::Zero(totalTransformedParams);
				OutputDir = outdir;
				
				std::string needlet_file = "../../ModelInputs/needlets_"+std::to_string(healpix_order)+"_"+std::to_string(needlet_order)+".csv";
				int i = 0;
			    forLineVectorInFile(needlet_file,',',
			 
					if (i > 0)
					{
				        needlet_u.push_back(std::stoi(FILE_LINE_VECTOR[0]));
				        needlet_v.push_back(std::stoi(FILE_LINE_VECTOR[1]));
				        needlet_w.push_back(std::stoi(FILE_LINE_VECTOR[2]));
					}
			        ++i;
			    );    
			    needletN = needlet_u.size();
				Value = 0;
				Gradient = VectorXd::Zero(totalRawParams);
		}
	    void DistributeCalculations(const TVector &y);
 
		double value(const TVector &x);
		void gradient(const TVector &x, TVector &grad);
		void Calculate(const VectorXd &x);
		void SavePosition(bool finalSave);
};


