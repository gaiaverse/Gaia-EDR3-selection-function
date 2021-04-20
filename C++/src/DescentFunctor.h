#pragma once
#include <vector>
#include <mpi.h>
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#define EIGEN_MPL2_ONLY
#include "libs/Eigen/Core"

#include "Star.h"
#include "LogLikelihoodPrior.h"

#include "timeCodes.h"


#include "libs/cppoptlib/meta.h"
#include "libs/cppoptlib/problem.h"
#include "libs/cppoptlib/solver/bfgssolver.h"
using Eigen::VectorXd;

using namespace cppoptlib;


//DescentFunctor is a function-like class which acts as a wrapper for the gradient descent algorithm. 
//The overloaded operator () allows the class to be called as a function by LBFGs, but the classlike nature allows the function to 
//access data without needing to continually reload it into the function


class DescentFunctor: public Problem<double> 
{
	private:
		int RunningID;
		std::vector<double> InterestVectors;


		//since the descent functor runs only once (i.e. on Root), we would still like root to use its CPU cycles to do some calculating, so we have a copy of 
		//the structures needed to do liklihood analysis stored within
		const std::vector<Star> &Data; 
		LogLikelihoodPrior L;

		int LoopID;
		std::chrono::time_point<std::chrono::system_clock> Start;

		
		
		
		//running values for the loglikelihood and gradient 

		
		std::string OutputDir;
		//to prevent double evaluations at the same point, prevlock saves current position against a threshold
		VectorXd PrevLock;
		const double lockLim = 1e-15;
		
		
		//holder for transformed values
		VectorXd TransformedPosition;
		VectorXd TransformedGradient;
		

		void ForwardTransform(const VectorXd &z);
		void BackwardTransform(bool print);		
		void ResetPosition();
		
		void TestGradient(const VectorXd position);
	public:
		using typename cppoptlib::Problem<double>::Scalar;
		using typename cppoptlib::Problem<double>::TVector;
	
		double CurrentValue;
		VectorXd CurrentGradient;
	
	    DescentFunctor(int n,const std::vector<Star> & data, std::vector<int> & bins, int nParams,std::string outdir): Data(data), L(data,bins, nParams,n)
	    {
				RunningID = n;
				LoopID = 0;
				Start = std::chrono::system_clock::now();
				CurrentValue = 0;
				CurrentGradient = VectorXd::Zero(totalRawParams);
				
				PrevLock = VectorXd::Random(totalRawParams);
				
				TransformedPosition = VectorXd::Zero(totalTransformedParams);
				TransformedGradient = VectorXd::Zero(totalTransformedParams);
				OutputDir = outdir;
		}
	    void DistributeCalculations(const TVector &y, bool printOn);
 
		double value(const TVector &x);
		void gradient(const TVector &x, TVector &grad);
		
		void SavePosition(bool finalSave);
};


