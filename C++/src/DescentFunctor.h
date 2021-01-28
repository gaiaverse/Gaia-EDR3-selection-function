#pragma once
#include <vector>
#include <mpi.h>
#include "libs/Eigen/Core"
#define EIGEN_MPL2_ONLY
#include "Star.h"
#include "Likelihood.h"
#include "timeCodes.h"
using Eigen::VectorXd;



//DescentFunctor is a function-like class which acts as a wrapper for the gradient descent algorithm. 
//The overloaded operator () allows the class to be called as a function by LBFGs, but the classlike nature allows the function to 
//access data without needing to continually reload it into the function
class DescentFunctor
{
	private:
		int RunningID;
		std::vector<double> InterestVectors;
		//since the descent functor runs only once (i.e. on Root), we would still like root to use its CPU cycles to do some calculating, so we have a copy of 
		//the structures needed to do liklihood analysis stored within
		const std::vector<Star> &Data; 
		Likelihood L;
		int LoopID;
		std::chrono::time_point<std::chrono::system_clock> Start;
		void ExamineInterestVectors(Eigen::VectorXd &position);
		
	public:
	    DescentFunctor(int n,const std::vector<Star> & data, std::vector<int> & bins, int nParams) : Data(data), L(data,bins, nParams,n) //initializer list (complicated, not really sure what it is, but it needs to be here)
	    {
				RunningID = n;
				LoopID = 0;
				Start = std::chrono::system_clock::now();
		}
	    double operator()(Eigen::VectorXd& x, Eigen::VectorXd& grad);
 
};
