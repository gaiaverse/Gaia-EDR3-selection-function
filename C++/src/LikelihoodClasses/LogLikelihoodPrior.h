#pragma once
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#include "LogLikelihood.h"

using Eigen::VectorXd;
using namespace Eigen;
class LogLikelihoodPrior : public LogLikelihood
{
	public:
		LogLikelihoodPrior(const std::vector<std::vector<Star>> & data,int id): LogLikelihood(data, id){};
	
			     
	    //additional functions
	    void PriorCalculate(const std::vector<double> & x, int effectiveBatchID, int effectiveBatches);
		

};
