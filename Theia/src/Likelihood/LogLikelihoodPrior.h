#pragma once
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#include "LogLikelihood.h"

using Eigen::VectorXd;
using namespace Eigen;
class LogLikelihoodPrior : public LogLikelihood
{
	public:
		LogLikelihoodPrior(const std::vector<std::vector<Star>> & data): LogLikelihood(data){};
	
		
	     
	    //additional functions
	    void TransformPrior(const std::vector<double> & TransformPosition, double * currentValue, std::vector<double> & TransformGradient, int effectiveBatches);
	    void RawPrior(const std::vector<double>& RawParams, double * currentValue, std::vector<double> * currentGradient, int effectiveBatches);
		

};
