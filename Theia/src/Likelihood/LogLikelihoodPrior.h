#pragma once
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#include "LogLikelihood.h"
#include "../Main/EfficiencyVector.h"
using Eigen::VectorXd;
using namespace Eigen;
class LogLikelihoodPrior : public LogLikelihood
{
	public:
		LogLikelihoodPrior(const std::vector<std::vector<Star>> & data): LogLikelihood(data){};
	
		
	     
	    //additional functions
	    double TransformPrior(EfficiencyVector & x, int effectiveBatches);
	    double RawPrior(EfficiencyVector & x, int effectiveBatches);
		

};
