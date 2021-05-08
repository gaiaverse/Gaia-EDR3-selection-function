#pragma once
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#include "LogLikelihood.h"

using Eigen::VectorXd;
using namespace Eigen;
class LogLikelihoodPrior : public LogLikelihood
{
	public:
		LogLikelihoodPrior(const std::vector<Star> & data, const std::vector<int> &offsets, int id): LogLikelihood(data, offsets, id){};
	
		
	     
	    //additional functions
	    void Prior(const Eigen::VectorXd& RawParams, double * currentValue, Eigen::VectorXd * currentGradient);
		void MakeCovarianceMatrix();
		
		//saved values 
		bool Kg_decomposed = false;
	    Eigen::Matrix<double, Nm, Nm> CholeskyKg;     
};
