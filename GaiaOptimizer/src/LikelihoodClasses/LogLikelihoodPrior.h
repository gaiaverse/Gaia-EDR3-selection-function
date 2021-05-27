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
	    void Prior(const Eigen::VectorXd& RawParams, double * currentValue, std::vector<double> * currentGradient, int effectiveBatches);
		void MakeCovarianceMatrix();
		
		//saved values 
		bool Kg_decomposed = false;
	    Eigen::Matrix<double, Nm, Nm> CholeskyKg;

	    double cholesky_tol = 1e-4;
	    int choleskyN;
	    std::vector<int> cholesky_u;
    	std::vector<int> cholesky_v;
    	std::vector<double> cholesky_w; 
};
