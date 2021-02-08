#pragma once
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#include "Likelihood.h"

using Eigen::VectorXd;
using namespace Eigen;
class LikelihoodPrior : public Likelihood
{
	public:
		LikelihoodPrior(const std::vector<Star> & data, std::vector<int> & magBins, int dimensionality, int id): Likelihood(data,magBins,dimensionality, id){};
	
		void Calculate(Eigen::VectorXd &position);
	private:
		
		//values saved for reuse in calculations
		bool Kg_decomposed = false;
	    Eigen::Matrix<double, Ng, Ng> Kg;
	    Eigen::Matrix<double, Ng, Ng> invKg;
	    double logdetKg;
	
	
		void Prior(Eigen::VectorXd& params);
        void PriorLengthscale(double lengthscale, int param_index);
        void PriorVariance(double variance, int param_index);
        void PriorMu( Map<VectorXd>& mu);
        void PriorX( Map<VectorXd>& x, Map<VectorXd>& mu, double lt, double sigma2);
};
