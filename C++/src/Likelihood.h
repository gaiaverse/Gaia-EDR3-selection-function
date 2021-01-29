#pragma once
#include <vector>
#include <iostream>

#define EIGEN_STACK_ALLOCATION_LIMIT 0 
#include "libs/Eigen/Core"
#include "libs/Eigen/Householder"
#include "libs/Eigen/QR"
#include <algorithm>

#define EIGEN_MPL2_ONLY

#include "Star.h"
#include "GlobalVariables.h"
using Eigen::VectorXd;
using namespace Eigen;
//Likelihood class acts as a container for the values of the log liklihood and its gradient
//Also contains the data necessary to update these values when Calculate(newPosition) is called
class Likelihood
{
	public:
		
		long double Value;
		Eigen::VectorXd Gradient;
		
		
		
		Likelihood(const std::vector<Star> & data, std::vector<int> & magBins, int dimensionality, int id);
		
		void Calculate(Eigen::VectorXd& position);
	private:
		std::vector<double> pmf;
		std::vector<double> subpmf;
		
	
		void Reset();
		void PerStarContribution(int id);
		//new member functions go here:
		
		//hardcoded parameters
		int MinVisits;
		

		
		int ID;
		const std::vector<Star> &Data;
		std::vector<int> MagBins;
		std::vector<std::vector<double>> perBinP;
		void GeneratePs(Eigen::VectorXd & position);
		
		
		void Prior(Eigen::VectorXd& params);
        void PriorLengthscale(double lengthscale, int param_index);
        void PriorVariance(double variance, int param_index);
        void PriorMu(Eigen::VectorXd& mu, double m, double tau2);
        void PriorX(Eigen::VectorXd& x, Eigen::VectorXd& mu, double lt, double lm, double sigma2);
        
};

