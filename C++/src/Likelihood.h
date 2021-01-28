#pragma once
#include <vector>
#include <iostream>
#include "libs/Eigen/Core"
#define EIGEN_MPL2_ONLY

#include "Star.h"
using Eigen::VectorXd;
//Likelihood class acts as a container for the values of the log liklihood and its gradient
//Also contains the data necessary to update these values when Calculate(newPosition) is called
class Likelihood
{
	public:
		
		long double Value;
		std::vector<double> Gradient;
		
		
		
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
		const int Nh = 5; // number of hyper-hyper-parameters
        const int Ng = 35; // number of magnitude bins
        const int Nt = 8967691; // number of time bins
        const Eigen::Vector<double, Ng> magnitudes << 2.5, 10.5, 16.5, 17.5, 18.1, 18.3, 18.5, 18.7, 18.9, 19.05, 19.15, 19.25, 19.35, 19.45, 19.55, 19.65, 19.75, 19.85, 19.95, 20.05, 20.15, 20.25, 20.35, 20.45, 20.55, 20.65, 20.75, 20.85, 20.95, 21.05, 21.15, 21.25, 21.35, 21.45, 23.25; // magniitude bins
		long double Value;
		
		int ID;
		const std::vector<Star> &Data;
		std::vector<int> MagBins;
		std::vector<std::vector<double>> perBinP;
		void GeneratePs(Eigen::VectorXd & position);
		
		
		void Prior(Eigen::VectorXd& params);
        void PriorLengthscale(double& lengthscale, int& param_index);
        void PriorVariance(double& variance, int& param_index);
        void PriorMu(Eigen::VectorXd& mu, double& m, double& tau2);
        void PriorX(Eigen::VectorXd& x, Eigen::VectorXd& mu, double& lt, double& lm, double& sigma2);
};

