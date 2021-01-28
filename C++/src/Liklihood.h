#pragma once
#include <vector>
#include "libs/Eigen/Core"
#define EIGEN_MPL2_ONLY

#include "customClasses.h"
using Eigen::VectorXd;
//Liklihood class acts as a container for the values of the log liklihood and its gradient
//Also contains the data necessary to update these values when Calculate(newPosition) is called
class Liklihood
{
	public:
		int ID;
        const int Nh = 5; // number of hyper-hyper-parameters
        const int Ng = 35; // number of magnitude bins
        const int Nt = 8967691; // number of time bins
        const Eigen::Vector<double, Ng> magnitudes << 2.5, 10.5, 16.5, 17.5, 18.1, 18.3, 18.5, 18.7, 18.9, 19.05, 19.15, 19.25, 19.35, 19.45, 19.55, 19.65, 19.75, 19.85, 19.95, 20.05, 20.15, 20.25, 20.35, 20.45, 20.55, 20.65, 20.75, 20.85, 20.95, 21.05, 21.15, 21.25, 21.35, 21.45, 23.25; // magniitude bins
		long double Value;
		std::vector<double> Gradient;
		const std::vector<Star> &Data;
		
		
		Liklihood(const std::vector<Star> & data, int nPoints, int id);
		
		void Calculate(Eigen::VectorXd& position);
	private:
		//std::vector<double> pmf
		//std::vector<double> subpmf
		
		void Prior();
        void PriorLengthscale(double& lengthscale, int& param_index);
        void PriorVariance(double& variance, int& param_index);
        void PriorMu(Eigen::VectorXd& mu, double& m, double& tau2);
        void PriorX(Eigen::VectorXd& x, Eigen::VectorXd& mu, double& lt, double& lm, double& sigma2);
		//new member functions go here:
};

