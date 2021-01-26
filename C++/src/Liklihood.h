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
        int Nh = 5; // number of hyper-hyper-parameters
        int Ng = 35; // number of magnitude bins
        int Nt = 8967691; // number of time bins
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
		//new member functions go here:
};

