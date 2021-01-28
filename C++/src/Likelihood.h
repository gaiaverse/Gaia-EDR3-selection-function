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
		
		void Prior();
		void Reset();
		void PerStarContribution(int id);
		//new member functions go here:
		
		//hardcoded parameters
		int MinVisits;
		int Nh;
		int Ng;
		int Nt;
		
		int ID;
		const std::vector<Star> &Data;
		std::vector<int> MagBins;
		std::vector<std::vector<double>> perBinP;
		void GeneratePs(Eigen::VectorXd & position);
};

