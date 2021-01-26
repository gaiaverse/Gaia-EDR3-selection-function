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
		long double Value;
		std::vector<double> Gradient;
		const std::vector<Star> &Data;
		
		int MinVisits;
		Liklihood(const std::vector<Star> & data, int nPoints, int id);
		
		void Calculate(Eigen::VectorXd& position);
	private:
		std::vector<double> pmf
		std::vector<double> subpmf
		
		void Prior();
		
		void PerStarContribution(int id,Eigen::VectorXd& position);
		//new member functions go here:
};

