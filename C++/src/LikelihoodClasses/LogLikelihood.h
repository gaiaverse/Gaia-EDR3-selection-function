#pragma once
#include <vector>
#include <iostream>
#include <string>
#include <iomanip>
#define EIGEN_STACK_ALLOCATION_LIMIT 0 
#include "../libs/Eigen/Core"
#include "../libs/Eigen/Householder"
#include "../libs/Eigen/QR"
#include <algorithm>

#define EIGEN_MPL2_ONLY

#include "../GenericFunctions/FileHandler.h"
#include "../DescentClasses/Star.h"
#include "../GlobalVariables.h"
#include "DataStorage.h"
#include "ProbabilityFunctions.h"
using Eigen::VectorXd;
using namespace Eigen;
//Likelihood class acts as a container for the values of the log liklihood and its gradient
//Also contains the data necessary to update these values when Calculate(newPosition) is called

class LogLikelihood
{
	public:
		double Value;
		Eigen::VectorXd Gradient;
		
		LogLikelihood(const std::vector<Star> & data, int id);
		void Calculate(Eigen::VectorXd& position);
		
	protected:
		
		//member data 
		LikelihoodData Data;

		//internal functions
		void Reset();
		void PerStarContribution(int id,Eigen::VectorXd & position);
		void inline CalculateSubPMF(int i, int n, int k,std::vector<double> & ps);
		void inline SubPMF_Forward(double p, int start, int end);
		void inline SubPMF_Backward(double p, int start, int end, int n);
		
		void GenerateDerivatives_Normal();
		void GenerateDerivatives_Long();
		
		
};
