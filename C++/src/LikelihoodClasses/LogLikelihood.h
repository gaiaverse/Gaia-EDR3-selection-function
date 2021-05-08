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
#include "../DataOperators/Star.h"
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
		int StarsUsed;
		Eigen::VectorXd Gradient;
		
		LogLikelihood(const std::vector<Star> & data, const std::vector<int> & offsets, int id);
		void Calculate(Eigen::VectorXd& position, int batchID, int effectiveBatches);
		
	protected:
		
		//member data 
		LikelihoodData Data;

		//internal functions
		void Reset();
		void PerStarContribution(int id,Eigen::VectorXd & position);
		
		void GeneratePs(const Star * candidate,Eigen::VectorXd & position);
		void GenerateContribution(const Star * candidate);
		void GenerateExactContribution(const Star * candidate);
		void AssignGradients(const Star * candidate);
};
