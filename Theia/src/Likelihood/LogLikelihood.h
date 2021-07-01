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


#include "../DataHandling/Star.h"
#include "../Main/GlobalVariables.h"
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
		std::vector<double> Gradient;
		
		LogLikelihood(const std::vector<std::vector<Star>> & data, int id);
		void Calculate(const std::vector<double> & position, int batchID, int effectiveBatches, int maxBatches);
		LikelihoodData Data;
	protected:
		
		//member data 
		

		//internal functions
		void Reset();
		void PerStarContribution(int batchId, int starID,const std::vector<double> & position);
		
		void GeneratePs(const Star * candidate,const std::vector<double> & position);
		void GenerateContribution(const Star * candidate);
		
		void AssignGradients(const Star * candidate);
		
		void NormalContribution(const Star * candidate);
		void PoissonContribution(const Star * candidate);
		void ExactPoissonContribution(const Star * candidate);
};
