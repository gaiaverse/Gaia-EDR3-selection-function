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
		std::vector<double> Gradient;
		
		LogLikelihood(const std::vector<std::vector<Star>> & data, int id);
		void Calculate(const std::vector<double> & position, int batchID, int effectiveBatches);
		
		
	protected:
		
		//member data 
		LikelihoodData Data;

		std::vector<double> TransformedPosition;
		std::vector<double> TransformedGradient;
		Eigen::Matrix<double, Nm, Nm> CholeskyKg;     

		void MakeCovarianceMatrix();
		int needletN;
		std::vector<int> needlet_u;
    	std::vector<int> needlet_v;
    	std::vector<double> needlet_w;
    	std::vector<double> bVector;


		//internal functions
		void Reset();
		void PerStarContribution(int batchId, int starID,const std::vector<double> & position);
		
		void GeneratePs(const Star * candidate,const std::vector<double> & position);
		void GenerateContribution(const Star * candidate);
		void GenerateExactContribution(const Star * candidate);
		void AssignGradients(const Star * candidate);
		
		void ForwardTransform(const std::vector<double> &z);
		void BackwardTransform();
};
