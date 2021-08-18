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
/*!
 This class holds the data and functions necessary to calculate log p(data | x), where x is the provided EfficiencyVector. The majority of the guts + internal workings of the class are hidden away in the LikelihoodData.
*/
class LogLikelihood
{
	public:
		
		//! The last-calculated value of the Calculate() function
		double Value; 
		//! The last-calculated gradient of the Calculate() function with respect to the associated `position`
		std::vector<double> Gradient;
		//! The number of stars within the last-called minibatch
		int StarsUsed;
		
		
		
		
		/*!
		 \brief Constructor function. 
		 \param data A vector of Star objects arranged according to the minibatching schedule. 
		 \param id the MPI ID of the running process. 
		*/
		LogLikelihood(const std::vector<std::vector<Star>> & data, int id);
		
		/*!The key function: executes a single minibatch calculation of the loglikelihood. 
		 \param position The current EfficiencyVector
		 \param batchID The minibatch to be executed
		 \param effectiveBatches The current number of active minibatches
		 \param maxBatches The original number of active minibatches
		 \returns Assigns the value of L to #Value, and the associated gradient to #Gradient
		*/
		void Calculate(const std::vector<double> & position, int batchID, int effectiveBatches, int maxBatches);
		
	protected:
		
		//! member data 
		LikelihoodData Data;

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
