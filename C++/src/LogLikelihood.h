#pragma once
#include <vector>
#include <iostream>
#include <string>

#define EIGEN_STACK_ALLOCATION_LIMIT 0 
#include "libs/Eigen/Core"
#include "libs/Eigen/Householder"
#include "libs/Eigen/QR"
#include <algorithm>

#define EIGEN_MPL2_ONLY

#include "FileHandler.h"
#include "Star.h"
#include "GlobalVariables.h"
using Eigen::VectorXd;
using namespace Eigen;
//Likelihood class acts as a container for the values of the log liklihood and its gradient
//Also contains the data necessary to update these values when Calculate(newPosition) is called

class LogLikelihood
{
	public:
		
		long double Value;
		Eigen::VectorXd Gradient;
		
		
		
		LogLikelihood(const std::vector<Star> & data, std::vector<int> & magBins, int dimensionality, int id);
		
		void Calculate(Eigen::VectorXd& position);
		
		int needletN;
		std::vector<int> needlet_u;
    	std::vector<int> needlet_v;
    	std::vector<double> needlet_w;
	protected:
			
		//member data 
		int ID;
		const std::vector<Star> &Data;

		//internal functions
		void Reset();
		void PerStarContribution(int id,Eigen::VectorXd & position);
	
		//compile-time structures for holding data
		int suitablyLargeNumber = 1024;
		
		std::vector<double> pmf = std::vector<double>(suitablyLargeNumber,0.0);
		std::vector<double> subpmf  = std::vector<double>(suitablyLargeNumber,0.0);

		std::string healpix_fov_file;
		std::string needlet_file;
		std::vector<int> healpix_fov_1;
    	std::vector<int> healpix_fov_2;

    	std::vector<int> time_mapping;
    
};


//function stolen from a git repo somewhere, see the implementation for more details
void direct_convolution_local(std::vector<double> & probs, int probslen, std::vector<double> & result);
