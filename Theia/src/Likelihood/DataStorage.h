#pragma once
#include "../Main/GlobalVariables.h"
#include "../Main/GlobalConstants.h"
#include "../DataHandling/Star.h"
#include "MiscFunctions.h"
#include <string>
#include <vector>
#include "VariancePopulation.h"
#include "../libs/JSL/JSL.h"
#include "../Optimizer/EfficiencyVector.h"

enum Probability { NormalApproximation, PoissonBinomial};


/*!
	A class which packages most of the gnarly bits of LogLikelihood. This is essentially just a container class for the requisite data 
*/
class LikelihoodData
{
	public:
		
		//! Constructor function \param data A vector of Star objects arranged according to the minibatching schedule. 
		LikelihoodData(const std::vector<std::vector<Star>> &data);
		
		//! The storage location for the reference to the large amount of data passed to this object
		const std::vector<std::vector<Star>> &Stars;
		
		//! A vector containing the healpix id that FOV_1 is looking at at each time index
		std::vector<int> healpix_fov_1;
		
		//! A vector containing the healpix id that FOV_2 is looking at at each time index
    	std::vector<int> healpix_fov_2;
    	
    	//! A vector mapping the #Nt coarse-grain times to the #TotalScanningTime bins that FOV_i uses and that the Star objects refer to
    	std::vector<int> time_mapping;

	
	
		//! A vector of length #NumberLargerThanMaxObservations, used to store intermediary results for the convolutions in the PoissonBinomial calculations
		std::vector<std::vector<double>> pmf_forward;
		
		//! A vector of length #NumberLargerThanMaxObservations, used to store intermediary results for the convolutions in the PoissonBinomial calculations
		std::vector<std::vector<double>> pmf_backward;
		
		//! A vector of length #NumberLargerThanMaxObservations, used to store intermediary results for the convolutions in the PoissonBinomial calculations
		std::vector<std::vector<double>> subpmf;
		
		//!A vector of length #NumberLargerThanMaxObservations, used to store the value of the derivative of the per-star likelihood with respect to the individual observation probabilities
		std::vector<double> dfdp_constantN;
		
		//! A vector of length #NumberLargerThanMaxObservations, used to store the value of the derivative of the per-star likelihood with respect to the variance-scaling parameter
		double dfdN_constantP;
		
		//! The contributions towards the hyperparameter gradients
		std::vector<double> hypergradient;
		
		//! A vector of length #NumberLargerThanMaxObservations, used to store the values of the proposed operating efficiency of Gaia at each of the times that the star was observed
		std::vector<double> pt;
		
		//! A vector of length #NumberLargerThanMaxObservations, used to store the values of the proposed operating efficiency of Gaia at each of the spatial points that the star was observed at
		std::vector<double> pml;
		
		//! A vector of length #NumberLargerThanMaxObservations, used to store the pointwise product of pt * pml --> the total likelihood that the star was seen at each of its given observations
		std::vector<double> p;
		
		//! A vector of length #NumberLargerThanMaxObservations, used to hold the spatial gradient components for FOV1
		std::vector<double> grad_elu_xml1;
		
		//! A vector of length #NumberLargerThanMaxObservations, used to hold the spatial gradient components for FOV2
		std::vector<double> grad_elu_xml2;
		
		
		//! A container populated by the GeneratePopulations() function
		std::vector<VariancePopulation> VariancePopulations;
		
		//! An enum which determines if the probability model used is PoissonBinomial, or the NormalApproximation
		Probability Mode;

		/*!
		 * Generates a new set of VariancePopulation objects to be re-used throughout the LogLikelihood::Calculate() loop, and stores them in #VariancePopulations.
		*/
		void GeneratePopulations(const EfficiencyVector & x);
};

