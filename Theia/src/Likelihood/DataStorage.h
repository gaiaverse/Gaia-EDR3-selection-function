#pragma once
#include "../Main/GlobalVariables.h"
#include "../Main/GlobalConstants.h"
#include "../DataHandling/Star.h"
#include "MiscFunctions.h"
#include <string>
#include <vector>
#include "VariancePopulation.h"
#include "../libs/JSL/JSL.h"

enum Probability { NormalApproximation, PoissonBinomial};


/*!
	A class which packages most of the gnarly bits of :doc:`likelihood`. This is essentially just a container class for the requisite data 
*/
class LikelihoodData
{
	public:
		
		//! Constructor function
		LikelihoodData(const std::vector<std::vector<Star>> &data);
		
		const std::vector<std::vector<Star>> &Stars;
		
		//Indexing data, allows us to index properly into the 
		//spatial and temporal parts of the vector
		std::vector<int> healpix_fov_1;
    	std::vector<int> healpix_fov_2;
    	std::vector<int> time_mapping;

	
	
		//frequently overwritten vectors
		//makes sense to initialise these vectors to a large size and be careful about memory access
		//rather than continually create and delete them 
		std::vector<std::vector<double>> pmf_forward;
		std::vector<std::vector<double>> pmf_backward;
		std::vector<std::vector<double>> subpmf;
		
		std::vector<double> dfdp_constantN;
		double dfdN_constantP;
		
		std::vector<double> hypergradient;
		std::vector<double> pt;
		std::vector<double> pml;
		std::vector<double> p;

		std::vector<double> grad_elu_xml1;
		std::vector<double> grad_elu_xml2;
		
		
		std::vector<VariancePopulation> VariancePopulations;
		Probability Mode;
		
		
		
		
		
		
		void GeneratePopulations(const std::vector<double> & x);
};

