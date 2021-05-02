#pragma once
#include "../GlobalVariables.h"
#include "../DescentClasses/Star.h"
#include "../GenericFunctions/FileHandler.h"
#include <string>
#include <vector>
class LikelihoodData
{
	public:
	
		LikelihoodData(const std::vector<Star> &data, int id);
		
		int ID;
		
		const std::vector<Star> &Stars;
		int NStars;
		
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
		std::vector<double> dfdp;
		std::vector<double> pt;
		std::vector<double> pml;
		std::vector<double> p;
	
	
};
