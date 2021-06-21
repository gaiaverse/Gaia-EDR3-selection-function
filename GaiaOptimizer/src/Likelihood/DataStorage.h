#pragma once
#include "../Main/GlobalVariables.h"
#include "../Main/GlobalConstants.h"
#include "../DataHandling/Star.h"
#include <string>
#include <vector>

#include "../libs/JSL/JSL.h"

enum Probability { NormalApproximation, PoissonBinomial};

struct VariancePopulation
{
	double Fraction;
	std::vector<double> PowerContributions;
	
	VariancePopulation(){};
	VariancePopulation(double fraction,std::vector<double> contributions)
	{
		Fraction = fraction;
		PowerContributions = contributions;
		
		
	}; 
	double Variance(double scaling)
	{
		double v = 0;
		for (int i =0; i < PowerContributions.size(); ++i)
		{
			v += PowerContributions[i] * pow(scaling,i);
		}
		return v;
	}
	double Gradient(double scaling,bool active)
	{
		double value = 0;
		if (active)
		{
			for (int i = 1; i < PowerContributions.size(); ++i)
			{
				value+= i * PowerContributions[i] * pow(scaling,i-1);
			}
		}
		return value;
	}
};
class LikelihoodData
{
	public:
	
		LikelihoodData(const std::vector<std::vector<Star>> &data, int id);
		
		int ID;
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
		std::vector<double> p_dfdp;
		std::vector<double> hypergradient;
		std::vector<double> pt;
		std::vector<double> pml;
		std::vector<double> p;

		std::vector<double> grad_elu_xml1;
		std::vector<double> grad_elu_xml2;
		
		
		double ExpectedActiveVisitations;
		std::vector<VariancePopulation> VariancePopulations;
		Probability Mode;
		
		void GeneratePopulations(const std::vector<double> & x);
};

