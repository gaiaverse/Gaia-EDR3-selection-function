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
	double BaselineVariance;
	double LinearVariance;
	double QuadraticVariance;
	
	//~ VariancePopulation(){};
	VariancePopulation(double fraction,double base,double linearscaling, double quadraticscaling)
	{
		Fraction = fraction;
		BaselineVariance = base;
		LinearVariance = linearscaling;
		QuadraticVariance = quadraticscaling;
		Print();
	}; 
	void Print()
	{
				
		std::cout << "I have: " << std::endl;
		std::cout << Fraction << "   " << BaselineVariance << "   " << LinearVariance << "   " << QuadraticVariance << std::endl;
	}
};
class LikelihoodData
{
	public:
	
		LikelihoodData(const std::vector<std::vector<Star>> &data, int id);
		
		int ID;
		
		const std::vector<std::vector<Star>> &Stars;
		//~ int NStars;
		
		//Indexing data, allows us to index properly into the 
		//spatial and temporal parts of the vector
		std::vector<int> healpix_fov_1;
    	std::vector<int> healpix_fov_2;
    	std::vector<int> time_mapping;
		//~ std::vector<int> magtime_mapping;

	
	
		//frequently overwritten vectors
		//makes sense to initialise these vectors to a large size and be careful about memory access
		//rather than continually create and delete them 
		std::vector<std::vector<double>> pmf_forward;
		std::vector<std::vector<double>> pmf_backward;
		std::vector<std::vector<double>> subpmf;
		std::vector<double> dfdp;
		std::vector<double> pt;
		//~ std::vector<double> ptm;
		std::vector<double> pml;
		std::vector<double> p;

		std::vector<double> grad_elu_xml1;
		std::vector<double> grad_elu_xml2;
		
		std::vector<VariancePopulation> VariancePopulations;
		Probability Mode;
};

