#pragma once
#include "../Main/GlobalVariables.h"
#include "../Main/GlobalConstants.h"
#include "../DataHandling/Star.h"
#include "MiscFunctions.h"
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
		double v = pow(PowerContributions[0],2);
		
		for (int i = 1; i <= hyperOrder/2; ++i)
		{
			int power = 2*i;
			double term = PowerContributions[power-1] + PowerContributions[power] * scaling;
			v += pow(term,power); 
		} 
		

		return v;
		
	}
	double dVariancedN(double scaling)
	{
		double v = 0;
		
		for (int i = 1; i <= hyperOrder/2; ++i)
		{
			int power = 2*i;
			double term = PowerContributions[power-1] + PowerContributions[power] * scaling;
			v += power * pow(term,power-1) * PowerContributions[power]; 
		} 
		return v;
	}
	double dVariancedAlpha(int term, double scaling)
	{
		double v;
		if (term > 0)
		{
			if (term % 2 == 0)
			{
				double bracket = PowerContributions[term-1] + PowerContributions[term] * scaling;
				v = term * pow(bracket,term - 1) * scaling;
			}
			else
			{
				double bracket = PowerContributions[term] + PowerContributions[term+1] * scaling;
				v = (term + 1 ) * pow(bracket,term);
			}
		}
		else
		{
			v = 2 * PowerContributions[0];
		}
		return v;
	}
	void Print(double scale)
	{
		std::cout << "I have fraction: " << Fraction << " and\n" ;
		for (int i = 0; i <= hyperOrder; ++i)
		{
			std::cout << "a" << i << " = " << PowerContributions[i] << ";" <<std::endl;
		}
		
		std::cout << "SO: V(n) = " << PowerContributions[0];
		for (int i = 1; i<= hyperOrder/2; ++i)
		{
			int power = 2*i-1;
			std::cout << "+ (" << PowerContributions[power] << " + " << PowerContributions[power+1] << "*n).^" << power+1;
		}
		std::cout<< "\n With n = " << scale << " I have: \n";
		
		
		std::cout << "V = " << Variance(scale) << std::endl;
		std::cout << "dVdN = " << dVariancedN(scale) << std::endl;
		for (int i = 0; i < hyperOrder + 1; ++i)
		{
			std::cout << "dVdA_" << i << " = " << dVariancedAlpha(i,scale) << std::endl;
		}
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

