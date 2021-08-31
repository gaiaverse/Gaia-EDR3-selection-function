#pragma once
#include <vector>


/*!
 * This structure uses the hyperparameters from the EfficiencyVector to generate the variance + associated derivatives.
*/
struct VariancePopulation
{
	//!The fraction of stars which belong to this population
	double Fraction;
	
	//!A list of the (unconstrained) prefactors within the polynomial expansion
	std::vector<double> Term;
	
	//!Default constructor
	VariancePopulation(){};
	
	/*!
	 * Constructor class
	 * \param fraction The value assigned to Fraction
	 * \param contributions The vector assigned to Term
	 */
	VariancePopulation(double fraction,std::vector<double> contributions)
	{
		Fraction = fraction;
		Term = contributions;				
	}; 
	
	/*!
	 *Calculates the Variance contribution from this population
	 * \param n The variance scaling parameter
	 * \returns v a guaranteed positive number 
	*/
	double Variance(double n)
	{
		double v = pow(Term[0],2);
		
		for (int i = 1; i <= hyperOrder/2; ++i)
		{
			int power = 2*i;
			double term = Term[power-1] + Term[power] * n;
			v += pow(term,power); 
		} 
		

		return v;
		
	}
	
	/*!
	 *Calculates the derivative of the Variance with respect to the scaling parameter
	 * \param n The variance scaling parameter
	 * \return The analytical derivative of the variance
	*/
	double dVariancedN(double n)
	{
		double v = 0;
		
		for (int i = 1; i <= hyperOrder/2; ++i)
		{
			int power = 2*i;
			double term = Term[power-1] + Term[power] * n;
			v += power * pow(term,power-1) * Term[power]; 
		} 
		return v;
	}
	
	/*!
	 *Calculates the derivative of the variance with respect to the ith-order Term parameter
	 * \param term the order of the parameter with which the derivative is taken
	 * \param n The variance scaling parameter
	 * \returns dv/d alpha_term
	*/
	double dVariancedAlpha(int term, double n)
	{
		double v;
		if (term > 0)
		{
			if (term % 2 == 0)
			{
				double bracket = Term[term-1] + Term[term] * n;
				v = term * pow(bracket,term - 1) * n;
			}
			else
			{
				double bracket = Term[term] + Term[term+1] * n;
				v = (term + 1 ) * pow(bracket,term);
			}
		}
		else
		{
			v = 2 * Term[0];
		}
		return v;
	}
};

