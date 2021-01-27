#include "Liklihood.h"

void direct_convolution_local(std::vector<double> probs, int probslen, std::vector<double> & result)
{
	//stolen from https://github.com/biscarri1/convpoibin/blob/master/src/convpoibin.c

	int oldlen = 2; // length of old kernel
	double signal[2];
	double t,tmp;
	
	
	// initialize (old kernel)
	result[0] = 1-probs[0];
	result[1] = probs[0];
	
	// loop through all other probs
	for(int i=1; i < probslen; ++i)
	{

		//~ // set signal
		signal[0] = probs[i];
		signal[1] = 1-probs[i];
		
		// initialize result and calculate the two edge cases
		result[oldlen] = signal[0] * result[oldlen-1];
		
		t = result[0];
		result[0] = signal[1]*t;
		  
		//calculate the interior cases
		for(int j=1; j < oldlen; ++j)
		{
			tmp=result[j];
			result[j] = signal[0] * t + signal[1] * result[j];
			t=tmp;
		}
		  
		oldlen++;
	}
}

Liklihood::Liklihood(const std::vector<Star> &data, const std::vector<int> & magBins, int dimension, int id): Data(data)
{
	ID = id;
	Value = 0.0;
	Gradient = std::vector<double>(dimension,0.0);
	MagBins = magBins;
	
	MinVisits = 5; //hard-coded parameter, the number of times a star has to be detected for it to enter gaia pipeline
	Nt = 9e6;
	Ng = 35;
	Nh = 5;
	
	int suitablyLargeNumber = 1024; // a number at least as large as the maximum number of entries in a single star's time series
	pmf = std::vector<double>(suitablyLargeNumber,0.0);
	subpmf = std::vector<double>(suitablyLargeNumber,0.0);
}

void Liklihood::Calculate(Eigen::VectorXd& x)
{

	Reset();	
	
	
	std::vector<double> probs = {0.3,0.6,0.2,0.1,0.9,0.5,0.3,0.2,0.8};
	if (ID == 0)
	{
		Prior();
	}

	for (int i = 0; i < Data.size(); ++i)
	{
		PerStarContribution(i,probs);
	}


}

void Liklihood::PerStarContribution(int star, std::vector<double> & probs)
{
	Star candidate = Data[star];

	int k = candidate.nMeasure;
	int n = candidate.nVisit;
	
	//copies in-place into pmf
	direct_convolution_local(probs,n,pmf);

	double likelihood = pmf[k];
	
	double correction = 1.0;
	for (int i = 0; i < MinVisits; ++i)
	{
		correction -= pmf[i];
	}
	
	Value += log(likelihood / correction);
	
	double gradient_first_term = 1.0;
	double gradient_second_term = 1.0;
	if (k == 0)
	{
		gradient_first_term = 0.0;
	}
	else
	{
		if (k == n)
		{
			gradient_second_term = 0.0;
		}
	}
	
	int offset = 0;//Nh + Ng + candidate.gBin * Nt;
	for (int i = 0; i < n; ++i)
	{
		double p = probs[i];
		double inv_p = 1.0/(1 - p);
		
		subpmf[0] = pmf[0] * inv_p;
		for (int j = 1; j < n; ++j)
		{
			subpmf[j] = (pmf[j] - subpmf[j-1]*p)*inv_p;
		}
		subpmf[n-1] = pmf[n]/p;
		
		double gradi = (gradient_first_term*subpmf[k-1]-gradient_second_term*subpmf[k])/likelihood - subpmf[MinVisits-1]/correction;
		int t= candidate.TimeSeries[i];

		Gradient[offset + t] = gradi;
	}
}

void Liklihood::Prior()
{
	//Value += whatever
	//Gradient[i] += etc
}

void Liklihood::Reset()
{
	Value = 0;
	for (int i = 0; i < Gradient.size(); ++i)
	{
		Gradient[i] = 0;
	}
}
