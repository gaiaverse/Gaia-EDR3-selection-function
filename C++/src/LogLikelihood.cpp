#include "LogLikelihood.h"

LogLikelihood::LogLikelihood(const std::vector<Star> &data, std::vector<int> & magBins, int dimension, int id): Data(data)
{
	ID = id;
	Value = 0.0;
	Gradient = Eigen::VectorXd::Zero(dimension);
	MagBins = magBins;
	
	
	
	int suitablyLargeNumber = 1024; // a number at least as large as the maximum number of entries in a single star's time series
	pmf = std::vector<double>(suitablyLargeNumber,0.0);
	subpmf = std::vector<double>(suitablyLargeNumber,0.0);
	
	//initialise an array of nbins per nt to store the modified values of the position vector each time step
	
	perBinP = std::vector<std::vector<double>>(magBins.size(),std::vector<double>(Nt,0.0));
}

void LogLikelihood::Calculate(Eigen::VectorXd& x)
{
	std::cout << "\t\t\tThe log-likelihood function hass been called" << std::endl;
	Reset();	
	

	
	GeneratePs(x);
	
	
	

	for (int i = 0; i < Data.size(); ++i)
	{
		if (ID == 0)
		{
			std::cout << "\t\tCalculating contribution from star " << i << std::endl;
		}
		PerStarContribution(i);
	}

	

}

void LogLikelihood::Reset()
{
	Value = 0;
	for (int i = 0; i < Gradient.size(); ++i)
	{
		Gradient[i] = 0;
	}
}

void LogLikelihood::GeneratePs(Eigen::VectorXd&x)
{
	for (int i = 0; i < MagBins.size(); ++i)
	{
		int bin = MagBins[i];
		int offset = Nh + Ng + bin* Nt;
		
		for (int j = 0; j < Nt; ++j)
		{
			perBinP[i][j] = 1.0/(1.0 - exp(-x[offset + j]));
		}
	}
}

void LogLikelihood::PerStarContribution(int star)
{
	Star candidate = Data[star];

	int k = candidate.nMeasure;
	int n = candidate.nVisit;
	
	//copies in-place into pmf
	direct_convolution_local(perBinP[candidate.gBin],candidate.TimeSeries,n,pmf);

	double likelihood = pmf[k];
	
	double correction = 1.0;
	for (int i = 0; i < PipelineMinVisits; ++i)
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
	
	int offset = Nh + Ng + candidate.gBin * Nt;
	for (int i = 0; i < n; ++i)
	{
		double p = perBinP[candidate.gBin][candidate.TimeSeries[i]];
		double inv_p = 1.0/(1 - p);
		
		subpmf[0] = pmf[0] * inv_p;
		for (int j = 1; j < n; ++j)
		{
			subpmf[j] = (pmf[j] - subpmf[j-1]*p)*inv_p;
		}
		subpmf[n-1] = pmf[n]/p;
		
		double dFdP = (gradient_first_term*subpmf[k-1]-gradient_second_term*subpmf[k])/likelihood - subpmf[PipelineMinVisits-1]/correction;
		double dPdX = p * (1 - p);
		int t= candidate.TimeSeries[i];

		Gradient[offset + t] = dFdP * dPdX;
	}
}



            
            
void direct_convolution_local(std::vector<double> & probsFull,std::vector<unsigned int> &probsIndex, int probslen, std::vector<double> & result)
{
	//stolen from https://github.com/biscarri1/convpoibin/blob/master/src/convpoibin.c

	int oldlen = 2; // length of old kernel
	double signal[2];
	double t,tmp;
	
	
	int index0 = probsIndex[0];
	// initialize (old kernel)
	result[0] = 1-probsFull[index0];
	result[1] = probsFull[index0];
	
	// loop through all other probs
	for(int i=1; i < probslen; ++i)
	{
		int index = probsIndex[i];
		//~ // set signal
		signal[0] = probsFull[index];
		signal[1] = 1-probsFull[index];
		
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

