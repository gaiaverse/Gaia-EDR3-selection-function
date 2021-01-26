#include "Liklihood.h"

double direct_convolution_local(double probs[], int probslen, double result[])
{
	//stolen from https://github.com/biscarri1/convpoibin/blob/master/src/convpoibin.c
  
  int i,j;
  int oldlen = 2; // length of old kernel
  double signal[2];
  double t,tmp;


  // initialize (old kernel)
  result[0] = 1-probs[0];
  result[1] = probs[0];

  // loop through all other probs
  for(i=1; i < probslen; i++){

    // set signal
    signal[0] = probs[i];
    signal[1] = 1-probs[i];

    // initialize result and calculate the two edge cases
    result[oldlen] = signal[0] * result[oldlen-1];

    t = result[0];
    result[0] = signal[1]*t;
      
    //calculate the interior cases
    for(j=1; j < oldlen; j++){
      tmp=result[j];
      result[j] = signal[0] * t + signal[1] * result[j];
      t=tmp;
    }
      
    oldlen++;
  }

}


Liklihood::Liklihood(const std::vector<Star> &data, int nPoints, int id): Data(data)
{
	ID = id;
	Value = 0.0;
	Gradient = std::vector<double>(nPoints,0.0);
	MinVists = 5; //hard-coded parameter, the number of times a star has to be detected for it to enter gaia pipeline
	
	int suitablyLargeNumber = 1024; // a number at least as large as the maximum number of entries in a single star's time series
	pmf = std::vector<double(suitablyLargeNumber,0.0);
	subpmf = std::vector<double(suitablyLargeNumber,0.0);
}


void Liklihood::Calculate(Eigen::VectorXd& x)
{

	Value = 0;
	Gradient[0] = 0;
	Gradient[1] = 0;
	Gradient[2] = 0;
	Gradient[3] = 0;	
	
	
	
	if (ID == 0)
	{
		Prior();
	}

	for (int i = 0; i < Data.size(); ++i)
	{
		PerStarContribution(i,x);
	}

	
}


void Liklihood::PerStarContribution(int star, Eigen::VectorXd& probs)
{
	Star candidate = Data[i];
	
	int k = candidate.nMeasure;
	int n = candidate.nVisit;
	
	//copies in-place into pmf
	direct_convolution_local(probs,n,pmf);
	double liklihood = pmf[k];
	
	double correction = 1.0;
	for (int i = 0; i < MinVisits; ++i)
	{
		correction -= pmf[i];
	}
	
	Value += log(likilhood / correction)
	
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
		
		double gradi = (gradient_first_term*subpmf[k-1]-gradient_second_term*subpmf[k])/likelihood - subpmf[c-1]/correction;
		int t= TimeSeries[i];
		
		int preguff = Nh + Ng + candidate.gBin * Nt;
		Gradient[preguff + t] = gradi;
	}
}

void Liklihood::Prior()
{
	//Value += whatever
	//Gradient[i] += etc
}
