#include "ProbabilityFunctions.h"


          
void poisson_binomial_pmf_forward(std::vector<double> &  probs, int probslen, std::vector<std::vector<double>> & result)
{
	//stolen from https://github.com/biscarri1/convpoibin/blob/master/src/convpoibin.c

	int oldlen = 2; // length of old kernel
	double p,q;
	
	// initialize (old kernel)
	result[0][0] = 1.0-probs[0];
	result[0][1] = probs[0];
	
	// loop through all other probs
	for(int i=1; i < probslen; ++i)
	{
		// set signal
		p = probs[i];
		q = 1.0 - p;
		
		// initialize result and calculate the two edge cases
		result[i][0] = q * result[i-1][0];
		result[i][oldlen] = p * result[i-1][oldlen-1];
		
		//calculate the interior cases
		for(int j=1; j < oldlen; ++j)
		{
			result[i][j] = p * result[i-1][j-1] + q * result[i-1][j];
		}
		  
		oldlen++;
	}
}

void  poisson_binomial_pmf_backward(std::vector<double> &  probs, int probslen, std::vector<std::vector<double>> & result)
{
	//stolen from https://github.com/biscarri1/convpoibin/blob/master/src/convpoibin.c

	int oldlen = 2; // length of old kernel
	double p,q;
	
	// initialize (old kernel)
	result[probslen-1][0] = 1.0-probs[probslen-1];
	result[probslen-1][1] = probs[probslen-1];
	
	// loop through all other probs
	for(int i=probslen-2; i >= 0; --i)
	{

		// set signal
		p = probs[i];
		q = 1.0 - p;
		
		// initialize result and calculate the two edge cases
		result[i][0] = q * result[i+1][0];
		result[i][oldlen] = p * result[i+1][oldlen-1];
		
		//calculate the interior cases
		for(int j=1; j < oldlen; ++j)
		{
			result[i][j] = p * result[i+1][j-1] + q * result[i+1][j];
		}
		  
		oldlen++;
	}
}

void  poisson_binomial_subpmf(int m, int probslen, std::vector<std::vector<double>> & pmf_forward, std::vector<std::vector<double>> & pmf_backward, std::vector<double> & result)
{
	double conv = 0;

	result[0] = pmf_backward[1][m];
	result[probslen-1] = pmf_forward[probslen-2][m];
	for(int i = 1; i < probslen - 1; ++i)
	{
		
		conv = 0; //pmf_forward[i-1][0]*pmf_backward[i+1][m];
		int lowerBound = std::max(0,m-probslen+i+1);
		int upperBound = std::min(m,i)+1;
		
		for(int j =lowerBound; j < upperBound; ++j)
		{
			conv += pmf_forward[i-1][j]*pmf_backward[i+1][m-j];
		}
		result[i] = conv;
	}

}

double  log_add_exp(double a, double b)
{
    if (a > b)
    {
        return a + log1p(exp(b-a));
    }
    else
    {
        return b + log1p(exp(a-b));
    }
}

void  poisson_binomial_lpmf_forward(std::vector<double> & probs, int probslen, std::vector<std::vector<double>> & result)
{
	//stolen from https://github.com/biscarri1/convpoibin/blob/master/src/convpoibin.c

	int oldlen = 2; // length of old kernel
	double log_p, log_q;
	
	// initialize (old kernel)
	result[0][0] = log1p(-probs[0]);
	result[0][1] = log(probs[0]);
	
	// loop through all other probs
	for(int i=1; i < probslen; ++i)
	{
		// set signal
		log_p = log(probs[i]);
		log_q = log1p(-probs[i]);
		
		// initialize result and calculate the two edge cases
		result[i][0] = log_q + result[i-1][0];
		result[i][oldlen] = log_p + result[i-1][oldlen-1];
		
		//calculate the interior cases
		for(int j=1; j < oldlen; ++j)
		{
			result[i][j] = log_add_exp( log_p + result[i-1][j-1], log_q + result[i-1][j]);
		}
		  
		oldlen++;
	}
}

void  poisson_binomial_lpmf_backward(std::vector<double> & probs, int probslen, std::vector<std::vector<double>> & result)
{
	//stolen from https://github.com/biscarri1/convpoibin/blob/master/src/convpoibin.c

	int oldlen = 2; // length of old kernel
	double log_p, log_q;
	
	// initialize (old kernel)
	result[probslen-1][0] = log1p(-probs[probslen-1]);
	result[probslen-1][1] = log(probs[probslen-1]);
	
	// loop through all other probs
	for(int i=probslen-2; i >= 0; --i)
	{
		// set signal
		log_p = log(probs[i]);
		log_q = log1p(-probs[i]);
		
		// initialize result and calculate the two edge cases
		result[i][0] = log_q + result[i+1][0];
		result[i][oldlen] = log_p + result[i+1][oldlen-1];
		
		//calculate the interior cases
		for(int j=1; j < oldlen; ++j)
		{
			result[i][j] = log_add_exp( log_p + result[i+1][j-1], log_q + result[i+1][j]);
		}
		  
		oldlen++;
	}
}

void  poisson_binomial_sublpmf(int m, int probslen, std::vector<std::vector<double>> & lpmf_forward, std::vector<std::vector<double>> & lpmf_backward, std::vector<double> & result)
{
	double conv;

	result[0] = lpmf_backward[1][m];
	result[probslen-1] = lpmf_forward[probslen-2][m];
	for(int i = 1; i < probslen - 1; ++i)
	{
		int lowerBound = std::max(0,m-probslen+i+1);
		int upperBound = std::min(m,i)+1;
		conv = -99999999;
		for(int j =lowerBound; j < upperBound; ++j)
		{
			conv = log_add_exp(conv,lpmf_forward[i-1][j]+lpmf_backward[i+1][m-j]);
		}
		result[i] = conv;
	}

}

void logphi(double z, double& f, double& df)
{
	
	
    if (z > 6)
    {
        // zeroth case: log(1+x) approx x
        f = -0.5 * erfc(z*one_over_root2);
        df = exp(-z*z/2-f)*one_over_root2pi;
    }
    else 
    {
		if (z*z < 0.0492)
	    {
	        // first case: close to zero
	        double lp0 = -z*one_over_root2pi;
	        f = 0;
	        for(int i = 0; i < 14; ++i)
	        {
	            f = lp0*(logphi_c[i]+f);
	        }
	        f = -2*f - log(2);
	        df = exp(-z*z/2-f)*one_over_root2pi;
	    }
		else 
		{
			if (z < -11.3137)
			{
				// second case: very small
				double num = 0.5641895835477550741;
				for(int i = 0; i < 5; ++i)
				{
					num = -z*num*one_over_root2 + logphi_r[i];
				}
			
				double den = 1.0;
				for(int i = 0; i < 6; ++i)
				{
					den = -z*den*one_over_root2 + logphi_q[i];
				}
			
				f = log(num/den/2) - z*z/2;
				df = abs(den/num) * root2_over_pi;
			}
			else
		    {
		        // third case: everything else
		        f = log(erfc(-z*one_over_root2)/2);
		        df = exp(-z*z/2-f)*one_over_root2pi;
		    }
		}
	}
		
}






double poisson_binomial_normal_lpmf(int k, const std::vector<double> & probs, int probslen, std::vector<double> & gradient, std::vector<VariancePopulation> & populations)
{
	
	double m_base = 0.0;
	double s2_base = 0;

	for(int i = 0; i < probslen; ++i)
	{
        m_base += probs[i];
        s2_base += probs[i]*(1.0-probs[i]);
	}
	
	std::vector<double> populationValues(populations.size(),0.0);
	std::vector<std::vector<double>> populationGradients(populations.size(), std::vector<double>(probslen,0.0));
	for (int i =0; i < populations.size(); ++i)
	{
		double m = m_base;
		
		bool mScaling = true;
		int scaling;
		int mGradientFactor;
		if (mScaling)
		{
			mGradientFactor = 1;
			scaling = m;
		}
		else
		{
			mGradientFactor = 0;
			scaling = probslen;
		}
		
		double s2 = s2_base + populations[i].BaselineVariance + populations[i].LinearVariance * scaling + populations[i].QuadraticVariance * scaling*scaling;
	    double s = sqrt(s2);

	    

		//~ double logPhiUpper, dlogPhiUpper;
		//~ double logPhiLower, dlogPhiLower;
		//~ double logPhiDifference;
		//~ double dlpmf_dm, dlpmf_ds2;

		//~ if ( (k-m)/s < 0 ){
		    //~ logphi((k-m+0.5)/s,logPhiUpper, dlogPhiUpper); 
		    //~ logphi((k-m-0.5)/s,logPhiLower, dlogPhiLower); 
			//~ logPhiDifference = logPhiUpper + log1p(-exp(logPhiLower-logPhiUpper));
			//~ dlpmf_dm = -(exp(logPhiUpper-logPhiDifference)*dlogPhiUpper - exp(logPhiLower-logPhiDifference)*dlogPhiLower)/s;
		    //~ dlpmf_ds2 = -0.5*((k-m+0.5)*exp(logPhiUpper-logPhiDifference)*dlogPhiUpper - (k-m-0.5)*exp(logPhiLower-logPhiDifference)*dlogPhiLower)/s/s2;
		//~ }
		//~ else{
			//~ logphi(-(k-m-0.5)/s,logPhiUpper, dlogPhiUpper); 
		    //~ logphi(-(k-m+0.5)/s,logPhiLower, dlogPhiLower); 
			//~ logPhiDifference = logPhiUpper + log1p(-exp(logPhiLower-logPhiUpper));
			//~ dlpmf_dm = (exp(logPhiUpper-logPhiDifference)*dlogPhiUpper - exp(logPhiLower-logPhiDifference)*dlogPhiLower)/s;
		    //~ dlpmf_ds2 = 0.5*((k-m-0.5)*exp(logPhiUpper-logPhiDifference)*dlogPhiUpper - (k-m+0.5)*exp(logPhiLower-logPhiDifference)*dlogPhiLower)/s/s2;
		//~ }

		//~ double logPhiMin, dlogPhiMin;
	    //~ logphi(-(PipelineMinVisits-m+0.5)/s,logPhiMin, dlogPhiMin);
	    //~ populationValues[i] = logPhiDifference - logPhiMin + log(populations[i].Fraction);
	    //~ dlpmf_dm -= dlogPhiMin/s;
	    //~ dlpmf_ds2 -= 0.5*(PipelineMinVisits-m+0.5)*dlogPhiMin/s/s2;

	
		double logPhi, dlogPhi;
	    logphi(-(PipelineMinVisits-m)/s,logPhi, dlogPhi); 
	    populationValues[i] = -0.5*log(2.0*M_PI*s2) - 0.5*(k-m)*(k-m)/s2 - logPhi + log(populations[i].Fraction);
	    double dlpmf_dm = (k-m)/s2 - dlogPhi/s;
	    double dlpmf_ds2 = 0.5*((k-m)*(k-m)/s2 - 1.0 - (PipelineMinVisits-m)*dlogPhi/s)/s2;
	   

	    for(int j = 0; j < probslen; ++j)
	    {
	        populationGradients[i][j] = dlpmf_dm + (1.0-2.0*probs[j] + mGradientFactor*(populations[i].LinearVariance + 2*m*populations[i].QuadraticVariance))*dlpmf_ds2;
	    }
    }
    
    double value = VerySmallLog;
    
    for (int j = 0; j < populations.size(); ++j)
    {
		
		value = log_add_exp(value, populationValues[j]);
    }

    for (int i = 0; i < probslen; ++i)
	{
		double temp = 0;
		for (int j = 0; j < populations.size(); ++j)
		{
			temp += exp(populationValues[j] - value) * populationGradients[j][i];
		}
		gradient[i] = temp;
	}
    return value;
}

double  sigmoid(double x)
{
   // return 0.5*(1.0+tanh(0.5*x));
    if (x < 0.0)
    {
        double a = exp(x);
        return a / (1.0 + a); 
    }
    else
    {
        return 1.0 / (1.0 + exp(-x));
    }
}

/*
double  elu(double x)
{
	return exp(-x);
}

double  elu_grad(double x, double elu_x)
{
	return -elu_x;
}
*/

double  elu(double x)
{
    if (x < density_cut)
    {
        return (1.0 + density_cut - x)*expm_density_cut; 
    }
    else
    {
        return exp(-x);
    }
}

double  elu_grad(double x, double elu_x)
{
    if (x < density_cut)
    {
        return -expm_density_cut; 
    }
    else
    {
        return -elu_x;
    }
}
