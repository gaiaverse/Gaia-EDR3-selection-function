#include "LogLikelihood.h"

LogLikelihood::LogLikelihood(const std::vector<Star> &data, int id): Data(data,id)
{

	Value = 0.0;
	Gradient = Eigen::VectorXd::Zero(totalTransformedParams);
}

void LogLikelihood::Calculate(Eigen::VectorXd& x)
{

	Reset();	

	for (int i = 0; i < Data.NStars; ++i)
	{
		PerStarContribution(i,x);
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


void LogLikelihood::PerStarContribution(int star, Eigen::VectorXd& x)
{

	Star candidate = Data.Stars[star];
	int k = candidate.nMeasure;
	int n = candidate.nVisit;
	

	//generate p vectors
	for (int i = 0; i < n; ++i)
	{
		int t= candidate.TimeSeries[i];

		double xt = x[Data.time_mapping[t]];
		int idx1 = Nt + Data.healpix_fov_1[t] * Nm + candidate.gBin;
		int idx2 = Nt + Data.healpix_fov_2[t] * Nm + candidate.gBin;
		double xml1 = x[idx1];
		double xml2 = x[idx2];
		
		Data.pt[i] = sigmoid(xt);
		Data.pml[i] = sigmoid(xml1 + xml2);
		Data.p[i] = Data.pt[i] *Data.pml[i];
	}
	
	bool poissonOverride = false;
	double zeroMeasureKiller = 0;
	double nMeasureKiller = 0;
	double likelihood = 0;
	double correction = 1.0;
	double log_likelihood = 0;
	double log_correction = 0.0;

	// probability black magic stuff
	poisson_binomial_pmf_forward(Data.p,n,Data.pmf_forward);
	for (int i = 0; i < PipelineMinVisits; ++i)
	{
		correction -= Data.pmf_forward[n-1][i];
	}
	
	if (correction < VerySmallNumber)
	{
		correction = 0;
		for (int i = PipelineMinVisits; i < n; ++i)
		{
			correction += Data.pmf_forward[n-1][i];
		}
	}

	bool longFlag = true;
	


	if ( (Data.pmf_forward[n-1][0] > VerySmallNumber) && (Data.pmf_forward[n-1][n] > VerySmallNumber))
	{
		longFlag = false;
		poisson_binomial_pmf_backward(Data.p,n,Data.pmf_backward);
		poisson_binomial_subpmf(PipelineMinVisits-1,n,Data.pmf_forward,Data.pmf_backward,Data.subpmf[0]);

		
		if (k > 0)
		{
			poisson_binomial_subpmf(k-1,n,Data.pmf_forward,Data.pmf_backward,Data.subpmf[1]);
			zeroMeasureKiller = 1;
		}
		
		if (k < n)
		{
			poisson_binomial_subpmf(k,n,Data.pmf_forward,Data.pmf_backward,Data.subpmf[2]);
			nMeasureKiller = 1;
		}
		
		likelihood = Data.pmf_forward[n-1][k];
		log_likelihood = log(likelihood);
		log_correction = log(correction);
		
		if (std::isnan(log_likelihood) || std::isnan(log_correction))
		{
			longFlag = true;
		}
		else
		{
			Value += log_likelihood - log_correction;
		}
	}

	
	if (longFlag)
	{
		poissonOverride = true;
		poisson_binomial_lpmf_forward(Data.p,n,Data.pmf_forward);
		poisson_binomial_lpmf_backward(Data.p,n,Data.pmf_backward);
		poisson_binomial_sublpmf(PipelineMinVisits-1,n,Data.pmf_forward,Data.pmf_backward,Data.subpmf[0]);

		if (k > 0)
		{
			poisson_binomial_sublpmf(k-1,n,Data.pmf_forward,Data.pmf_backward,Data.subpmf[1]);
			zeroMeasureKiller = 1;
		}
		
		if (k < n)
		{
			//~ std::cout << "trigger2" << std::endl;
			poisson_binomial_sublpmf(k,n,Data.pmf_forward,Data.pmf_backward,Data.subpmf[2]);
			nMeasureKiller = 1;
		}
		

		log_likelihood = Data.pmf_forward[n-1][k];

		for (int i = 0; i < PipelineMinVisits; ++i)
		{
			log_correction += log1p(-exp(Data.pmf_forward[n-1][i]-log_correction));
		}

		if (std::isnan(log_correction))
		{
			log_correction = -9999999;
			for (int i=PipelineMinVisits; i < n; ++i)
			{
				log_correction += log1p(exp(Data.pmf_forward[n-1][i] - log_correction));
			}
		}
		Value += log_likelihood - log_correction;
		likelihood = exp(log_likelihood);
		correction = exp(log_correction);
	}
	
			
	if (std::isnan(Value))
	{
		std::cout << "\n\n Error detected! NaN found in Value calculation on core " << Data.ID << " whilst calculating star " << star <<"\n";
		std::cout << "n = " << n << "k = " << k << std::endl;
		std::cout << "likelihood = " << likelihood << "  correction = " << correction << "\n";
		std::cout << "Longform calculation flag = " << longFlag << "\n";
		std::cout << "p = (";
		for (int i = 0; i < n; ++i)
		{
			std::cout << Data.p[i] << ", ";
		}
		std::cout << ")\n\npmf = (";
		for (int i = 0; i < n; ++i)
		{
			double v = Data.pmf_forward[n-1][i];
			if (longFlag)
			{
				v= exp(v);
			}
			std::cout << v << ",";
		}
		std::cout << ") \n\n\n";
		ERROR(100, "See above output");
	}
	
	
	for (int i = 0; i < n; ++i)
	{
		
		double dFdP_p;
		if (poissonOverride)
		{
			dFdP_p = Data.p[i] * ( exp(Data.subpmf[1][i]-log_likelihood)*zeroMeasureKiller - exp(Data.subpmf[2][i]-log_likelihood)*nMeasureKiller - exp(Data.subpmf[0][i] - log_correction));
		}
		else
		{
			dFdP_p = Data.p[i] * (   (Data.subpmf[1][i]*zeroMeasureKiller-Data.subpmf[2][i]*nMeasureKiller)/likelihood - Data.subpmf[0][i]/correction );
		}

		int t = candidate.TimeSeries[i];

		Gradient[Data.time_mapping[t]] += dFdP_p * (1.0 - Data.pt[i]);
		
		int offset = Nt + candidate.gBin;
		double mlGrad = dFdP_p * (1.0 - Data.pml[i]);
		
		int index1 = offset +  Data.healpix_fov_1[t] * Nm;
		int index2 = offset +  Data.healpix_fov_2[t] * Nm;
		
		Gradient[index1] += mlGrad;
		Gradient[index2] += mlGrad;
	}
	

}


          
void inline poisson_binomial_pmf_forward(std::vector<double> &  probs, int probslen, std::vector<std::vector<double>> & result)
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
