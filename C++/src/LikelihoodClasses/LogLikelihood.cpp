#include "LogLikelihood.h"

LogLikelihood::LogLikelihood(const std::vector<Star> &data, std::vector<int> & magBins, int dimension, int id): Data(data)
{
	//is magbins still needed?
	ID = id;
	Value = 0.0;
	Gradient = Eigen::VectorXd::Zero(dimension);

	
	//~ int suitablyLargeNumber = 1024; // a number at least as large as the maximum number of entries in a single star's time series
	//~ pmf = std::vector<double>(suitablyLargeNumber,0.0);
	//~ subpmf = std::vector<double>(suitablyLargeNumber,0.0);

	
	pmf_forward = std::vector<std::vector<double>>(suitablyLargeNumber,std::vector<double>(suitablyLargeNumber,0));
	pmf_backward =  std::vector<std::vector<double>>(suitablyLargeNumber,std::vector<double>(suitablyLargeNumber,0));
	subpmf =  std::vector<std::vector<double>>(3,std::vector<double>(suitablyLargeNumber,0));
	
    
    
    std::string healpix_fov_file = "../../ModelInputs/scanninglaw_to_healpix_"+std::to_string(healpix_order)+".csv";
    std::string needlet_file = "../../ModelInputs/needlets_"+std::to_string(healpix_order)+"_"+std::to_string(needlet_order)+".csv";
	healpix_fov_1 = std::vector<int>(TotalScanningTime,0);
	healpix_fov_2 = std::vector<int>(TotalScanningTime,0);


	time_mapping = std::vector<int>(TotalScanningTime,0);
	
    double time_ratio = 1;
    if (Nt < TotalScanningTime)
    {
		time_ratio = (double)Nt/TotalScanningTime;
	}

	for (int i = 0; i < TotalScanningTime; ++i)
	{
		time_mapping[i] = floor(time_ratio*i);
	}

    int i = 0;
    
    forLineVectorInFile(healpix_fov_file,',',
    
		if (i > 0)
		{
	        healpix_fov_1[i] = std::stoi(FILE_LINE_VECTOR[1]);
	        healpix_fov_2[i] = std::stoi(FILE_LINE_VECTOR[2]);  
	        
	    }
	     ++i;
    );
    
    i = 0;
    forLineVectorInFile(needlet_file,',',
 
		if (i > 0)
		{
	        needlet_u.push_back(std::stoi(FILE_LINE_VECTOR[0]));
	        needlet_v.push_back(std::stoi(FILE_LINE_VECTOR[1]));
	        needlet_w.push_back(std::stoi(FILE_LINE_VECTOR[2]));
		}
        ++i;
    );    
    needletN = needlet_u.size();
	pt = std::vector<double>(suitablyLargeNumber,0);
	pml = std::vector<double>(suitablyLargeNumber,0);
	p = std::vector<double>(suitablyLargeNumber,0);
}

void LogLikelihood::Calculate(Eigen::VectorXd& x)
{

	Reset();	

	for (int i = 0; i < Data.size(); ++i)
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

	Star candidate = Data[star];
	int k = candidate.nMeasure;
	int n = candidate.nVisit;
	

	//generate p vectors
	for (int i = 0; i < n; ++i)
	{
		int t= candidate.TimeSeries[i];

		double xt = x[time_mapping[t]];
		int idx1 = Nt + healpix_fov_1[t] * Nm + candidate.gBin;
		int idx2 = Nt + healpix_fov_2[t] * Nm + candidate.gBin;
		double xml1 = x[idx1];
		double xml2 = x[idx2];
		
		pt[i] = sigmoid(xt);
		pml[i] = sigmoid(xml1 + xml2);
		p[i] = pt[i] *pml[i];
	}
	
	bool poissonOverride = false;
	double zeroMeasureKiller = 0;
	double nMeasureKiller = 0;
	double likelihood = 0;
	double correction = 1.0;
	double log_likelihood = 0;
	double log_correction = 0.0;

	// probability black magic stuff
	poisson_binomial_pmf_forward(p,n,pmf_forward);
	for (int i = 0; i < PipelineMinVisits; ++i)
	{
		correction -= pmf_forward[n-1][i];
	}
	
	if (correction < verySmallNumber)
	{
		correction = 0;
		for (int i = PipelineMinVisits; i < n; ++i)
		{
			correction += pmf_forward[n-1][i];
		}
	}

	bool longFlag = true;
	


	if ( (pmf_forward[n-1][0] > verySmallNumber) && (pmf_forward[n-1][n] > verySmallNumber))
	{
		longFlag = false;
		poisson_binomial_pmf_backward(p,n,pmf_backward);
		poisson_binomial_subpmf(PipelineMinVisits-1,n,pmf_forward,pmf_backward,subpmf[0]);

		
		if (k > 0)
		{
			poisson_binomial_subpmf(k-1,n,pmf_forward,pmf_backward,subpmf[1]);
			zeroMeasureKiller = 1;
		}
		
		if (k < n)
		{
			poisson_binomial_subpmf(k,n,pmf_forward,pmf_backward,subpmf[2]);
			nMeasureKiller = 1;
		}
		
		likelihood = pmf_forward[n-1][k];
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
		poisson_binomial_lpmf_forward(p,n,pmf_forward);
		poisson_binomial_lpmf_backward(p,n,pmf_backward);
		poisson_binomial_sublpmf(PipelineMinVisits-1,n,pmf_forward,pmf_backward,subpmf[0]);

		if (k > 0)
		{
			poisson_binomial_sublpmf(k-1,n,pmf_forward,pmf_backward,subpmf[1]);
			zeroMeasureKiller = 1;
		}
		
		if (k < n)
		{
			//~ std::cout << "trigger2" << std::endl;
			poisson_binomial_sublpmf(k,n,pmf_forward,pmf_backward,subpmf[2]);
			nMeasureKiller = 1;
		}
		

		log_likelihood = pmf_forward[n-1][k];

		for (int i = 0; i < PipelineMinVisits; ++i)
		{
			log_correction += log1p(-exp(pmf_forward[n-1][i]-log_correction));
		}

		if (std::isnan(log_correction))
		{
			log_correction = -9999999;
			for (int i=PipelineMinVisits; i < n; ++i)
			{
				log_correction += log1p(exp(pmf_forward[n-1][i] - log_correction));
			}
		}
		Value += log_likelihood - log_correction;
		likelihood = exp(log_likelihood);
		correction = exp(log_correction);
	}
	
			
	if (std::isnan(Value))
	{
		std::cout << "\n\n Error detected! NaN found in Value calculation on core " << ID << " whilst calculating star " << star <<"\n";
		std::cout << "n = " << n << "k = " << k << std::endl;
		std::cout << "likelihood = " << likelihood << "  correction = " << correction << "\n";
		std::cout << "Longform calculation flag = " << longFlag << "\n";
		std::cout << "p = (";
		for (int i = 0; i < n; ++i)
		{
			std::cout << p[i] << ", ";
		}
		std::cout << ")\n\npmf = (";
		for (int i = 0; i < n; ++i)
		{
			double v = pmf_forward[n-1][i];
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
			dFdP_p = p[i] * ( exp(subpmf[1][i]-log_likelihood)*zeroMeasureKiller - exp(subpmf[2][i]-log_likelihood)*nMeasureKiller - exp(subpmf[0][i] - log_correction));
		}
		else
		{
			dFdP_p = p[i] * (   (subpmf[1][i]*zeroMeasureKiller-subpmf[2][i]*nMeasureKiller)/likelihood - subpmf[0][i]/correction );
		}
		GlobalDebug(1,
			std::cout << subpmf[1][i]-log_likelihood << "  " << subpmf[2][i] - log_likelihood << "   " << subpmf[0][i] - log_correction << "   " << dFdP_p << "\n";
		);
		int t = candidate.TimeSeries[i];

		Gradient[time_mapping[t]] += dFdP_p * (1.0 - pt[i]);
		
		int offset = Nt + candidate.gBin;
		double mlGrad = dFdP_p * (1.0 - pml[i]);
		
		int index1 = offset +  healpix_fov_1[t] * Nm;
		int index2 = offset +  healpix_fov_2[t] * Nm;
		
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
