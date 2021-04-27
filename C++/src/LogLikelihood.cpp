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

double inline sigmoid(double x)
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

void LogLikelihood::PerStarContribution(int star, Eigen::VectorXd& x)
{

	//~ pmf_forward = std::vector<std::vector<double>>(suitablyLargeNumber,std::vector<double>(suitablyLargeNumber,999));
	//~ pmf_backward =  std::vector<std::vector<double>>(suitablyLargeNumber,std::vector<double>(suitablyLargeNumber,999));
	//~ subpmf =  std::vector<std::vector<double>>(3,std::vector<double>(suitablyLargeNumber,999));
	Star candidate = Data[star];

	int k = candidate.nMeasure;
	int n = candidate.nVisit;
	
	//copies in-place into pmf
	
	
	//~ std::cout << "\t\t\t p = (";
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

		//~ std::cout << p[i] << ", ";
	}
	//~ std::cout << ")"  << std::endl;
	
	// probability black magic stuff
	//direct_convolution_local(p,n,pmf);
	poisson_binomial_pmf_forward(p,n,pmf_forward);
	poisson_binomial_pmf_backward(p,n,pmf_backward);
	poisson_binomial_subpmf(PipelineMinVisits-1,n,pmf_forward,pmf_backward,subpmf[0]);

	
	double zeroMeasureKiller = 1;
	if (k > 0)
	{
		poisson_binomial_subpmf(k-1,n,pmf_forward,pmf_backward,subpmf[1]);
	}
	else
	{
		zeroMeasureKiller = 0;
	}

	double nMeasureKiller = 1;
	if (k < n)
	{
		//~ std::cout << "trigger2" << std::endl;
		poisson_binomial_subpmf(k,n,pmf_forward,pmf_backward,subpmf[2]);
	}
	else
	{
		nMeasureKiller = 0;
	}
	
	
	


	// Might as well use both
	double likelihood = pmf_forward[n-1][k];
	double correction = 1.0;
	
	for (int i = 0; i < PipelineMinVisits; ++i)
	{
		correction -= pmf_forward[n-1][i];
	}
	
	if (correction < 0)
	{
		std::cout << "ERROR! Correction went negative" << std::endl;
		exit(1);
	}
	
	
	
	Value += log(likelihood / correction);
	//~ int w = 20;
	//~ std::cout << std::setw(w) << star <<  std::setw(w) << Value << std::setw(w) << likelihood << std::setw(w)<< correction << std::endl;
	
	//~ std::cout << "\n i \t pmf	\t subpmf0 \t subpmf1 \t subpmf2 \n";
	//~ for (int i = 0; i < n; ++i)
	//~ {
		//~ std::cout << std::setw(w) << i << std::setw(w) << pmf_forward[n-1][i] << std::setw(w)<< subpmf[0][i] << std::setw(w)<< subpmf[1][i]<< std::setw(w) << subpmf[2][i] << "\n"; 
	//~ }
	//~ std::cout << "table ends\n\n";
	for (int i = 0; i < n; ++i)
	{
		double dFdP_p = p[i] * (   (subpmf[1][i]*zeroMeasureKiller-subpmf[2][i]*nMeasureKiller)/likelihood - subpmf[0][i]/correction );
	
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

void inline poisson_binomial_pmf_backward(std::vector<double> &  probs, int probslen, std::vector<std::vector<double>> & result)
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

void inline poisson_binomial_subpmf(int m, int probslen, std::vector<std::vector<double>> & pmf_forward, std::vector<std::vector<double>> & pmf_backward, std::vector<double> & result)
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
