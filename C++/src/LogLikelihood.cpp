#include "LogLikelihood.h"

LogLikelihood::LogLikelihood(const std::vector<Star> &data, std::vector<int> & magBins, int dimension, int id): Data(data)
{
	//is magbins still needed?
	ID = id;
	Value = 0.0;
	Gradient = Eigen::VectorXd::Zero(dimension);

	
	int suitablyLargeNumber = 1024; // a number at least as large as the maximum number of entries in a single star's time series
	pmf = std::vector<double>(suitablyLargeNumber,0.0);
	subpmf = std::vector<double>(suitablyLargeNumber,0.0);
	
	
    
    
    std::string healpix_fov_file = "../../../ModelInputs/scanninglaw_to_healpix_"+std::to_string(healpix_order)+".csv";
    std::string needlet_file = "../../../ModelInputs/needlets_"+std::to_string(healpix_order)+"_"+std::to_string(needlet_order)+".csv";
	std::vector<int> healpix_fov_1 = std::vector<int>(0,NT);
	std::vector<int> healpix_fov_2 = std::vector<int>(0,NT);

    int i = 0;
    forLineVectorInFile(healpix_fov_file,',',
    
        healpix_fov_1[i] = std::stoi(FILE_LINE_VECTOR[1]);
        healpix_fov_2[i] = std::stoi(FILE_LINE_VECTOR[2]);
        ++i;
    
    );
    
    forLineVectorInFile(needlet_file,',',
 
        needlet_u.push_back(std::stoi(FILE_LINE_VECTOR[0]));
        needlet_v.push_back(std::stoi(FILE_LINE_VECTOR[1]));
        needlet_w.push_back(std::stoi(FILE_LINE_VECTOR[2]));
    );    
    needletN = needlet_u.size();

    std::vector<int> time_mapping = std::vector<int>(0,NT);
    double time_ratio = (double)Nt / (double)NT;
    if (Nt == NT){
        for (int i = 0; i < NT; ++i)
        {
            time_mapping[i] = i;
        }
    }
    else{
        for (int i = 0; i < NT; ++i)
        {
            time_mapping[i] = time_ratio*i;
        }
    }
    
}

void LogLikelihood::Calculate(Eigen::VectorXd& x)
{
	std::cout << "\t\t\tThe log-likelihood function has been called" << std::endl;
	Reset();	
	

	for (int i = 0; i < Data.size(); ++i)
	{
		if (ID == 0)
		{
			std::cout << "\t\tCalculating contribution from star " << i << std::endl;
		}
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

double sigmoid(double x)
{
    return 0.5*(1.0+tanh(0.5*x));
}

void LogLikelihood::PerStarContribution(int star, Eigen::VectorXd& x)
{
	Star candidate = Data[star];

	int k = candidate.nMeasure;
	int n = candidate.nVisit;
	
	//copies in-place into pmf
	
	std::vector<double> pt = std::vector<double>(0,n); ///what is p?
	std::vector<double> pml = std::vector<double>(0,n);
	std::vector<double> p = std::vector<double>(0,n);
	
	for (int i = 0; i < n; ++i)
	{
		int t= candidate.TimeSeries[i];
		double xt = x[time_mapping[t]];
		double xml1 = x[Nt + healpix_fov_1[t] * Nm + candidate.gBin];
		double xml2 = x[Nt + healpix_fov_2[t] * Nm + candidate.gBin];
		
		pt[i] = sigmoid(xt);
		pml[i] = sigmoid(xml1 + xml2);
		
		p[i] = pt[i] *pml[i];
	}
	

	//modifies pmf and subpmf in place to set them to nice values
	direct_convolution_local(p,n,pmf);

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
	
	
	for (int i = 0; i < n; ++i)
	{
		
		double inv_p = 1.0/(1 - p[i]);
		
		subpmf[0] = pmf[0] * inv_p;
		for (int j = 1; j < n; ++j)
		{
			subpmf[j] = (pmf[j] - subpmf[j-1]*p[i])*inv_p;
		}
		subpmf[n-1] = pmf[n]/p[i];
		
		double dFdP = (gradient_first_term*subpmf[k-1]-gradient_second_term*subpmf[k])/likelihood - subpmf[PipelineMinVisits-1]/correction;
		
		
		
				
		
		int t= candidate.TimeSeries[i];

		Gradient[time_mapping[t]] += dFdP * p[i] * (1.0 - pt[i]);
		
		int offset = Nt + candidate.gBin;
		double mlGrad = dFdP * p[i] * (1.0 - pml[i]);
		Gradient[offset +  healpix_fov_1[t] * Nm] += mlGrad;
		Gradient[offset +  healpix_fov_2[t] * Nm] += mlGrad;
	}
}



            
            
void direct_convolution_local(std::vector<double> &  probs, int probslen, std::vector<double> & result)
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

