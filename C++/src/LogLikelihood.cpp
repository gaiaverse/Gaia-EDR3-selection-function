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

}

void LogLikelihood::Calculate(Eigen::VectorXd& x)
{
	//std::cout << "\t\t\tThe log-likelihood function has been called" << std::endl;
	Reset();	
	
	
	for (int i = 0; i < Data.size(); ++i)
	{
		PerStarContribution(i,x);
		
		
		//~ for (int j = 0; j <x.size(); ++j)
		//~ {
			//~ double d = x[j] - 5;
			//~ Value += -0.5*d*d;
			//~ Gradient[j] -= d;
		//~ }
	}
	
	//~ VectorXd xNudge;
	//~ double OldVal = Value;
	//~ VectorXd OldGrad = Gradient;
	//~ double delta = 1e-8;
	//~ for (int j = 0; j < x.size(); ++j)
	//~ {
		//~ Reset();
		
		//~ xNudge = x;
		//~ xNudge[j] += delta;
		
		//~ for (int i = 0; i < Data.size(); ++i)
		//~ {
			//~ PerStarContribution(i,xNudge);
		//~ }
		
		//~ double testGrad = (Value - OldVal)/delta;
		
		//~ std::cout << "\t\t\t\tInternal test: " << testGrad << " (direct calc = " << OldGrad[j] << ") \n"; 
	//~ }
	//~ Value = OldVal;
	//~ Gradient = OldGrad;
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
    //return 0.5*(1.0+tanh(0.5*x));
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
	Star candidate = Data[star];

	int k = candidate.nMeasure;
	int n = candidate.nVisit;
	
	//copies in-place into pmf
	
	std::vector<double> pt = std::vector<double>(n,0); ///what is p?
	std::vector<double> pml = std::vector<double>(n,0);
	std::vector<double> p = std::vector<double>(n,0);
	double p_sum = 0;

	
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
		p_sum += p[i];
		//~ std::cout << p[i] << ", ";
	}
	//~ std::cout << "), sum = " << p_sum << std::endl;
	// probability black magic stuff
	direct_convolution_local(p,n,pmf);

	double likelihood = pmf[k];
	
	double correction = 1.0;
	for (int i = 0; i < PipelineMinVisits; ++i)
	{
		correction -= pmf[i];
	}
	
	//~ std::cout << "pmf = (";
	//~ for (int i = 0; i < n; ++i)
	//~ {
		//~ std::cout << pmf[i] << ", ";
	//~ }
	//~ std::cout << ")\n\n";

	if (correction < 0)
	{
		std::cout << "CORRECTION ERROR in STAR " << star << " ( " << correction << ") " << std::endl;
		
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
	
	VectorXd perStarGrad = VectorXd::Zero(Gradient.size());
	int peak = n/2;
	
	double prevP = -99999999999;
	double prevdFdP = 0;
	for (int i = 0; i < n; ++i)
	{
		double dFdP;
		

		if ( abs(p[i] - prevP) < 1e-10)
		{
			dFdP = prevdFdP;
		}
		else
		{
			CalculatePMF(i,n,k,p);
			dFdP = (gradient_first_term*subpmf[k-1]-gradient_second_term*subpmf[k])/likelihood - subpmf[PipelineMinVisits-1]/correction;
			prevP = p[i];
			prevdFdP = dFdP;
		}
	
		int t= candidate.TimeSeries[i];

		Gradient[time_mapping[t]] += dFdP * p[i] * (1.0 - pt[i]);
		
		int offset = Nt + candidate.gBin;
		double mlGrad = dFdP * p[i] * (1.0 - pml[i]);
		
		int index1 = offset +  healpix_fov_1[t] * Nm;
		int index2 = offset +  healpix_fov_2[t] * Nm;
		
		Gradient[index1] += mlGrad;
		Gradient[index2] += mlGrad;
		
		perStarGrad[time_mapping[t]] = dFdP * p[i] * (1.0 - pt[i]);
		perStarGrad[index1] = mlGrad;
		perStarGrad[index2] = mlGrad;
	}
	
	if (perStarGrad.norm() > 1e10)
	{
		std::cout << "\n\nERROR: Star " << star << " has a gradient: " << perStarGrad.transpose() << std::endl;
		std::cout << "This is usually indicative of an error. ";
		
		if (QuitOnLargeGrad)
		{
			std::cout << "User options request termination on this condition.\n";
			exit(1);
		}
		std::cout << "Please proceed with caution " << std::endl;
	}

}



void inline LogLikelihood::CalculatePMF(int i,int n, int k,std::vector<double> & ps)
{
	double p = ps[i];
	
	
	double inv_p = 1.0/p;
	double inv_1mp = 1.0/(1.0 - p);
	
	
	bool needsExplicitCalculation = false;
	
	if (p*inv_1mp < 1)
	{
		
		subpmf[0] = pmf[0] * inv_1mp;
		for (int j = 1; j < n; ++j)
		{
			subpmf[j] = (pmf[j] - subpmf[j-1]*p)*inv_1mp;
		}
	}
	else 
	{
		subpmf[n-1] = pmf[n]*inv_p;
		for (int j = n-1; j > 0; --j)
		{
			subpmf[j-1] = (pmf[j] - subpmf[j]*(1.0-p))*inv_p;
		}
	}
	
}

std::vector<double> LogLikelihood::LikelihoodGivenP(std::vector<double> p, int n, int k)
{
	
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

