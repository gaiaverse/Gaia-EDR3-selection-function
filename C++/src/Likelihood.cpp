#include "Likelihood.h"

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




Likelihood::Likelihood(const std::vector<Star> &data, std::vector<int> & magBins, int dimension, int id): Data(data)
{
	ID = id;
	Value = 0.0;
	Gradient = Eigen::VectorXd::Zero(dimension);
	MagBins = magBins;
	
	MinVisits = 5; //hard-coded parameter, the number of times a star has to be detected for it to enter gaia pipeline
    
    mu_mean = -3.0;
    mu_variance = 1.0;
    lg = 0.1;

    Kg_decomposed = false;
	
	int suitablyLargeNumber = 1024; // a number at least as large as the maximum number of entries in a single star's time series
	pmf = std::vector<double>(suitablyLargeNumber,0.0);
	subpmf = std::vector<double>(suitablyLargeNumber,0.0);
	
	//initialise an array of nbins per nt to store the modified values of the position vector each time step
	
	perBinP = std::vector<std::vector<double>>(magBins.size(),std::vector<double>(Nt,0.0));
}

void Likelihood::Calculate(Eigen::VectorXd& x)
{

	Reset();	
	

	
	GeneratePs(x);
	
	
	//~ std::vector<double> probs;
	if (ID == 0)
	{
		Prior(x);

	}

	for (int i = 0; i < Data.size(); ++i)
	{
		if (ID == 0)
		{
			std::cout << "\t\tCalculating contribution from star " << i << std::endl;
		}
		PerStarContribution(i);
	}

	

}

void Likelihood::GeneratePs(Eigen::VectorXd&x)
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
	std::cout << "\tProbs generated" << std::endl;
}
void Likelihood::PerStarContribution(int star)
{
	Star candidate = Data[star];

	int k = candidate.nMeasure;
	int n = candidate.nVisit;
	
	//copies in-place into pmf
	direct_convolution_local(perBinP[candidate.gBin],candidate.TimeSeries,n,pmf);

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
		
		double dFdP = (gradient_first_term*subpmf[k-1]-gradient_second_term*subpmf[k])/likelihood - subpmf[MinVisits-1]/correction;
		double dPdX = p * (1 - p);
		int t= candidate.TimeSeries[i];

		Gradient[offset + t] = dFdP * dPdX;
	}
}



void Likelihood::Reset()
{
	Value = 0;
	for (int i = 0; i < Gradient.size(); ++i)
	{
		Gradient[i] = 0;
	}
	std::cout << "\t Rest values + ready to go!" << std::endl;
}

void Likelihood::Prior(Eigen::VectorXd& params)
{

    // Unpack parameters
    //double lt = exp(params(0));
    //double lg = exp(params(1));
    //double sigma2 = exp(params(2));
    //double m = exp(params(3));
    //double tau2 = exp(params(4));
    double lt = exp(params(0));
    double sigma2 = exp(params(1));

    
    VectorXd mu = params.segment(Nh, Ng);
    VectorXd x = params.segment(Nh+Ng, Ng*Nt);
    
    // Apply the priors on the hyper-hyper-parameters
    
    //old priors
    //~ PriorLengthscale(lt,  0);
    //~ PriorLengthscale(lg,  1);
    //~ PriorVariance(sigma2, 2);
    //~ PriorLengthscale(m,   3);
    //~ PriorVariance(tau2,   4);
    
    //new priors
    PriorLengthscale(lt,  0);
    PriorVariance(sigma2, 1);
    
    std::cout << "\t Hyper-Priors Calculated" << std::endl;
    // Apply the prior on the hyper-parameters
    PriorMu(mu);
	std::cout << "\t Prior Mu Calculated" << std::endl;
	
		std::cout << "\t Beginning PriorX" << std::endl;
    // Apply the prior on the parameters
	PriorX(x, mu, lt, sigma2);
	std::cout << "\t Prior X Calculated" << std::endl;

}

void Likelihood::PriorLengthscale(double lengthscale, int param_index)
{
    // Implements an InverseGamma(1,2) prior for the lengthscales
    
    // We use this three times - division is evil
    double two_over_lengthscale = 2.0/lengthscale;
    
    // lnA = np.log(2.0)-2.0*np.log(l)-2.0/l
    Value += log(two_over_lengthscale/lengthscale) - two_over_lengthscale;
    
    // dlnAdl = 2.0*(1.0-l)/l/l, return dlnAdlnl = l*dlnAdl
    Gradient[param_index] += two_over_lengthscale*(1.0-lengthscale);
}

void Likelihood::PriorVariance(double variance, int param_index)
{
    // Implements a Gamma(1,1) prior for the variances
    
    // lnS = - sigma2
    Value -= variance;
    
    // dlnSdsigma2 = - 1.0, return dlnSdlnsigma2 = sigma2*dlnSdsigma2
    Gradient[param_index] -= variance;
}

void Likelihood::PriorMu(Eigen::VectorXd& mu)
{
    // Implements the Gaussian Process prior on mu ~ N(mu_mean,mu_variance)
    
    VectorXd diff = mu.array() - mu_mean;
    Value += - 0.5 * Ng * log( 2.0 * M_PI * mu_variance ) - 0.5 * diff.squaredNorm() / mu_variance;

    for (int ig = 0; ig < Ng; ig++) 
    {
		
        Gradient[Nh+ig] -= diff[ig] / mu_variance;
    }
 
    
}

void Likelihood::PriorX(Eigen::VectorXd& x, Eigen::VectorXd& mu, double lt, double sigma2)
{
    // Implements the Gauss-Markov prior on x
    std::cout << "\t Beginning priorX" << std::endl;
    // Useful shorthands
    double u = exp(-1.0/lt);
    double u2 = u*u;
    double oneplusu2 = 1.0 + u2;
    double oneoveroneminusu2 = 1.0 / ( 1.0 - u2 );
    std::cout << u2 << oneplusu2 << std::endl;
    // Reshape to form Y
    // Y = x.reshape((Ng,Nt)) - mu.reshape((Ng,1))
    
    /*
    std::cout << "\t Attempting to initialise a stupendously big matrix" << std::endl;
    
    Eigen::Matrix<double, Ng, Nt> Y;
    int signpost = 0;
     std::cout << "\t Did it!" << std::endl;
    for (int i = 0; i < Ng; i++) 
    {
        Y.row(i) = x.segment(i*Nt,Nt).array() - mu[i];
    }
    std::cout << "Signpost " << signpost << std::endl; ++signpost;
    if (Kg_decomposed == false){
        
        std::cout << "Building Kg" << std::endl;
        // Build covariance matrix
        for (int i = 0; i < Ng; i++) 
        {
            for (int j = 0; j < i; j++) 
            {
                Kg(i,j) = Kg(j,i) = exp(-pow(i - j,2)/(2.0*lg*lg));
            }
            Kg(i,i) = 1.0 + SingularityPreventer;
        }
         std::cout << "Signpost " << signpost << std::endl; ++signpost;
        // Householder decomposition (with pivoting!)
        Eigen::FullPivHouseholderQR<Matrix<double, Ng,Ng>> decomp  = Kg.fullPivHouseholderQr();
        
        // Compute quantities we will need later
        invKg = decomp.inverse();
        logdetKg = decomp.logAbsDeterminant();
		std::cout << "Signpost " << signpost << std::endl; ++signpost;
        // Set flag so we don't do this again
        Kg_decomposed = true;
    }

	 std::cout << "Signpost " << signpost << std::endl; ++signpost;
    MatrixXd invKgY = invKg*Y;
    double logdetinvKt = (Nt-1.0)*log( oneoveroneminusu2 );
    double TrJt = -2.0*(Nt-1.0)*u2*oneoveroneminusu2/(lt*lt);
    

    // Compute invKgYinvKt
    Matrix<double, Ng, Nt> invKgYinvKt;
    for (int ig = 0; ig < Ng; ig++) 
    {
        invKgYinvKt(ig,0) = ( invKgY(ig,0) - u * invKgY(ig,1) )*oneoveroneminusu2;
        invKgYinvKt(ig,Nt-1) = ( invKgY(ig,Nt-1) - u * invKgY(ig,Nt-2) )*oneoveroneminusu2;
        for (int it = 1; it < Nt-1; it++) 
        {
            invKgYinvKt(ig,it) = ( oneplusu2 * invKgY(ig,it) - u * ( invKgY(ig,it-1) + invKgY(ig,it+1) ) )*oneoveroneminusu2;
        }
    }
     std::cout << "Signpost " << signpost << std::endl; ++signpost;
    double Y_invKgYinvKt = (Y.array()*invKgYinvKt.array()).sum();
    
    
    // Compute YJt
    Matrix<double, Ng, Nt> YJt;
    int M = std::min(10+ceil(-lt*log(1e-16)),(double)Nt);
	 std::cout << "Signpost " << signpost << std::endl; ++signpost;
    std::vector<double> power_u = std::vector<double>(M);
    power_u[0] = u;
    for (int i = 1; i < M; i++) 
    {
        power_u[i] = u*power_u[i-1];
    }
    
	 std::cout << "Signpost " << signpost << std::endl; ++signpost;
    double res;
    for (int i = 0; i < Ng; i++) 
    {
        for (int l = 0; l < Nt; l++) 
        {
            
            res =  - oneplusu2*Y(i,l)*oneoveroneminusu2;
            
            if (l < M)
            {
                res += u2*Y(i,0)*power_u[l]*oneoveroneminusu2;
            }   
            if (Nt-l-1 < M)
            {
                res += u2*Y(i,Nt-1)*power_u[Nt-l-1]*oneoveroneminusu2;
            }
            for (int j = std::max(0,l-M); j < std::min(Nt,l+M); j++) 
            {
                res += Y(i,j)*power_u[abs(l-j)];
            }
            YJt(i,l) = res/(lt*lt);
        }
    }
                
	 std::cout << "Signpost " << signpost << std::endl; ++signpost;
    double YJt_invKgYinvKt = (YJt.array()*invKgYinvKt.array()).sum();
                
    // lnP = -Ng*Nt*np.log(2.0*np.pi*sigma2)/2.0 + Ng*logdetinvKt/2.0 - Nt*logdetKg/2.0 - Y_invKgYinvKt/2.0/sigma2
    Value += -Ng*Nt*log(2.0*M_PI*sigma2)/2.0 + Ng*logdetinvKt/2.0 - Nt*logdetKg/2.0 - Y_invKgYinvKt/(2.0*sigma2);
    
    //Matrix<double, Ng, Nt> dlnP_dX = -invKgYinvKt/sigma2;
    //dlnP_dx = dlnP_dX.ravel()
    // dlnP_dmu = -dlnP_dX.sum(axis=1)
    for (int ig = 0; ig < Ng; ig++) 
    {
		for (int t =0; t < Nt; ++t)
		{
			//Gradient.segment(Nh+Ng+ig*Nt,Nt).array() += -invKgYinvKt.row(ig)/sigma2;
			Gradient[Nh + Ng + ig*Nt + t] -= invKgYinvKt(ig,t)/sigma2;
        }
        Gradient[Nh+ig] += -Gradient.segment(Nh+Ng+ig*Nt,Nt).sum();
    }

    //dlnP_dlt = -Ng*TrJt/2.0 + YJt_invKgYinvKt/2.0/sigma2
    Gradient[0] += lt*(-Ng*TrJt/2.0 + YJt_invKgYinvKt/(2.0*sigma2));
    
    //dlnP_dlg = -Nt*TrJg/2.0 + JgTY_invKgYinvKt/2.0/sigma2
    //Gradient[1] += lg*(-Nt*TrJg/2.0 + JgTY_invKgYinvKt/(2.0*sigma2));
            
    //dlnP_dsigma2 = -Ng*Nt/2.0/sigma2 + Y_invKgYinvKt/2.0/sigma2/sigma2
    Gradient[1] += -Ng*Nt/2.0 + Y_invKgYinvKt/2.0/sigma2;
	 std::cout << "Signpost " << signpost << std::endl; ++signpost;
	 
	 */
}


            
            
