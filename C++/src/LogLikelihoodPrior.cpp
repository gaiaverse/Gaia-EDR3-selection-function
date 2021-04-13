#include "LogLikelihoodPrior.h"

void LogLikelihoodPrior::Calculate(Eigen::VectorXd& x)
{
	std::cout << "\t\t\tThe prior-log-likelihood function hass been called" << std::endl;
	Reset();		
	GeneratePs(x);

	//special bit for child class - call the Prior before the data likelihood. 
	//calls the prior first (that way if it crashes it happens quickly!)
	Prior(x);

	for (int i = 0; i < Data.size(); ++i)
	{
		PerStarContribution(i);
	}
}


void LogLikelihoodPrior::Prior(Eigen::VectorXd& params)
{
    // Unpack parameters

    Map<VectorXd> mu(params.segment(0, Nm).data(), Nm);

    Map<VectorXd> zt(params.segment(Nm, Nt).data(), Nt);

    Map<VectorXd> zml(params.segment(Nm + Nt, Nm*Nl).data(), Nm*Nl);

    
    

    PriorMu(mu);
	Prior(x, mu, lt, sigma2);

}

void LogLikelihoodPrior::PriorMu( Map<VectorXd> & mu)
{
    // Implements the Gaussian Process prior on mu ~ N(mu_mean,mu_variance)
    //mu_mean and mu_variance are currently global parameters
    
    VectorXd diff = mu.array() - mu_mean;
  
    Value += - 0.5 * Ng * log( 2.0 * M_PI * mu_variance ) - 0.5 * diff.squaredNorm() / mu_variance;

    for (int ig = 0; ig < Ng; ig++) 
    {
		
        Gradient[Nh+ig] -= diff[ig] / mu_variance;
    }
 
    
}

void LogLikelihoodPrior::PriorX( Map<VectorXd> & x,  Map<VectorXd> & mu, double lt, double sigma2)
{
    // Implements the Gauss-Markov prior on x
   
    // Useful shorthands
    double u = exp(-1.0/lt);
    double u2 = u*u;
    double onePlusu2 = 1.0 + u2;
    double oneOverOneMinusu2 = 1.0 / ( 1.0 - u2 );

    // M is a truncation length which reduces the order of the calculation to manageable sizes
    int M = std::min(10+ceil(-lt*log(1e-16)),(double)Nt);
   
    std::vector<double> power_u = std::vector<double>(M);
    power_u[0] = u;
    for (int i = 1; i < M; i++) 
    {
        power_u[i] = u*power_u[i-1];
    }

	 // Inverse Covariance matrix is constant + is small, so calculate once and store in memory
    if (Kg_decomposed == false)
    {
		MakeCovarianceMatrix();
        Kg_decomposed = true;
    }
	
    double logdetinvKt = (Nt-1.0)*log( oneOverOneMinusu2 );
    double TrJt = -2.0*(Nt-1.0)*u2*oneOverOneMinusu2/(lt*lt);
	
    // Compute invKgYinvKt
    //This is a slightly hacky computation, written to minimise memory usage + therefore no matrix library operations
    VectorXd Y(Nt);
    VectorXd invKgY(Nt);
    VectorXd invKgmu = invKg*mu;
    double invKgYinvKt;
    double YJt;
    double Y_invKgYinvKt = 0.0;
    double YJt_invKgYinvKt = 0.0;

    for (int ig = 0; ig < Ng; ig++) 
    {
		int offset = ig*Nt;
		for (int id = 0; id < Nt; ++id)
		{
			Y[id] = x[offset + id] - invKgmu[ig];
			
			invKgY[id] = -invKgmu[ig];
			for (int k = 0; k < Ng; ++k)
			{
				int xId = k*(Nt) + id;
				invKgY[id] += invKg(ig,k) * x[xId];
			}
		}

        for (int it = 0; it < Nt; it++)
        {

            switch (it)
            {
                case 0: 

					invKgYinvKt = ( invKgY[0] - u * invKgY[1] )*oneOverOneMinusu2;
					break;
                case Nt-1: 
					invKgYinvKt = ( invKgY[Nt-1] - u * invKgY[Nt-2] )*oneOverOneMinusu2;
					break;
                default: 
					invKgYinvKt = ( onePlusu2 * invKgY[it] - u * ( invKgY[it-1] + invKgY[it+1] ) )*oneOverOneMinusu2;
					break;
            }

            // Increment Y_invKgYinvKt
            Y_invKgYinvKt += Y[it] * invKgYinvKt;

            // Increment gradient on x and mu
            Gradient[Nh + Ng + ig*Nt + it] -= invKgYinvKt/sigma2;
            Gradient[Nh + ig] += invKgYinvKt/sigma2;

            // Compute Yjt
            YJt =  - onePlusu2*Y[it]*oneOverOneMinusu2;
            
            if (it < M)
            {
                YJt += u2*Y[0]*power_u[it]*oneOverOneMinusu2;
            }   
            if (Nt-it-1 < M)
            {
                YJt += u2*Y[Nt-1]*power_u[Nt-it-1]*oneOverOneMinusu2;
            }
            for (int jt = std::max(0,it-M); jt < std::min(Nt,it+M); jt++) 
            { 
                YJt += Y[jt]*power_u[abs(it-jt)];
            }

            // Increment YJt_invKgYinvKt, the lt*lt is because YJt is off by that factor
            YJt_invKgYinvKt += YJt*invKgYinvKt/(lt*lt);

        }
    }

    Value += -Ng*Nt*log(2.0*M_PI*sigma2)/2.0 + Ng*logdetinvKt/2.0 - Nt*logdetKg/2.0 - Y_invKgYinvKt/(2.0*sigma2);

    Gradient[0] += lt*(-Ng*TrJt/2.0 + YJt_invKgYinvKt/(2.0*sigma2));
    Gradient[1] += -Ng*Nt/2.0 + Y_invKgYinvKt/2.0/sigma2;

}

void LogLikelihoodPrior::MakeCovarianceMatrix()
{
	Eigen::Matrix<double, Ng, Ng> Kg;
	for (int i = 0; i < Ng; i++) 
	{
		for (int j = 0; j < i; j++) 
		{
			Kg(i,j) = Kg(j,i) = exp(-pow(i - j,2)/(2.0*lg*lg));
		}
		Kg(i,i) = 1.0 + SingularityPreventer;
	}

	// Householder decomposition (with pivoting!)
	Eigen::FullPivHouseholderQR<Matrix<double, Ng,Ng>> decomp  = Kg.fullPivHouseholderQr();
	
	// Compute quantities we will need later
	invKg = decomp.inverse();
	logdetKg = decomp.logAbsDeterminant();
}
