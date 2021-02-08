#include "LikelihoodPrior.h"

void LikelihoodPrior::Calculate(Eigen::VectorXd& x)
{
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


void LikelihoodPrior::Prior(Eigen::VectorXd& params)
{
	std::cout << "Prior Called" << std::endl;
    // Unpack parameters
    //double lt = exp(params(0));
    //double lg = exp(params(1));
    //double sigma2 = exp(params(2));
    //double m = exp(params(3));
    //double tau2 = exp(params(4));
    double lt = exp(params(0));
    double sigma2 = exp(params(1));

    
    //VectorXd mu = params.segment(Nh, Ng);
    Map<VectorXd> mu(params.segment(Nh, Ng).data(), Ng);
    //VectorXd x = params.segment(Nh+Ng, Ng*Nt);
    Map<VectorXd> x(params.segment(Nh + Ng, Ng * Nt).data(), Ng * Nt);
    
    // Apply the priors on the hyper-hyper-parameters
    
    //old priors
    //~ PriorLengthscale(lt,  0);
    //~ PriorLengthscale(lg,  1);
    //~ PriorVariance(sigma2, 2);
    //~ PriorLengthscale(m,   3);
    //~ PriorVariance(tau2,   4);
    
    //new priors
    //PriorLengthscale(lt,  0);
    //PriorVariance(sigma2, 1);
    
    //~ std::cout << "\t Hyper-Priors Calculated" << std::endl;
    // Apply the prior on the hyper-parameters
	PriorMu(mu);
	//~ std::cout << "\t Prior Mu Calculated" << std::endl;
	
		//~ std::cout << "\t Beginning PriorX" << std::endl;
    // Apply the prior on the parameters
	PriorX(x, mu, lt, sigma2);
	//~ std::cout << "\t Prior X Calculated" << std::endl;
	
	//~ std::cout << Value << std::endl;
	std::cout << "Prior Finished" << std::endl;
}

void LikelihoodPrior::PriorLengthscale(double lengthscale, int param_index)
{
    // Implements an InverseGamma(1,2) prior for the lengthscales
    
    // We use this three times - division is evil
    double two_over_lengthscale = 2.0/lengthscale;
    
    // lnA = np.log(2.0)-2.0*np.log(l)-2.0/l
    Value += log(two_over_lengthscale/lengthscale) - two_over_lengthscale;
    
    // dlnAdl = 2.0*(1.0-l)/l/l, return dlnAdlnl = l*dlnAdl
    Gradient[param_index] += two_over_lengthscale*(1.0-lengthscale);
}

void LikelihoodPrior::PriorVariance(double variance, int param_index)
{
    // Implements a Gamma(1,1) prior for the variances
    
    // lnS = - sigma2
    Value -= variance;
    
    // dlnSdsigma2 = - 1.0, return dlnSdlnsigma2 = sigma2*dlnSdsigma2
    Gradient[param_index] -= variance;
}

void LikelihoodPrior::PriorMu( Map<VectorXd> & mu)
{
    // Implements the Gaussian Process prior on mu ~ N(mu_mean,mu_variance)
    std::cout << "\tPrior Mu Called" << std::endl;
    VectorXd diff = mu.array() - mu_mean;
  
    Value += - 0.5 * Ng * log( 2.0 * M_PI * mu_variance ) - 0.5 * diff.squaredNorm() / mu_variance;

    for (int ig = 0; ig < Ng; ig++) 
    {
		
        Gradient[Nh+ig] -= diff[ig] / mu_variance;
    }
 
    
}

void LikelihoodPrior::PriorX( Map<VectorXd> & x,  Map<VectorXd> & mu, double lt, double sigma2)
{
	std::cout << "Prior X Called" << std::endl;
    // Implements the Gauss-Markov prior on x
   
    // Useful shorthands
    double u = exp(-1.0/lt);
    double u2 = u*u;
    double oneplusu2 = 1.0 + u2;
    double oneoveroneminusu2 = 1.0 / ( 1.0 - u2 );


    // Reshape to form Y
    // Y = x.reshape((Ng,Nt)) - mu.reshape((Ng,1))

    // Required for YJt
    int M = std::min(10+ceil(-lt*log(1e-16)),(double)Nt);
   
    std::vector<double> power_u = std::vector<double>(M);
    power_u[0] = u;
    for (int i = 1; i < M; i++) 
    {
        power_u[i] = u*power_u[i-1];
    }

    if (Kg_decomposed == false){
        

        // Build covariance matrix
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
	
        // Set flag so we don't do this again
        Kg_decomposed = true;
    }
	

    double logdetinvKt = (Nt-1.0)*log( oneoveroneminusu2 );
    double TrJt = -2.0*(Nt-1.0)*u2*oneoveroneminusu2/(lt*lt);
	
    // Compute invKgYinvKt
    
    VectorXd vec_Y(Nt);
    
    VectorXd vec_invKgY(Nt);
    VectorXd vec_invKgmu = invKg*mu;
    double sca_invKgYinvKt;
    double sca_YJt;
    double Y_invKgYinvKt = 0.0;
    double YJt_invKgYinvKt = 0.0;

	std::cout << "\t\tInitialisig stonking big vector" << std::endl;
	Eigen::Map<Eigen::Matrix<double, Ng,Nt, Eigen::RowMajor>> X(x.data());
    
    std::cout << "\t\tBeginning big loop" << std::endl;
    for (int ig = 0; ig < Ng; ig++) 
    {
		std::cout << "\t\t\t" << ig << std::endl;
        vec_Y = X.row(ig).array() - vec_invKgmu[ig];
        vec_invKgY = (invKg.row(ig)*X).array() - vec_invKgmu[ig];

        for (int it = 0; it < Nt; it++)
        {
			std::cout << "\t\t\t\t" << it << std::endl;
            // Compute invKgYinvKt(ig,it)
            switch (it)
            {
                case 0: 

					sca_invKgYinvKt = ( vec_invKgY[0] - u * vec_invKgY[1] )*oneoveroneminusu2;
					break;
                case Nt-1: 
					sca_invKgYinvKt = ( vec_invKgY[Nt-1] - u * vec_invKgY[Nt-2] )*oneoveroneminusu2;
					break;
                default: 
					sca_invKgYinvKt = ( oneplusu2 * vec_invKgY[it] - u * ( vec_invKgY[it-1] + vec_invKgY[it+1] ) )*oneoveroneminusu2;
					break;
            }

            // Increment Y_invKgYinvKt
            Y_invKgYinvKt += vec_Y[it] * sca_invKgYinvKt;

            // Increment gradient on x and mu
            Gradient[Nh + Ng + ig*Nt + it] -= sca_invKgYinvKt/sigma2;
            Gradient[Nh + ig] += sca_invKgYinvKt/sigma2;

            // Compute Yjt
            sca_YJt =  - oneplusu2*vec_Y[it]*oneoveroneminusu2;
            
            if (it < M)
            {
                sca_YJt += u2*vec_Y[0]*power_u[it]*oneoveroneminusu2;
            }   
            if (Nt-it-1 < M)
            {
                sca_YJt += u2*vec_Y[Nt-1]*power_u[Nt-it-1]*oneoveroneminusu2;
            }
            for (int jt = std::max(0,it-M); jt < std::min(Nt,it+M); jt++) 
            { 
                sca_YJt += vec_Y[jt]*power_u[abs(it-jt)];
            }

            // Increment YJt_invKgYinvKt, the lt*lt is because YJt is off by that factor
            YJt_invKgYinvKt += sca_YJt*sca_invKgYinvKt/(lt*lt);

        }
    }

    // Compute YJt
    
    // lnP = -Ng*Nt*np.log(2.0*np.pi*sigma2)/2.0 + Ng*logdetinvKt/2.0 - Nt*logdetKg/2.0 - Y_invKgYinvKt/2.0/sigma2
    Value += -Ng*Nt*log(2.0*M_PI*sigma2)/2.0 + Ng*logdetinvKt/2.0 - Nt*logdetKg/2.0 - Y_invKgYinvKt/(2.0*sigma2);

    //dlnP_dlt = -Ng*TrJt/2.0 + YJt_invKgYinvKt/2.0/sigma2
    Gradient[0] += lt*(-Ng*TrJt/2.0 + YJt_invKgYinvKt/(2.0*sigma2));
    
    //dlnP_dlg = -Nt*TrJg/2.0 + JgTY_invKgYinvKt/2.0/sigma2
    //Gradient[1] += lg*(-Nt*TrJg/2.0 + JgTY_invKgYinvKt/(2.0*sigma2));
            
    //dlnP_dsigma2 = -Ng*Nt/2.0/sigma2 + Y_invKgYinvKt/2.0/sigma2/sigma2
    Gradient[1] += -Ng*Nt/2.0 + Y_invKgYinvKt/2.0/sigma2;
	std::cout << "Prior X Completed" << std::endl;
}
