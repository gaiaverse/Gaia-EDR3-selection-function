#include "Liklihood.h"

Liklihood::Liklihood(const std::vector<Star> &data, int nPoints, int id): Data(data)
{
	ID = id;
	Value = 0.0;
	Gradient = std::vector<double>(nPoints,0.0);
	//NVists = 5
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

	}

	
}

void Liklihood::Prior(Eigen::VectorXd& params)
{
    
    // Unpack parameters
    double lt = exp(params(0));
    double lg = exp(params(1));
    double sigma2 = exp(params(2));
    double m = exp(params(3));
    double tau2 = exp(params(4));
    
    VectorXd mu = params.segment(Nh, Ng);
    VectorXd x = params.segment(Nh+Ng, Ng*Nt);
    
    PriorLengthscale(lt, 0);
    PriorLengthscale(lg, 1);
    PriorLengthscale(m,  3);
    
	//Value += whatever
	//Gradient[i] += etc
}

void Liklihood::PriorLengthscale(double& lengthscale, int& param_index)
{
    // Implements an InverseGamma(1,2) prior for the lengthscales
    
    // We use this three times - division is evil
    double two_over_lengthscale = 2.0/lengthscale;
    
    // lnA = np.log(2.0)-2.0*np.log(l)-2.0/l
    Value += log(two_over_lengthscale/lengthscale) - two_over_lengthscale;
    
    // dlnAdl = 2.0*(1.0-l)/l/l, return dlnAdlnl = l*dlnAdl
    Gradient[param_index] += two_over_lengthscale*(1.0-lengthscale);
}


def log_prior_lengthscale(l):
    
    # InvGamma(1,2)
    
    lnA = np.log(2.0)-2.0*np.log(l)-2.0/l
    
    dlnAdl = 2.0*(1.0-l)/l/l
    
    return lnA, dlnAdl