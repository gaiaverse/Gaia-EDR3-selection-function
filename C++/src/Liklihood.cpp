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
    
    // Apply the priors on the hyper-hyper-parameters
    PriorLengthscale(lt,  0);
    PriorLengthscale(lg,  1);
    PriorVariance(sigma2, 2);
    PriorLengthscale(m,   3);
    PriorVariance(tau2,   4);
    
    // Apply the prior on the hyper-parameters
    PriorMu(mu, m, tau2);
    
    // Apply the prior on the parameters
    //PriorX(mu, m, tau2);
    
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

void Liklihood::PriorVariance(double& variance, int& param_index)
{
    // Implements a Gamma(1,1) prior for the variances
    
    // lnS = - sigma2
    Value -= variance;
    
    // dlnSdsigma2 = - 1.0, return dlnSdlnsigma2 = sigma2*dlnSdsigma2
    Gradient[param_index] -= variance;
}

void Liklihood::PriorMu(Eigen::VectorXd& mu, double& m, double& tau2)
{
    // Implements a Gamma(1,1) prior for the variances
    
    Matrix<double, Ng, Ng> J;
    Matrix<double, Ng, Ng> dJdm;
    double magnitude_distance;
    
    // Build covariance matrix
    for (int i = 0; i < Ng; i++) {
        for (int j = 0; j < Ng; j++) {
            magnitude_distance = pow(magnitude(i)-magnitude(j),2);
            J(i,j) = J(j,i) = exp(-magnitude_distance/(2.0*m*m));
            dJdm(i,j) = dJdm(j,i) = J(i,j)*magnitude_distance/(m*m*m);
        }
    }
    
    // Householder decomposition (with pivoting!)
    // It's possible that this is all broken.
    Eigen::fullPivHouseholderQr()< Matrix<double, Ng, Ng> > decomp(J);
    
    // Compute quantities we will need later
    Vector<double, Ng> invJmu = decomp.solve(mu);
    Matrix<double, Ng, Ng> invJdJdm = decomp.solve(dJdm);
    double muinvJmu = invJmu.dot(mu);
    
    // lnQ = +0.5*np.linalg.slogdet(J_inv)[1]-0.5*np.dot(mu.T,J_inv_mu) - (M/2)*np.log(2.0*np.pi)
    Value += -0.5*Ng*log(2.0*M_PI*tau2) -0.5*decomp.logAbsDeterminant() -0.5*muinvJmu/tau2;
    
    // dlnQdm = -0.5*np.trace(np.dot(J_inv,dJdm))+0.5*np.dot(J_inv_mu.T,np.dot(dJdm,J_inv_mu))
    Gradient[3] += m*(-0.5*invJdJdm.trace() + 0.5*invJmu.adjoint()*dJdm*invJmu/tau2);
    
    // dlnQ_dlntau2, correcting for log factor
    Gradient[4] += -0.5*Ng + 0.5*muinvJmu/tau2;
    
    Gradient.segment(Nh,Ng) -= invJmu/tau2;
}

/*
def log_prior_mu(mu,m,tau2):
    
    J = tau2*np.exp(-G_squared/2.0/m**2)
    J_inv = np.linalg.pinv(J)
    J_inv_mu = np.dot(J_inv,mu)
    
    lnQ = +0.5*np.linalg.slogdet(J_inv)[1]-0.5*np.dot(mu.T,J_inv_mu) - (M/2)*np.log(2.0*np.pi)
    
    dJdm = J*G2/m**3
    dJdtau2 = J/tau2

    dlnQdm = -0.5*np.trace(np.dot(J_inv,dJdm))+0.5*np.dot(J_inv_mu.T,np.dot(dJdm,J_inv_mu))
    dlnQdtau2 = -0.5*np.trace(np.dot(J_inv,dJdtau2))+0.5*np.dot(J_inv_mu.T,np.dot(dJdtau2,J_inv_mu))
    dlnQdmu = -J_inv_mu
    
    return lnQ, dlnQdmu, dlnQdm, dlnQdtau2
    
*/