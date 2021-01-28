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
    PriorX(x, mu, lt, lm, sigma2);
    
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
    // Implements the Gaussian Process prior on mu
    
    Matrix<double, Ng, Ng> K;
    Matrix<double, Ng, Ng> dK_dm;
    double magnitude_distance;
    
    // Build covariance matrix
    for (int i = 0; i < Ng; i++) {
        for (int j = 0; j < Ng; j++) {
            magnitude_distance = pow(magnitude(i)-magnitude(j),2);
            K(i,j) = K(j,i) = exp(-magnitude_distance/(2.0*m*m));
            dK_dm(i,j) = dK_dm(j,i) = K(i,j)*magnitude_distance/(m*m*m);
        }
    }
    
    // Householder decomposition (with pivoting!)
    // It's possible that this is all broken.
    Eigen::fullPivHouseholderQr()< Matrix<double, Ng, Ng> > decomp(K);
    
    // Compute quantities we will need later
    Vector<double, Ng> invKmu = decomp.solve(mu);
    Matrix<double, Ng, Ng> J = decomp.solve(dK_dm);
    double muinvJmu = invKmu.dot(mu);
    
    // lnQ = +0.5*np.linalg.slogdet(J_inv)[1]-0.5*np.dot(mu.T,J_inv_mu) - (M/2)*np.log(2.0*np.pi)
    Value += -0.5*Ng*log(2.0*M_PI*tau2) -0.5*decomp.logAbsDeterminant() -0.5*muinvKmu/tau2;
    
    // dlnQdm = -0.5*np.trace(np.dot(J_inv,dJdm))+0.5*np.dot(J_inv_mu.T,np.dot(dJdm,J_inv_mu))
    Gradient[3] += m*(-0.5*J.trace() + 0.5*invKmu.adjoint()*dJdm*invKmu/tau2);
    
    // dlnQ_dlntau2, correcting for log factor
    Gradient[4] += -0.5*Ng + 0.5*muinvKmu/tau2;
    
    Gradient.segment(Nh,Ng).array() -= invKmu/tau2;
}

void Liklihood::PriorX(Eigen::VectorXd& x, Eigen::VectorXd& mu, double& lt, double& lm, double& sigma2)
{
    // Implements the Gauss-Markov prior on x
    
    // Useful shorthands
    double u = exp(-1.0/lt);
    double u2 = u*u;
    double oneplusu2 = 1.0 + u2;
    double oneoveroneminusu2 = 1.0 / ( 1.0 - u2 );
    
    // Reshape to form Y
    // Y = x.reshape((Ng,Nt)) - mu.reshape((Ng,1))
    Matrix<double, Ng, Nt> Y;
    for (int i = 0; i < Ng; i++) {
        Y.row(i) = x.segment(i*Nt,Nt) - mu[i];
    }
    
    // Create quantities relating to Kg
    Matrix<double, Ng, Ng> Kg;
    Matrix<double, Ng, Ng> dKg_dlg;
    double magnitude_distance;
    
    // Build covariance matrix
    for (int i = 0; i < Ng; i++) {
        for (int j = 0; j < Ng; j++) {
            magnitude_distance = pow(magnitude(i)-magnitude(j),2);
            Kg(i,j) = Kg(j,i) = exp(-magnitude_distance/(2.0*lg*lg));
            dKg_dlg(i,j) = dKg_dlg(j,i) = Kg(i,j)*magnitude_distance/(lg*lg*lg);
        }
    }
    
    // Householder decomposition (with pivoting!)
    // It's possible that this is all broken.
    Eigen::fullPivHouseholderQr()< Matrix<double, Ng, Ng> > decomp(Kg);
    
    // Compute quantities we will need later
    Matrix<double, Ng, Ng> Jg = decomp.solve(dKgdlg);
    Matrix<double, Ng, Ng> invKgY = decomp.solve(Y);
    double logdetKg = decomp.logAbsDeterminant();
    double TrJg =  Jg.trace();
    
    double logdetinvKt = (Nt-1.0)*log( oneoveroneminusu2 );
    double TrJt = -2.0*(Nt-1.0)*u2*oneoveroneminusu2/(lt*lt);
    

        
    // Compute invKgYinvKt
    Matrix<double, Ng, Nt> invKgYinvKt;
    for (int ig = 0; ig < Ng; ig++) {
        invKgYinvKt[ig,0] = ( invKgY(ig,0) - u * invKgY(ig,1) )*oneoveroneminusu2;
        invKgYinvKt[ig,-1] = ( invKgY(ig,Nt-1) - u * invKgY(ig,Nt-2) )*oneoveroneminusu2;
        for (int it = 1; it < Nt-1; it++) {
            invKgYinvKt(ig,it) = ( oneplusu2 * invKgY(ig,it) - u * ( invKgY(ig,it-1) + invKgY(ig,it+1) ) )*oneoveroneminusu2;
        }
    }
    
    double Y_invKgYinvKt = (Y.array()*invKgYinvKt.array()).sum();
    
    double JgTY_invKgYinvKt = ((Jg.transpose()*Y).array()*invKgYinvKt.array()).sum();
    
    
    // Compute YJt
    Matrix<double, Ng, Nt> YJt;
    int M = 10+ceil(-lt*log(1e-16));
    Vector<double, M> power_u;
    power_u[0] = u;
    for (int i = 1; i < M; i++) {
        power_u[i] = u*power_u[i-1];
    }
    
    double res;
    for (int i = 0; i < Ng; i++) {
        for (int l = 0; l < Nt; l++) {
            
            res =  - oneplusu2*Y[i,l]*oneoveroneminusu2;
            
            if l < M:
                res += u2*Y[i,0]*power_u[l]*oneoveroneminusu2;
                
            if Nt-l-1 < M:
                res += u2*Y[i,-1]*power_u[Nt-l-1]*oneoveroneminusu2;
            
            for (int j = max(0,l-M); j < min(Nt,l+M); j++) {
                res += Y(i,j)*power_u[abs(l-j)];
            }
            YJt[i,l] = res/(lt*lt);
        }
    }
                
    Matrix<double, Ng, Nt> YJt;        
    double YJt_invKgYinvKt = (YJt.array()*invKgYinvKt.array()).sum();
                
    // lnP = -Ng*Nt*np.log(2.0*np.pi*sigma2)/2.0 + Ng*logdetinvKt/2.0 - Nt*logdetKg/2.0 - Y_invKgYinvKt/2.0/sigma2
    Value += -Ng*Nt*log(2.0*M_PI*sigma2)/2.0 + Ng*logdetinvKt/2.0 - Nt*logdetKg/2.0 - Y_invKgYinvKt/(2.0*sigma2);
    
    //Matrix<double, Ng, Nt> dlnP_dX = -invKgYinvKt/sigma2;
    //dlnP_dx = dlnP_dX.ravel()
    // dlnP_dmu = -dlnP_dX.sum(axis=1)
    for (int ig = 0; ig < Ng; ig++) {
        Gradient.segment(Nh+Ng+ig*Nt,Nt).array() += -invKgYinvKt.row(ig)/sigma2;
        Gradient[Nh+ig] += -Gradient.segment(Nh+Ng+ig*Nt,Nt).sum();
    }
    
    //dlnP_dlt = -Ng*TrJt/2.0 + YJt_invKgYinvKt/2.0/sigma2
    Gradient[0] += lt*(-Ng*TrJt/2.0 + YJt_invKgYinvKt/(2.0*sigma2));
    
    //dlnP_dlg = -Nt*TrJg/2.0 + JgTY_invKgYinvKt/2.0/sigma2
    Gradient[1] += lg*(-Nt*TrJg/2.0 + JgTY_invKgYinvKt/(2.0*sigma2));
            
    //dlnP_dsigma2 = -Ng*Nt/2.0/sigma2 + Y_invKgYinvKt/2.0/sigma2/sigma2
    Gradient[2] += -Ng*Nt/2.0 + Y_invKgYinvKt/2.0/sigma2;

}


            
            
