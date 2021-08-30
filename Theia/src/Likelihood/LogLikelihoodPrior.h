#pragma once
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#include "LogLikelihood.h"
#include "../Main/EfficiencyVector.h"
using Eigen::VectorXd;
using namespace Eigen;


/*!
 * A holder struct for the value and derivative of functions of a single variable. Useful in non-correlated priors
*/
struct F_dF
{
	//! Functional value
	double F;
	//! Functional derivative
	double dF;
	
	//! Default constructor: zero-initialised values 
	F_dF()
	{
		F = 0;
		dF = 0;
	}
	//! In-place constructor function. Sets the passed values to their stored counterparts.
	F_dF(double f, double df) : F(f), dF(df){};
};

/*!
 * Implementation of (the logarithm of) the Gaussian-normal prior N(X = x | mu, sigma) 
 * \param x The parameter value
 * \param mu The mean value of the prior on the parameter
 * \param sigma The standard deviation of the prior on the parameter
*/
F_dF inline Normal(double x, double mu, double sigma)
{
	double d = (x - mu);
	double s2 = sigma*sigma;
	double v = -0.5 * d*d/(s2);
	double dv = - d/s2;
	
	return F_dF(v,dv);
};

/*!
 * Implementation of (logarithm of) the Student-T prior S(X = x | mu, nu) 
 * \param x The parameter value
 * \param mu The mean value of the prior on the parameter
 * \param nu The student-T parameter
*/
F_dF inline StudentT(double x, double mu, double nu)
{
	double div = 2;
	double d = (x - mu)/div;
	double v = - (nu + 1)/2 * log(1 + d*d/nu);
	double dv = - (nu + 1) * d/ ( nu * (1 + d*d/nu) * div);
	
	return F_dF(v,dv); 
};

/*!
 * Implementation of (the logarithm of) the Beta-Distribution prior B(P = p | alpha, beta) where p is the probability associated with the temporal efficiency parameter x_t
 * \param x The temporal efficiency parameter x_t (given by p_t = sigmoid(double x_t) )
 * \param alpha The Beta distribution  alpha-parameter
 * \param sigma The Beta distribution beta-parameter
*/
F_dF inline TemporalBetaPrior(double x)
{
	double v = gapPriorAlpha * x - (gapPriorBeta + gapPriorAlpha)* log(1.0 + exp(x));
	double dv = (gapPriorAlpha - gapPriorBeta* exp(x))/(1 + exp(x));
	return F_dF(v,dv);
};



/*!
 * A subclass of the LogLikelihood. This version is run by the Root process from within the DescentFunctor object. This object is almost identical to the base class, with the exception of its ability to call Prior functions. 
*/
class LogLikelihoodPrior : public LogLikelihood
{
	public:
		LogLikelihoodPrior(const std::vector<std::vector<Star>> & data): LogLikelihood(data){};
	
		
	     
	    /*! 
	     *Executes the Prior in transform space. Calls the TemporalBetaPrior on all x_t which lie within one of the pre-catalogued gaps, and hence ensures that our variance model expands to capture those gaps, rather than attempting to smooth over them.
	     \param x The current efficiency vector 
	     * \param effectiveBatches the current number of minibatches
	     * \returns The value of the prior. The associated gradients are loaded into #EfficiencyVector::TransformedGradient
	    */
	    double TransformPrior(EfficiencyVector & x, int effectiveBatches);
	    
	    /*!
	     * Executes the Prior in Raw space. By design, the priors here are very simple (either Normal() or StudentT()) with zero mean and unit variance. The transform enforces any correlations or smoothing.
	     * \param x The current efficiency vector (which has had EfficiencyVector::BackwardTransform() already called on it)
	     * \param effectiveBatches the current number of minibatches
	     * \returns The value of the prior. The associated gradients are loaded into #EfficiencyVector::RawGradient
	    */
	    double RawPrior(EfficiencyVector & x, int effectiveBatches);
		

};
