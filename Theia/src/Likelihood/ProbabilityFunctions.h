#pragma once
#include "../Main/GlobalVariables.h"
#include "../Main/GlobalConstants.h"
#include "DataStorage.h"
#include "MiscFunctions.h"
#include <cmath>



/*!
 * Called by LogLikelihood::PoissonContribution(). Uses a truncated convolution to calculate p(K =k | probs). The implementation is heavily based on the one \verbatim embed:rst:inline `we stole from github <https://github.com/biscarri1/convpoibin/blob/master/src/convpoibin.c>`_ \endverbatim 
 * \param probs The vector of probabilities to calculate the Poisson Binomial on. This vector has length #NumberLargerThanMaxObservations, to prevent continual re-initialisation
 * \param probslen The number of elements of ``probs`` which are initialised. Referred to as ``n`` in our theory work
 * \param result A reference to a matrix / vector-of-vectors used to store the output result in-place. The value of p(K = k| probs) is at element ``[n-1][k]``. The rest of the object is populated with intermediary results needed to calculate the derivative of the result w.r.t. the probabilities.
*/
void  poisson_binomial_pmf_forward(std::vector<double> &  probs, int probslen, std::vector<std::vector<double>> & result);

/*!
 * Called by LogLikelihood::PoissonContribution(). Acts almost identically to poisson_binomial_pmf_forward(), but the convolution happens in reverse order, therefore populating ``result`` with a different set of intermediary results.
 * \param probs The vector of probabilities to calculate the Poisson Binomial on. This vector has length #NumberLargerThanMaxObservations, to prevent continual re-initialisation
 * \param probslen The number of elements of ``probs`` which are initialised. Referred to as ``n`` in our theory work
 * \param result A reference to a matrix / vector-of-vectors used to store the output result in-place. The value of p(K = k| probs) is at element ``[k][n-1]``. The rest of the object is populated with intermediary results needed to calculate the derivative of the result w.r.t. the probabilities.
*/
void  poisson_binomial_pmf_backward(std::vector<double> &  probs, int probslen, std::vector<std::vector<double>> & result);

/*!
 * Uses the output of poisson_binomial_pmf_forward() and poisson_binomial_pmf_forward() to calculate the correct contribution to the derivative of the probabiity (needed for the gradient descent). 
 * \param m An offset parameter either #PipelineMinVisists -1 , k-1 or k, depending on the properties of the star which have been observed.
 * \param pmf_forward the output of poisson_binomial_pmf_forward()
 * \param pmf_backward the output of poisson_binomial_pmf_backward()
 * \param result The vector into which the results are inserted in-place
*/
void  poisson_binomial_subpmf(int m, int probslen, std::vector<std::vector<double>> & pmf_forward, std::vector<std::vector<double>> & pmf_backward, std::vector<double> & result);

/*!
 * Called by LogLikelihood::ExactPoissonContribution(). Almost identical to poisson_binomial_pmf_forward(), but operates in log-space and has a far higher level of precision, and crucially is **much slower**. This should only be called when the results of the original function are in doubt (p ~1 or p~0, or result[x][y] --> infty) 
  *\param probs The vector of probabilities to calculate the Poisson Binomial on. This vector has length #NumberLargerThanMaxObservations, to prevent continual re-initialisation
 * \param probslen The number of elements of ``probs`` which are initialised. Referred to as ``n`` in our theory work
 * \param result A reference to a matrix / vector-of-vectors used to store the output result in-place. The value of p(K = k| probs) is at element ``[n-1][k]``. The rest of the object is populated with intermediary results needed to calculate the derivative of the result w.r.t. the probabilities.
*/
void  poisson_binomial_lpmf_forward(std::vector<double> & probs, int probslen, std::vector<std::vector<double>> & result);


/*!
 * Called by LogLikelihood::ExactPoissonContribution(). Acts almost identically to poisson_binomial_lpmf_forward(), but the convolution happens in reverse order, therefore populating ``result`` with a different set of intermediary results. Again, note that this is considerably slower than poisson_binomial_pmf_backward().
 * \param probs The vector of probabilities to calculate the Poisson Binomial on. This vector has length #NumberLargerThanMaxObservations, to prevent continual re-initialisation
 * \param probslen The number of elements of ``probs`` which are initialised. Referred to as ``n`` in our theory work
 * \param result A reference to a matrix / vector-of-vectors used to store the output result in-place. The value of p(K = k| probs) is at element ``[k][n-1]``. The rest of the object is populated with intermediary results needed to calculate the derivative of the result w.r.t. the probabilities.
*/
void  poisson_binomial_lpmf_backward(std::vector<double> & probs, int probslen, std::vector<std::vector<double>> & result);

/*!
 * Uses the output of poisson_binomial_lpmf_forward() and poisson_binomial_lpmf_forward() to calculate the correct (exact) contribution to the derivative of the probabiity (needed for the gradient descent). 
 * \param m An offset parameter either #PipelineMinVisists -1 , k-1 or k, depending on the properties of the star which have been observed.
 * \param pmf_forward the output of poisson_binomial_pmf_forward()
 * \param pmf_backward the output of poisson_binomial_pmf_backward()
 * \param result The vector into which the results are inserted in-place
*/
void  poisson_binomial_sublpmf(int m, int probslen, std::vector<std::vector<double>> & lpmf_forward, std::vector<std::vector<double>> & lpmf_backward, std::vector<double> & result);

/*!
 * Completely different implementation of the probabilitity model. Allows for variance away from the strictness of the PoissonBinomial (whilst still being an approximation). This is good for us as although k is from Gaia, n is estimated from ScanningLaw stuff, so might be horribly wrong: hence the inclulsion of the VariancePopulation. 
 * \param k The number of detections within the Gaia population for the target star
 * \param probslen The number of visitation of the star (called ``n`` in theory work)
 * \param data A LikelihoodData object, including an initialised LikelihoodData::p vector.
 * \returns The probability p(K = k | p), and populates the gradient vectors LikelihoodData::dfdp_constantN, LikelihoodData::dfdN_constantP and LikelihoodData::hypergradient
*/
double poisson_binomial_normal_lpmf(int k, int probslen, LikelihoodData & data);


