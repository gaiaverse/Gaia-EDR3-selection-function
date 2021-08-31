#pragma once
#include "../Main/GlobalVariables.h"
#include "../Main/GlobalConstants.h"

/*! 
 * Implements an expit sigmoid. Some cleverness was used to avoid numerical overflows for exp(large numbers).
 * \param x Value to be expit-ed
 * \returns  \verbatim embed:rst:inline :math:`\sigma(x) = \frac{1}{1 + \exp(-x)}` \endverbatim
*/
double  sigmoid(double x);

/*!
 * A truncated exponential (apparently elu is a mahcine learning term?). For x > #elu_transitionPoint, returns exp(-x), otherwise returns the linear function which makes elu(x) both continuous and smooth. 
 * \param x: Value to be elu-ed
 * \returns \verbatim embed:rst:inline :math:`\text{elu}(x) =\begin{cases} \exp(-x) ~~~& x > \rho \\ (1 + \rho - x) \exp(-\rho) & x \leq \rho  \end{cases}` \endverbatim
*/
double  elu(double x);

/*!
 * Calculate the gradient of elu() evaluated at x. 
 * \param x The value at which gradient is to be evaluated at
 * \param elu_x the functional value of elu() at the chosen point. Can avoid having to recalulate exp(x) if possible
 * \returns The requisite gradient
*/
double  elu_grad(double x, double elu_x);

/*!
 * Given two (potentially very large) numbers a = log(x) ,b = log(y) we often want to know log(x+y). To compute this we need to take exponentials and then add/subtract them: even if the answer *should* be a finite number, the intermediary values can often lead to numerical overflows. This function calculates it in a much safer way. 
 * \param a The logarithm of a large number
 * \param b The logarithm of another large number
 * \returns  \verbatim embed:rst:inline :math:`\log\left( e^a + e^b\right)` \endverbatim, calculated such that the intermediary results do not overflow
*/
double  log_add_exp(double a, double b);


/*!
 * Some numbers we'd rather not calculate again and again if possible
*/
const double one_over_root2 = 1.0/sqrt(2.0);
const double one_over_root2pi = 1.0/sqrt(2.0*M_PI);
const double root2_over_pi = sqrt(2.0/M_PI);

