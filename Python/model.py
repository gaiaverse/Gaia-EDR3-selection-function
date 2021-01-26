from numba import njit
import numpy as np
            
@njit
def poisson_binomial_log_likelihood_truncated(k,c,n,probs,log_likelihood,gradient,pmf,subpmf):
    
    # Assumes that c <= k <= n
    # k is the number of measurements reported by Gaia
    # c is the minimum number of measurements for a star to be reported by Gaia
    # n is the predicted number of visits by Gaia of that star
    # p is the list of measurement probabilities at each visit
    # log_likelihood is where you want the log likelihood to be stored (we can have this be an output of the function instead, if you prefer)
    # gradient is where you want the gradient to be stored (same as above)
    # pmf and subpmf are used for intermediate storage to avoid excessive memory allocation. Set these to be 1000 long and use them for every star.
    
    # Compute the pmf and log_likelihood
    direct_convolution_local(probs,n,pmf)
    
    likelihood = pmf[k]
    correction = 1.0 - np.sum(pmf[0:c])
    log_likelihood[:] = np.log(likelihood) - np.log(correction)
    
    # Branching to set terms in the loop
    if k == 0:
        gradient_first_term,gradient_second_term=0.0,1.0
    elif k == n:
        gradient_first_term,gradient_second_term=1.0,0.0
    else:
        gradient_first_term,gradient_second_term=1.0,1.0
    
    # Loop over the p's calculating the gradient
    for i in range(n):
        
        p = probs[i]
        oneoveroneminusp = 1.0/(1.0-p)
        
        subpmf[0] = pmf[0]*oneoveroneminusp
        for j in range(1,n):
            subpmf[j] = (pmf[j]-subpmf[j-1]*p)*oneoveroneminusp
        subpmf[n-1] = pmf[n]/p
        
        gradient[i] = (gradient_first_term*subpmf[k-1]-gradient_second_term*subpmf[k])/likelihood - subpmf[c-1]/correction
        
p = np.array([0.3,0.6,0.2,0.1,0.9,0.5,0.3,0.2,0.8])
n = p.size
k = 6
c = 5
Ns = 500
subpmf = np.zeros(Ns)
gradient = np.zeros(Ns)
log_likelihood = np.zeros(1)
poisson_binomial_log_likelihood_truncated(k,c,n,p,log_likelihood,gradient,pmf,subpmf)

dp = 1e-8
print('Predicted log-likelihood at p   ',log_likelihood[0])
print('Predicted log-likelihood at p+dp',log_likelihood[0]+dp*gradient[3])

# Should print:
# Predicted log-likelihood at p    -1.3512295610051956
# Predicted log-likelihood at p+dp -1.3512295548589086