from numba import njit
import numpy as np

# This calculates the full pmf and stores it in result.
# It is much, much quicker to pre-create result and then re-use it!



from math import log, log1p, exp
@njit
def log_sum_exp(a,b):
    if a > b:
        return a + log1p(exp(b-a))
    else:
        return b + log1p(exp(a-b))

@njit
def poisson_binomial_lpmf_forward(probs,probslen,result):
    
    result[0,0] = log1p(-probs[0])
    result[0,1] = log(probs[0])
    
    oldlen = 2
    for i in range(1,probslen):

        # set signal
        log_p = log(probs[i])
        log_q = log1p(-probs[i])

        # initialize result and calculate the two edge cases
        result[i,0] = result[i-1,0] + log_q
        result[i,oldlen] = result[i-1,oldlen-1] + log_p #p * result[i-1,oldlen-1]

        # calculate the interior cases
        for j in range(1,oldlen):
            result[i,j] = log_sum_exp(result[i-1,j-1]+log_p,result[i-1,j]+log_q)
      
        oldlen += 1
        
@njit
def poisson_binomial_lpmf_backward(probs,probslen,result):
    
    result[probslen-1,0] = log1p(-probs[probslen-1])
    result[probslen-1,1] = log(probs[probslen-1])
    
    oldlen = 2
    for i in range(probslen-2,-1,-1):

        # set signal
        log_p = log(probs[i])
        log_q = log1p(-probs[i])

        # initialize result and calculate the two edge cases
        result[i,0] = result[i+1,0] + log_q
        result[i,oldlen] = result[i+1,oldlen-1] + log_p #p * result[i-1,oldlen-1]

        # calculate the interior cases
        for j in range(1,oldlen):
            result[i,j] = log_sum_exp(result[i+1,j-1]+log_p,result[i+1,j]+log_q)
      
        oldlen += 1
        
@njit
def poisson_binomial_sublpmf(m,probs,probslen,lpmf_forward,lpmf_backward,result):
    conv = 0
    result[0] = lpmf_backward[1,m]
    result[probslen-1] = lpmf_forward[probslen-2,m]
    for i in range(1,probslen-1):
        conv = -99999 #pmf_forward[i-1,0]*pmf_backward[i+1,m]
        lower = max(0,m-probslen+i+1)
        upper = min(m,i)+1
        for j in range(lower,upper):
            conv = log_sum_exp(conv,lpmf_forward[i-1,j]+lpmf_backward[i+1,m-j])
        result[i] = conv
        
@njit
def poisson_binomial_pmf(probs,probslen,result):

    result[0] = 1.0-probs[0]
    result[1] = probs[0]
    signal_0,signal_1 = 0.0,0.0
    
    
    oldlen = 2
    for i in range(1,probslen):

        # set signal
        signal_0 = probs[i]
        signal_1 = 1.0-probs[i]

        # initialize result and calculate the two edge cases
        result[oldlen] = signal_0 * result[oldlen-1]

        t = result[0]
        result[0] = signal_1*t

        # calculate the interior cases
        for j in range(1,oldlen):
            tmp=result[j]
            result[j] = signal_0 * t + signal_1 * result[j]
            t=tmp
      
        oldlen += 1
    
@njit
def poisson_binomial_pmf_direct(probs):

    probslen = probs.size
    result = np.zeros(probslen+1)
    poisson_binomial_pmf(probs,probslen,result)
    return result

@njit
def poisson_binomial_pmf_forward(probs,probslen,result):

    result[0,0] = 1.0-probs[0]
    result[0,1] = probs[0]
    
    oldlen = 2
    for i in range(1,probslen):

        # set signal
        p = probs[i]
        q = 1.0-probs[i]

        # initialize result and calculate the two edge cases
        result[i,0] = q * result[i-1,0]
        result[i,oldlen] = p * result[i-1,oldlen-1]

        # calculate the interior cases
        for j in range(1,oldlen):
            result[i,j] = p * result[i-1,j-1] + q * result[i-1,j]
      
        oldlen += 1
        
@njit
def poisson_binomial_pmf_backward(probs,probslen,result):

    result[probslen-1,0] = 1.0-probs[probslen-1]
    result[probslen-1,1] = probs[probslen-1]
    
    oldlen = 2
    for i in range(probslen-2,-1,-1):
        # set signal
        p = probs[i]
        q = 1.0-probs[i]

        # initialize result and calculate the two edge cases
        result[i,0] = q * result[i+1,0]
        result[i,oldlen] = p * result[i+1,oldlen-1]

        # calculate the interior cases
        for j in range(1,oldlen):
            result[i,j] = p * result[i+1,j-1] + q * result[i+1,j]
      
        oldlen += 1
        

@njit
def poisson_binomial_subpmf(m,probs,probslen,pmf_forward,pmf_backward,result):
    conv = 0
    result[0] = pmf_backward[1,m]
    result[probslen-1] = pmf_forward[probslen-2,m]
    for i in range(1,probslen-1):
        conv = 0 #pmf_forward[i-1,0]*pmf_backward[i+1,m]
        for j in range(max(0,m-probslen+i+1),min(m,i)+1):
            conv += pmf_forward[i-1,j]*pmf_backward[i+1,m-j]
        result[i] = conv
        
        
@njit
def poisson_binomial_selection(probs,probslen,result,k,c):
    
    if k < c:
        return 0

    result[0] = 1.0-probs[0]
    result[1] = probs[0]
    signal_0,signal_1 = 0.0,0.0
    
    
    oldlen = 2
    for i in range(1,probslen):

        # set signal
        signal_0 = probs[i]
        signal_1 = 1.0-probs[i]

        # initialize result and calculate the two edge cases
        result[oldlen] = signal_0 * result[oldlen-1]

        t = result[0]
        result[0] = signal_1*t

        # calculate the interior cases
        for j in range(1,min(oldlen,k+1)):
            tmp=result[j]
            result[j] = signal_0 * t + signal_1 * result[j]
            t=tmp
      
        oldlen += 1
        
    correction = 0
    for i in range(c):
        correction += result[i]

    return result[k]/(1.0-correction)

@njit
def poisson_binomial_log_likelihood_truncated(k,c,n,probs,gradient,pmf,subpmf):
    
    # Assumes that c <= k <= n
    # k is the number of measurements reported by Gaia
    # c is the minimum number of measurements for a star to be reported by Gaia
    # n is the predicted number of visits by Gaia of that star
    # p is the list of measurement probabilities at each visit
    # log_likelihood is where you want the log likelihood to be stored (we can have this be an output of the function instead, if you prefer)
    # gradient is where you want the gradient to be stored (same as above)
    # pmf and subpmf are used for intermediate storage to avoid excessive memory allocation. Set these to be 1000 long and use them for every star.
    
    # Compute the pmf and log_likelihood
    poisson_binomial_pmf(probs,n,pmf)
    
    likelihood = pmf[k]
    correction = 1.0 - np.sum(pmf[0:c])
    log_likelihood = np.log(likelihood) - np.log(correction)
    
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
        
        gradient[i] += (gradient_first_term*subpmf[k-1]-gradient_second_term*subpmf[k])/likelihood - subpmf[c-1]/correction
    return log_likelihood