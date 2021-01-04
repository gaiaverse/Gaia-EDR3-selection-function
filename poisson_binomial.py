from numba import njit
import numpy as np

# This calculates the full pmf and stores it in result.
# It is much, much quicker to pre-create result and then re-use it!

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

    return result

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