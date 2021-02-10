import emcee
from numba import njit
import scipy.special
import numba_special


@njit
def log_prior(theta):
    alpha,beta = theta
    if -10.0 < alpha < -1.0 and 0.1 < beta < 10.0:
        return 0.0
    return -np.inf

@njit
def log_likelihood(theta, g):
    alpha, beta = theta
    y = (21.0-G)/np.sqrt(0.01**2+np.exp(2.0*(alpha+beta*G)))
    return special.log_ndtr(y)

@njit
def log_probability(theta, g):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, g)

n_stars = 3*17*6791 # number of stars in each block
n_block = 5231 # number of blocks

# Compute percentiles of sigma5dmax
keys = ['phot_g_mean_mag']
n_blocks = 10
phot_g_mean_mag = np.zeros((n_blocks,n_stars))
block_choice = np.random.choice(np.arange(n_block),n_blocks,replace=False).astype(np.int)


with h5py.File('./gaiaedr3.h5', 'r') as f:
    
    for idx_order,idx_block in tqdm.tqdm(enumerate(block_choice)):

        # Load in data
        f['phot_g_mean_mag'].read_direct(phot_g_mean_mag[idx_order],np.s_[idx_block*n_stars:(idx_block+1)*n_stars])

phot_g_mean_mag = phot_g_mean_mag.ravel()