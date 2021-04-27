import emcee
from numba import njit
from scipy import special
#import numba_special
import numpy as np
import h5py
import tqdm

from numba.extending import get_cython_function_address
from numba import vectorize
import ctypes

addr = get_cython_function_address("scipy.special.cython_special", "__pyx_fuse_1log_ndtr")
functype = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
log_ndtr_fn = functype(addr)


@vectorize('float64(float64)')
def vec_log_ndtr(x):
    return log_ndtr_fn(x)

@njit
def log_ndtr_in_njit(x):
    return vec_log_ndtr(x)


@njit
def log_prior(theta):
    alpha,beta = theta
    if -100.0 < alpha < -1.0 and 0.001 < beta < 100.0:
        return 0.0
    return -np.inf

@njit
def log_likelihood(theta, g):
    alpha, beta = theta
    y = (21.0-g)/np.sqrt(0.01**2+np.exp(2.0*(alpha+beta*(g-20.5))))
    return log_ndtr_in_njit(y).sum()

@njit
def log_probability(theta, g):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, g)

n_stars = 3*17*6791 # number of stars in each block
n_block = 5231 # number of blocks

# Compute percentiles of sigma5dmax
n_blocks = 1
phot_g_mean_mag = np.zeros((n_blocks,n_stars))
block_choice = np.random.choice(np.arange(n_block),n_blocks,replace=False).astype(np.int)


with h5py.File('./gaiaedr3.h5', 'r') as f:
    
    for idx_order,idx_block in tqdm.tqdm(enumerate(block_choice)):

        # Load in data
        f['phot_g_mean_mag'].read_direct(phot_g_mean_mag[idx_order],np.s_[idx_block*n_stars:(idx_block+1)*n_stars])

phot_g_mean_mag = phot_g_mean_mag.ravel()
_keep = np.where((np.isnan(phot_g_mean_mag) == False) & (phot_g_mean_mag > 20.5))
phot_g_mean_mag = phot_g_mean_mag[_keep][:1000]
print('kept',_keep[0].size)

pos = np.array([-40.0,1.0]) + 1e-4 * np.random.randn(32, 2)
nwalkers, ndim = pos.shape

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(phot_g_mean_mag,))
sampler.run_mcmc(pos, 5000, progress=True);

flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
print(np.percentile(flat_samples,[16,50,84],axis=0))
