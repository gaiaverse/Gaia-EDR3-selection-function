""" Script to extract healpix maps of the scanning law and the Gaia source density. """

import numpy as np
import h5py
import tqdm
import multiprocessing
import healpy as hp
from fast_histogram import histogram2d
from numba import njit

# First, load in all of the phot_g_mean_mag's and source_id's
n_stars = 3*17*6791 # number of stars in each block
n_block = 5231 #
n_total = n_block*n_stars
nside = 1024
npix = hp.nside2npix(nside)

keys = ['phot_g_mean_mag','source_id']
box = {k:np.zeros(n_total) for k in keys}
box['source_id'] = box['source_id'].astype(np.int64)

@njit
def incrementer(array,i,j):
    for _i,_j in zip(i,j):
        array[_i,_j] += 1

with h5py.File('./gaiaedr3.h5', 'r') as f:
    for idx_block in tqdm.tqdm(range(n_block)):
        
        idx_low, idx_upp = idx_block*n_stars, (idx_block+1)*n_stars
        
        # Load in data
        for _k in keys:
            
            f[_k].read_direct(box[_k][idx_low:idx_upp],np.s_[idx_low:idx_upp])
                
# Compute percentiles
g_percentiles = np.nanpercentile(box['phot_g_mean_mag'],[0,10,20,30,40,50,60,70,80,90,100])

# Compute maps
sky_map = np.zeros((11,npix))
for idx_block in tqdm.tqdm(range(n_block)):
    idx_low, idx_upp = idx_block*n_stars, (idx_block+1)*n_stars
    
    hpx_idx = np.floor(box['source_id'][idx_low:idx_upp]/549755813888).astype(np.int)
    g_idx = np.digitize(box['phot_g_mean_mag'][idx_low:idx_upp],g_percentiles)-1
    incrementer(sky_map,g_idx,hpx_idx)
    #sky_map += np.histogram2d(g_idx,hpx_idx,bins=(g_percentiles,np.arange(-0.5,npix)))[0]

# Load scanning law

with h5py.File('./gaiaedr3_times_healpix_1024_16.h5', 'r') as f:
        
    scanninglaw_map = f['fov_1_n'][:]+f['fov_2_n'][:]

# Output
with h5py.File(f'./powermaps.h5', 'w') as f:
        
    # Create datasets
    f.create_dataset('phot_g_mean_mag_percentiles', data=g_percentiles, dtype=np.float64)
    f.create_dataset('sky_map', data=sky_map, compression="gzip", compression_opts=9, chunks = True, shape = (11,npix,), dtype = np.uint64, fletcher32 = False, shuffle = True, scaleoffset=0)
    f.create_dataset('scanninglaw_map', data=scanninglaw_map, compression="gzip", compression_opts=9, chunks = True, shape = (npix,), dtype = np.uint64, fletcher32 = False, shuffle = True, scaleoffset=0)