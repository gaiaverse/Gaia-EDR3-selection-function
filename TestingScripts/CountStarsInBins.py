""" Open the HDF5 file and histogram how many stars are in each magnitude bin. """

import numpy as np
import h5py
import tqdm
import multiprocessing

N_core = 20
N_sources = 1811709771
N_chunk = 3*17*6791
N_block = 5231
N_maxobs = 256

g_bins = np.arange(1.7,23.05,0.1)
N_g = g_bins.size-1
n_bins = np.arange(-0.5,2*N_maxobs)
N_n = n_bins.size-1

def mp_worker(args):
    
    block_idx = args

    with h5py.File('./gaiaedr3.h5', 'r') as g, h5py.File('./gaiaedr3_times.h5', 'r') as t:
        
        # Load in source_id's and verify that they are the same
        source_id_g = g['source_id'][block_idx*N_chunk:(block_idx+1)*N_chunk]
        source_id_t = t['source_id'][block_idx*N_chunk:(block_idx+1)*N_chunk]
        assert (source_id_g == source_id_t).all()
        
        # Load in data
        phot_g_mean_mag = g['phot_g_mean_mag'][block_idx*N_chunk:(block_idx+1)*N_chunk]
        matched_transits = g['astrometric_matched_transits'][block_idx*N_chunk:(block_idx+1)*N_chunk]
        fov_1_n = t['fov_1_n'][block_idx*N_chunk:(block_idx+1)*N_chunk].astype(np.uint16)
        fov_2_n = t['fov_2_n'][block_idx*N_chunk:(block_idx+1)*N_chunk].astype(np.uint16)
        
    return np.histogramdd((phot_g_mean_mag, fov_1_n + fov_2_n, matched_transits), bins = (g_bins, n_bins, n_bins))[0]

pool = multiprocessing.Pool(N_core)
grid = np.zeros((N_g,N_n,N_n))
for result in tqdm.tqdm(pool.imap(mp_worker, range(N_block)),total=N_block):
    grid += result

with h5py.File('./gaiaedr3_grid.h5', 'w') as f:
    f.create_dataset('grid', data=grid, compression = "lzf", dtype = np.uint64, fletcher32 = False, shuffle = True, scaleoffset=0)
    f.create_dataset('extent', data=np.array([g_bins[0],g_bins[-1],n_bins[0],n_bins[-1],n_bins[0],n_bins[-1]]))