""" Open the HDF5 file and histogram how many stars are in each magnitude bin. """

import numpy as np
import h5py
import tqdm
import pickle

n_stars = 3*17*6791 # number of stars in each block
n_block = 1 # number of blocks in each group
n_group = 5231

# Compute percentiles of sigma5dmax
keys = ['phot_g_mean_mag']
output_box = {}
output_box['bins'] = np.array([0.0, 5.0, 16.0, 17.0, 18.0, 18.2, 18.4, 18.6, 18.8, 19.0, 19.1, 19.2, 19.3, 19.4, 19.5, 19.6, 19.7, 19.8, 19.9, 20.0, 20.1, 20.2, 20.3, 20.4, 20.5, 20.6, 20.7, 20.8, 20.9, 21.0, 21.1, 21.2, 21.3, 21.4, 21.5, 25.0])
output_box['counts'] = np.zeros((output_box['bins'].size-1))

with h5py.File('./gaiaedr3.h5', 'r') as f:
    for idx_group in tqdm.tqdm(range(n_group)):

        # Load in data

        box = {k:np.zeros(n_block*n_stars) for k in keys}
        for _k in keys:
            f[_k].read_direct(box[_k],np.s_[idx_group*n_block*n_stars:(idx_group+1)*n_block*n_stars])
                
        # Add to histogram
        output_box['counts'] += np.histogram(box['phot_g_mean_mag'],bins=output_box['bins'])[0]

# Output
with open('./magnitude_bin_counts.p', 'wb') as fp:
    pickle.dump(output_box,fp)