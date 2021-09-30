""" For each Gaia chunk, compute n observations. """

import numpy as np
import sys
import pickle
import tqdm
import pandas as pd
from scipy import spatial
import multiprocessing
import os
import time
import h5py
from scipy import interpolate
from math import sqrt, atan2, asin
from numba import njit
            
##### Load in sources and assign to bins
N_core = 5
N_chunk = int(2e5)
N_block = 50
N_sources = N_chunk*N_block
N_maxobs = 256
print('Using',N_core,'cores.')

##### Load in data that all will need
    
# Scanning law
_data = pd.read_csv('./CommandedScanLaw_001.csv')
_columns = ['jd_time', 'ra_fov1', 'dec_fov1', 'ra_fov2', 'dec_fov2']
_keys = ['tcb_at_gaia', 'ra_fov_1', 'dec_fov_1', 'ra_fov_2', 'dec_fov_2']
_box = {}
_order = np.argsort(_data['jd_time'].values)
for j,k in zip(_columns,_keys):
    _box[k] = _data[j].values[_order]
del _data
print('There are',_box['tcb_at_gaia'].size,'time points.')
    
fov_1_xyz = np.stack([np.cos(np.deg2rad(_box['ra_fov_1']))*np.cos(np.deg2rad(_box['dec_fov_1'])), np.sin(np.deg2rad(_box['ra_fov_1']))*np.cos(np.deg2rad(_box['dec_fov_1'])), np.sin(np.deg2rad(_box['dec_fov_1']))]).T
fov_2_xyz = np.stack([np.cos(np.deg2rad(_box['ra_fov_2']))*np.cos(np.deg2rad(_box['dec_fov_2'])), np.sin(np.deg2rad(_box['ra_fov_2']))*np.cos(np.deg2rad(_box['dec_fov_2'])), np.sin(np.deg2rad(_box['dec_fov_2']))]).T
fov_1_tree = spatial.cKDTree(fov_1_xyz)
fov_2_tree = spatial.cKDTree(fov_2_xyz)
fov_times = _box['tcb_at_gaia']
del _box
print('Loaded scanning law.')

# Emphemeris
speed_of_light_AU_per_day = 299792458.0*(86400.0/149597870.700/1e3)
gaia_ephem_data = pd.read_csv('./horizons_results_gaia.txt',skiprows=64)
gaia_ephem_box = {k:gaia_ephem_data[k].values for k in ['JDTDB','X','Y','Z','VX','VY','VZ']}
gaia_ephem_velocity = interpolate.interp1d( gaia_ephem_box['JDTDB']-2455197.5, np.stack([gaia_ephem_box['VX'],gaia_ephem_box['VY'],gaia_ephem_box['VZ']])/speed_of_light_AU_per_day, kind='cubic')
gaia_velocity = gaia_ephem_velocity(fov_times).T
print('Loaded emphemeris data.')

# Compute rotation matrices
_xaxis = fov_1_xyz
_zaxis = np.cross(fov_1_xyz,fov_2_xyz)
_yaxis = -np.cross(_xaxis,_zaxis)
_yaxis /= np.linalg.norm(_yaxis,axis=1)[:,np.newaxis]
_zaxis /= np.linalg.norm(_zaxis,axis=1)[:,np.newaxis]

_uaxis = np.array([1,0,0])
_vaxis = np.array([0,1,0])
_waxis = np.array([0,0,1])

_matrix = np.moveaxis(np.stack([np.dot(_xaxis,_uaxis),np.dot(_xaxis,_vaxis),np.dot(_xaxis,_waxis),
                                np.dot(_yaxis,_uaxis),np.dot(_yaxis,_vaxis),np.dot(_yaxis,_waxis),
                                np.dot(_zaxis,_uaxis),np.dot(_zaxis,_vaxis),np.dot(_zaxis,_waxis)]).reshape((3,3,_xaxis.shape[0])),2,0)
print('Computed rotation matrices.')

def sample_spherical(npoints):
    vec = np.random.randn(3, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    ra = np.arctan2(vec[1],vec[0])
    dec = np.arctan2(vec[2],np.sqrt(vec[0]*vec[0]+vec[1]*vec[1]))
    return ra,dec

def mp_worker(args):
    
    block_idx = args
    
    # Randomly generate stars
    source_id = np.arange(block_idx*N_chunk,(block_idx+1)*N_chunk).astype(int)
    source_ra, source_dec = sample_spherical(N_chunk)
    
    # Compute tree
    source_xyz = np.stack([np.cos(source_ra)*np.cos(source_dec), np.sin(source_ra)*np.cos(source_dec), np.sin(source_dec)]).T
    source_tree = spatial.cKDTree(source_xyz)
    
    ##### Find intersections
    source_N = N_chunk
    source_fov_1_n = np.zeros(source_N).astype(np.uint8)
    source_fov_2_n = np.zeros(source_N).astype(np.uint8)
    source_fov_1_times = np.zeros((source_N,N_maxobs)).astype(np.uint32)
    source_fov_2_times = np.zeros((source_N,N_maxobs)).astype(np.uint32)
    t_diff = 1.0/24.0 # 1 hours
    time_idx_diff = 360
    r_search = np.tan(np.deg2rad((0.345+0.05+0.05)*np.sqrt(2)))
    b_fov = 0.345
    zeta_origin_1 = +220.997922127812/3600
    zeta_origin_2 = -220.997922127812/3600
    
    #l_fov_1 = 0.0
    #l_fov_2 = 106.5
    rad2deg = np.rad2deg(1.0)
    
    source_fov_1_pairs = [np.array(_pairs) for _pairs in source_tree.query_ball_tree(fov_1_tree, r = r_search)]
    source_fov_2_pairs = [np.array(_pairs) for _pairs in source_tree.query_ball_tree(fov_2_tree, r = r_search)]
    
    @njit
    def test_validity(xyz,source_fov_pairs,source_fov_times,zeta_origin):
        
        b_upp = zeta_origin+b_fov
        b_low = zeta_origin-b_fov

        time_previous_idx = -99999
        _n = 0
            
        for time_idx in source_fov_pairs:
                
            _xyz = xyz + gaia_velocity[time_idx]
            _norm = sqrt(_xyz[0]**2+_xyz[1]**2+_xyz[2]**2)
        
            _rot = _matrix[time_idx]
            #_U = (_rot[0,0] * _xyz[0] + _rot[0,1] * _xyz[1] + _rot[0,2] * _xyz[2])/_norm
            #_V = (_rot[1,0] * _xyz[0] + _rot[1,1] * _xyz[1] + _rot[1,2] * _xyz[2])/_norm
            _W = (_rot[2,0] * _xyz[0] + _rot[2,1] * _xyz[1] + _rot[2,2] * _xyz[2])/_norm
                
            #_l = rad2deg*atan2(_V,_U)
            _b = rad2deg*asin(_W)
                
            if _b < b_upp and _b > b_low and time_idx - time_previous_idx > time_idx_diff:
                source_fov_times[_n] = time_idx
                time_previous_idx = time_idx
                _n += 1
        
        return _n
    
    for source_idx in range(source_N):
        source_fov_1_n[source_idx] = test_validity(source_xyz[source_idx],source_fov_1_pairs[source_idx],source_fov_1_times[source_idx],zeta_origin_1)
        source_fov_2_n[source_idx] = test_validity(source_xyz[source_idx],source_fov_2_pairs[source_idx],source_fov_2_times[source_idx],zeta_origin_2)

        
    with h5py.File(f'./store/simulated_times_{block_idx}.h5', 'w') as f:
        
        # Create datasets
        f.create_dataset('source_id', data=source_id, compression = "lzf", chunks = (N_chunk,), shape = (N_chunk,), dtype = np.uint64, fletcher32 = False, shuffle = True, scaleoffset=0)
        f.create_dataset('fov_1_n', data=source_fov_1_n, compression = "lzf", chunks = (N_chunk,), shape = (N_chunk,), dtype = np.uint8, fletcher32 = False, shuffle = True, scaleoffset=0)
        f.create_dataset('fov_2_n', data=source_fov_2_n, compression = "lzf", chunks = (N_chunk,), shape = (N_chunk,), dtype = np.uint8, fletcher32 = False, shuffle = True, scaleoffset=0)
        f.create_dataset('fov_1_times', data=source_fov_1_times, compression = "lzf", chunks = (N_chunk, N_maxobs, ), shape = (N_chunk, N_maxobs, ), dtype = np.uint32, fletcher32 = False, shuffle = True, scaleoffset=0)
        f.create_dataset('fov_2_times', data=source_fov_2_times, compression = "lzf", chunks = (N_chunk, N_maxobs, ), shape = (N_chunk, N_maxobs, ), dtype = np.uint32, fletcher32 = False, shuffle = True, scaleoffset=0)    
    return block_idx

if True:
    
    with h5py.File('simulated_times.h5', 'w') as f:
        
        # Create datasets
        f.create_dataset('source_id', compression = "lzf", chunks = (N_chunk,), shape = (N_sources,), dtype = np.uint64, fletcher32 = False, shuffle = True, scaleoffset=0)
        f.create_dataset('fov_1_n', compression = "lzf", chunks = (N_chunk,), shape = (N_sources,), dtype = np.uint8, fletcher32 = False, shuffle = True, scaleoffset=0)
        f.create_dataset('fov_2_n', compression = "lzf", chunks = (N_chunk,), shape = (N_sources,), dtype = np.uint8, fletcher32 = False, shuffle = True, scaleoffset=0)
        f.create_dataset('fov_1_times', compression = "lzf", chunks = (N_chunk, N_maxobs, ), shape = (N_sources, N_maxobs, ), dtype = np.uint32, fletcher32 = False, shuffle = True, scaleoffset=0)
        f.create_dataset('fov_2_times', compression = "lzf", chunks = (N_chunk, N_maxobs, ), shape = (N_sources, N_maxobs, ), dtype = np.uint32, fletcher32 = False, shuffle = True, scaleoffset=0)


with h5py.File('simulated_times.h5', 'a') as f:
    #flag = f['source_id'][::N_chunk]
    
    
    pool = multiprocessing.Pool(N_core)
    for result in tqdm.tqdm(pool.imap(mp_worker, range(N_block)),total=N_block):
            
        block_idx = result
        
        with h5py.File(f'./store/simulated_times_{block_idx}.h5', 'r') as g:
        
            f['source_id'][block_idx*N_chunk:(block_idx+1)*N_chunk] = g['source_id'][:]
            f['fov_1_n'][block_idx*N_chunk:(block_idx+1)*N_chunk] = g['fov_1_n'][:]
            f['fov_2_n'][block_idx*N_chunk:(block_idx+1)*N_chunk] = g['fov_2_n'][:]
            f['fov_1_times'][block_idx*N_chunk:(block_idx+1)*N_chunk] = g['fov_1_times'][:]
            f['fov_2_times'][block_idx*N_chunk:(block_idx+1)*N_chunk] = g['fov_2_times'][:]
        
        os.remove(f'./store/simulated_times_{block_idx}.h5')
