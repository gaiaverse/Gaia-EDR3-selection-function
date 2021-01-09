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

##### Script for saving data to pickle
def save_as_pickled_object(obj, filepath):
    """
    This is a defensive way to write pickle.write, allowing for very large files on all platforms
    """
    max_bytes = 2**31 - 1
    bytes_out = pickle.dumps(obj)
    n_bytes = sys.getsizeof(bytes_out)
    with open(filepath, 'wb') as f_out:
        for idx in range(0, n_bytes, max_bytes):
            f_out.write(bytes_out[idx:idx+max_bytes])
            
##### Load in sources and assign to bins
N_core = 1

bins = np.array([0.0, 5.0, 16.0, 17.0, 18.0, 18.2, 18.4, 18.6, 18.8, 19.0, 19.1, 19.2, 19.3, 19.4, 19.5, 19.6, 19.7, 19.8, 19.9, 20.0, 20.1, 20.2, 20.3, 20.4, 20.5, 20.6, 20.7, 20.8, 20.9, 21.0, 21.1, 21.2, 21.3, 21.4, 21.5, 25.0])
N_bins = bins.size-1

source_keys = ['ra','dec','phot_g_mean_mag','matched_transits']
source_box = {}
with h5py.File('./gaiaedr3.h5', 'r') as f:
    for key in source_keys:
        source_box[key] = f[key][:1000000]
source_box['matched_transits'] = source_box['matched_transits'].astype(np.int)
print('Loaded sources')
        
input_box = {}
#source_index = np.searchsorted(bins,source_box['G'])-1
for i_bin in range(N_bins):
    in_bin = np.where( (source_box['phot_g_mean_mag'] >= bins[i_bin]) & (source_box['phot_g_mean_mag'] < bins[i_bin+1]) )
    input_box[i_bin] = [i_bin]+[source_box[key][in_bin] for key in ['ra','dec','matched_transits']]
del source_box
print('Created inputs')

##### Loop through sources

def mp_worker(args):
    
    i_bin, ra_source, dec_source, k_source = args
    
    N_source = ra_source.size
    
    ##### Compute tree
    xyz_source = np.stack([np.cos(np.deg2rad(ra_source))*np.cos(np.deg2rad(dec_source)), np.sin(np.deg2rad(ra_source))*np.cos(np.deg2rad(dec_source)), np.sin(np.deg2rad(dec_source))]).T
    tree_source = spatial.cKDTree(xyz_source)
    
    from scipy import interpolate

    # Define units
    speed_of_light_AU_per_day = 299792458.0*(86400.0/149597870.700/1e3)

    # Prepare ephem data
    gaia_ephem_data = pd.read_csv('./horizons_results_gaia.txt',skiprows=64)
    gaia_ephem_box = {k:gaia_ephem_data[k].values for k in ['JDTDB','X','Y','Z','VX','VY','VZ']}
    gaia_ephem_velocity = interpolate.interp1d( gaia_ephem_box['JDTDB']-2455197.5, np.stack([gaia_ephem_box['VX'],gaia_ephem_box['VY'],gaia_ephem_box['VZ']])/speed_of_light_AU_per_day, kind='cubic')
    
    ##### Load in scanning law
    _data = pd.read_csv('./CommandedScanLaw_001.csv')
    _columns = ['jd_time', 'bjd_fov1', 'bjd_fov2','ra_fov1', 'dec_fov1', 'scan_angle_fov1', 'ra_fov2', 'dec_fov2', 'scan_angle_fov2']
    _keys = ['tcb_at_gaia','tcb_at_bary1','tcb_at_bary2','ra_fov_1','dec_fov_1','angle_fov_1','ra_fov_2','dec_fov_2','angle_fov_2']
    #print('There are',_box['tcb_at_gaia'].size,'time points.')
    _box = {}
    _order = np.argsort(_data['jd_time'].values)
    for j,k in zip(_columns,_keys):
        _box[k] = _data[j].values[_order]
    
    xyz_fov_1 = np.stack([np.cos(np.deg2rad(_box['ra_fov_1']))*np.cos(np.deg2rad(_box['dec_fov_1'])),np.sin(np.deg2rad(_box['ra_fov_1']))*np.cos(np.deg2rad(_box['dec_fov_1'])),np.sin(np.deg2rad(_box['dec_fov_1']))]).T
    xyz_fov_2 = np.stack([np.cos(np.deg2rad(_box['ra_fov_2']))*np.cos(np.deg2rad(_box['dec_fov_2'])),np.sin(np.deg2rad(_box['ra_fov_2']))*np.cos(np.deg2rad(_box['dec_fov_2'])),np.sin(np.deg2rad(_box['dec_fov_2']))]).T
    #print('Loaded scanning law.')

    ##### Compute rotation matrices
    _xaxis = xyz_fov_1
    _zaxis = np.cross(xyz_fov_1,xyz_fov_2)
    _yaxis = -np.cross(_xaxis,_zaxis)
    _yaxis /= np.linalg.norm(_yaxis,axis=1)[:,np.newaxis]
    _zaxis /= np.linalg.norm(_zaxis,axis=1)[:,np.newaxis]

    _uaxis = np.array([1,0,0])
    _vaxis = np.array([0,1,0])
    _waxis = np.array([0,0,1])

    _matrix = np.moveaxis(np.stack([np.dot(_xaxis,_uaxis),np.dot(_xaxis,_vaxis),np.dot(_xaxis,_waxis),
          np.dot(_yaxis,_uaxis),np.dot(_yaxis,_vaxis),np.dot(_yaxis,_waxis),
          np.dot(_zaxis,_uaxis),np.dot(_zaxis,_vaxis),np.dot(_zaxis,_waxis)]).reshape((3,3,_xaxis.shape[0])),2,0)

    #print('Computed rotation matrix.')


    
    ##### Find intersections
    n_source = np.zeros(N_source)
    t_previous = -99999*np.ones(N_source)
    t_diff = 1.0/24.0 # 1 hours
    r_search = np.tan(np.deg2rad((0.345+0.05)*np.sqrt(2)))
    
    from numba.typed import List
    t_lists_1, t_lists_2 = List(), List()
    for _idx in range(N_source):
        l_1 = List()
        l_2 = List()
        l_1.append(-1.0)
        l_2.append(-1.0)
        t_lists_1.append(l_1)
        t_lists_2.append(l_2)

    b_fov = 0.345
    zeta_origin_1 = +220.997922127812/3600
    zeta_origin_2 = -220.997922127812/3600
    b_upp_1 = b_fov+zeta_origin_1
    b_low_1 = -b_fov+zeta_origin_1
    b_upp_2 = b_fov+zeta_origin_2
    b_low_2 = -b_fov+zeta_origin_2
    l_fov_1 = 0.0
    l_fov_2 = 106.5
    
    from numba import njit

    @njit
    def test_validity(_tidx,_in_this_fov,_n_source,_t_previous,_t_lists_1,_t_lists_2,_t_now,_gaia_velocity):
        #_xyz = xyz_source[_in_this_fov].T
        #_uvw = np.dot(_matrix[_tidx],_xyz)
        #_l = np.rad2deg(np.arctan2(_uvw[1],_uvw[0]))
        #_b = np.rad2deg(np.arctan2(_uvw[2],np.sqrt(_uvw[0]**2.0+_uvw[1]**2.0)))

        _xyz = xyz_source[_in_this_fov] + _gaia_velocity
        #_xyz = _xyz/np.linalg.norm(_xyz,axis=1)[:,np.newaxis]
        _norm = np.sqrt(_xyz[:,0]**2+_xyz[:,1]**2+_xyz[:,2]**2)
        
        #_xyz = xyz_source[_in_this_fov]
        _rot = _matrix[_tidx]
        _U = (_rot[0,0] * _xyz[:,0] + _rot[0,1] * _xyz[:,1] + _rot[0,2] * _xyz[:,2])/_norm
        _V = (_rot[1,0] * _xyz[:,0] + _rot[1,1] * _xyz[:,1] + _rot[1,2] * _xyz[:,2])/_norm
        _W = (_rot[2,0] * _xyz[:,0] + _rot[2,1] * _xyz[:,1] + _rot[2,2] * _xyz[:,2])/_norm
        _l = np.rad2deg(np.arctan2(_V,_U))
        _b = np.rad2deg(np.arcsin(_W))

        _valid_fov_1 =  _in_this_fov[np.where(((np.abs(_l-l_fov_1)<1.0)&(_b<b_upp_1)&(_b>b_low_1))&(_t_now-_t_previous[_in_this_fov]>t_diff))[0]]
        _n_source[_valid_fov_1] += 1
        _t_previous[_valid_fov_1] = _t_now
        for _sidx in _valid_fov_1:
            _t_lists_1[_sidx].append(_tidx)

        _valid_fov_2 =  _in_this_fov[np.where(((np.abs(_l-l_fov_2)<1.0)&(_b<b_upp_2)&(_b>b_low_2))&(_t_now-_t_previous[_in_this_fov]>t_diff))[0]]
        _n_source[_valid_fov_2] += 1
        _t_previous[_valid_fov_2] = _t_now
        for _sidx in _valid_fov_2:
            _t_lists_2[_sidx].append(_tidx)

        return _n_source, _t_previous

    for _tidx in tqdm.tqdm(range(0,xyz_fov_1.shape[0])):
        _in_fov = tree_source.query_ball_point([xyz_fov_1[_tidx].copy(order='C'),xyz_fov_2[_tidx].copy(order='C')],r_search)
        _in_fov = np.array(_in_fov[0]+_in_fov[1])
        if _in_fov.size < 1:
            continue
        _gaia_velocity = gaia_ephem_velocity(_box['tcb_at_gaia'][_tidx])
        n_source, t_previous = test_validity(_tidx,_in_fov,n_source,t_previous,t_lists_1,t_lists_2,_box['tcb_at_gaia'][_tidx],_gaia_velocity)
    
    
    
    ##### Output
    filename = 'ObservationTimes_'+str(i_bin)+'.csv'
    #print(t_lists_1)
    with open(filename, 'a') as f:
        for i in range(N_source):
            _times = list(t_lists_1[i][1:])+list(t_lists_2[i][1:])
            _times.sort()
            if len(_times) == 0:
                continue
            line = str(k_source[i])+','+str(int(n_source[i]))+','+','.join(map(str, map(int,_times)))+'\n'
            f.write(line)
    return 1

def mp_handler():
    pool = multiprocessing.Pool(N_core)
    _input = [input_box[i_bin] for i_bin in [2]]#range(N_bins)]
    #result = list(tqdm.tqdm(pool.imap_unordered(mp_worker, _input), total=N_bins))
    result = list(pool.imap(mp_worker, _input))
    return result

_output = mp_handler()
