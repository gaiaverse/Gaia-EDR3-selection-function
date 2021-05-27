""" For each Gaia chunk, compute n observations. """

import numpy as np
import sys
import tqdm
import pandas as pd
import os
import healpy as hp
            
##### Load in sources and assign to bins
N_chunk = 11941
N_block = 751
directory = './ModelInputs/'

# Check it exists, if not then create
if not os.path.exists(directory):
    os.makedirs(directory)

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
del _box
print('Loaded scanning law.')

fov_1_hpx = np.zeros(N_chunk*N_block,dtype=int)
fov_2_hpx = np.zeros(N_chunk*N_block,dtype=int)
fov_times = np.arange(N_chunk*N_block,dtype=int)

for order in range(8):
    nside = hp.order2nside(order)
    print(f'Working on order={order}, nside={nside}.')
    
    for i_block in tqdm.tqdm(range(N_block)):
        start = i_block * N_chunk
        end = start + N_chunk
        fov_1_hpx[start:end] = hp.vec2pix(x=fov_1_xyz[start:end,0],y=fov_1_xyz[start:end,1],z=fov_1_xyz[start:end,2],nside=nside,nest=False)
        fov_2_hpx[start:end] = hp.vec2pix(x=fov_2_xyz[start:end,0],y=fov_2_xyz[start:end,1],z=fov_2_xyz[start:end,2],nside=nside,nest=False)
        
    df = pd.DataFrame({"fov_times" : fov_times, "fov_1_hpx" : fov_1_hpx, "fov_2_hpx": fov_2_hpx})
    df.to_csv(directory+f'scanninglaw_to_healpix_{order}.csv', index=False)