""" For each Gaia chunk, compute n observations. """

import numpy as np
import sys
import tqdm
import pandas as pd
import os
import healpy as hp
from numba import njit
from math import sin, cos
from scipy import sparse
import copy
            
##### Load in sources and assign to bins
B = 2.0
p = 1.0
needle_sparse_tol = 1e-4
directory = './ModelInputs/'

# Check it exists, if not then create
if not os.path.exists(directory):
    os.makedirs(directory)

##### Load in data that all will need

# Function to compute needlet across sky
@njit
def pixel_space (Y, cos_gamma, window, start, end, legendre):
    '''Return the value of a needlet at gamma radians from the needlet centre.'''

    legendre[0] = 1.0
    legendre[1] = cos_gamma
    for cur_l in range(2, end + 1):
        legendre[cur_l] = ((cos_gamma * (2 * cur_l - 1) * legendre[cur_l - 1] - (cur_l - 1) * legendre[cur_l - 2])) / cur_l

    Y[:] = np.dot(window,legendre[start:end+1])

# Needlet weighting function
class chisquare:
    
    def __init__(self, p = 1.0, B = 2.0, F = 1e-6):
        self.p = p
        self.B = B
        self.F = F

    def window_function(self, l, j):
        u = l*(l+1) / np.power(self.B,2.0*j)
        return np.power(u,self.p)*np.exp(-u)*(2.0*l+1.0)/(4.0*np.pi)
            
    def start(self, j):
        return 1
    
    def end(self, j):
        from scipy import special
        G = -self.p*special.lambertw(-np.power(self.F,1.0/self.p)/np.e,k=-1).real*np.power(self.B,2.0*j)
        return int(np.ceil(0.5*(-1.0+np.sqrt(1.0+4.0*G))))
    
    def normalise(self, j):
        N = 0
        for l in np.arange(self.start(j), self.end(j) + 1, dtype = 'float'):
            N += self.window_function(l,j)
        return N

weighting = chisquare(p = p, B = B)

for healpix_order in range(8):
    
    # Compute locations of pixels
    nside = hp.order2nside(healpix_order)
    print(f'Working on order={healpix_order}, nside={nside}.')
    npix = hp.nside2npix(nside)
    colat, lon = np.array(hp.pix2ang(nside=nside,ipix=np.arange(npix),lonlat=False))
    cos_colat, sin_colat = np.cos(colat), np.sin(colat)
    cos_lon, sin_lon = np.cos(lon), np.sin(lon)
    J = [-1] + [j for j in range(healpix_order+1)]

    # Initialise variables
    running_index = 0
    needlet_w, needlet_v, needlet_u = [], [], []
    Y = np.zeros(npix)
    legendre = np.zeros((1+weighting.end(max(J)),npix))

    for needlet_order in J:

        print(f'Working on order {needlet_order}.')

        if needlet_order == -1:
            needlet_w.append(np.ones(npix))
            needlet_v.append(np.arange(npix))
            needlet_u.append(0)
            running_index += npix
        else:

            needlet_nside = hp.order2nside(needlet_order)
            needlet_npix = hp.nside2npix(needlet_nside)

            start = weighting.start(needlet_order)
            end = weighting.end(needlet_order)
            modes = np.arange(start, end + 1, dtype = 'float')
            window = np.sqrt(4.0*np.pi/needlet_npix)*weighting.window_function(modes,needlet_order)/weighting.normalise(needlet_order)

            for needlet_ipix in tqdm.tqdm(range(needlet_npix),file=sys.stdout):

                colat_needle, lon_needle = hp.pix2ang(nside=needlet_nside,ipix=needlet_ipix,lonlat=False)

                cos_gamma = cos(colat_needle) * cos_colat + sin(colat_needle) * sin_colat * (cos(lon_needle) * cos_lon + sin(lon_needle) * sin_lon)

                pixel_space(Y, cos_gamma = cos_gamma, window = window, start = start, end = end, legendre = legendre)

                _significant = np.where(np.abs(Y) > Y.max()*needle_sparse_tol)[0]
                needlet_w.append(Y[_significant])
                needlet_v.append(_significant)
                needlet_u.append(running_index)
                running_index += _significant.size
        
        # Add the ending index to u
        needlet_u.append(running_index)
        
        # Concatenate the lists
        wavelet_w = np.concatenate(copy.deepcopy(needlet_w))
        wavelet_v = np.concatenate(copy.deepcopy(needlet_v))
        wavelet_u = np.array(copy.deepcopy(needlet_u))
        
        # Flip them round
        csr_matrix = sparse.csr_matrix((wavelet_w,wavelet_v,wavelet_u)).transpose().tocsr()
        wavelet_w, wavelet_v, wavelet_u = csr_matrix.data, csr_matrix.indices, csr_matrix.indptr
        
        # Convert u to long-form
        wavelet_u = np.concatenate([np.repeat(i,wavelet_u[i+1]-wavelet_u[i]) for i in range(wavelet_u.size-1)])
        
        # Save the needlets
        df = pd.DataFrame({"u" : wavelet_u, "v" :wavelet_v, "w": wavelet_w})
        df.to_csv(directory+f'needlets_{healpix_order}_{needlet_order}.csv', index=False)
        
        # Remove the ending index
        needlet_u.pop()        
    