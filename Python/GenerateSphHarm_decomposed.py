""" Compute the real-valued spherical harmonics for nested Healpix maps. """

import sys
import numpy as np, healpy as hp, scipy.stats
import h5py, tqdm

nside = int(sys.argv[1])
lmax = int(sys.argv[2])
Npix = hp.nside2npix(nside)

# Form the l's and m's
l = np.concatenate([np.repeat(_l,2*_l+1) for _l in range(lmax+1)]).astype(np.int)
m = np.concatenate([np.arange(-_l,_l+1e-9) for _l in range(lmax+1)]).astype(np.int)
l_hp,m_hp = hp.sphtfunc.Alm.getlm(lmax=lmax)
Nmodes = int((lmax+1)**2)
Nmodes_healpy = int((lmax+1)*(lmax+2)/2)

# Ring idxs of pixels with phi=0
theta, phi = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
theta_ring = np.unique(theta)

# Generate Lambda
alm_hp = np.zeros(Nmodes_healpy)
_lambda = np.zeros((Nmodes, 4*nside-1))
for i,(_l,_m) in enumerate(zip(tqdm.tqdm(l),m)):
    hp_idx = np.where((l_hp==_l) & (m_hp == _m))[0]
    alm_hp[hp_idx] = 1.0
    _lambda[i] = np.real( scipy.special.sph_harm(np.abs(_m), _l, theta_ring*0., theta_ring) )
    alm_hp[hp_idx] = 0.0

# Generate Exponential
azimuth = np.ones((2*lmax+1,Npix))
for _m in range(-lmax, lmax+1):
    if _m<0:   azimuth[_m+lmax] = np.sin(-_m*phi)
    elif _m>0: azimuth[_m+lmax] = np.cos(_m*phi)
    else: pass

save_kwargs = {'compression':"lzf", 'chunks':True, 'fletcher32':False, 'shuffle':True}
with h5py.File('/data/asfe2/Projects/gaia_edr3/sphericalharmonics_decomposed_nside{0}_lmax{1}.h5'.format(nside,lmax), 'w') as f:
    # Create datasets
    f.create_dataset('lambda', data = _lambda, shape = (Nmodes, 4*nside-1,), dtype = np.float64, **save_kwargs)
    f.create_dataset('azimuth',data = azimuth, shape = (2*lmax+1, Npix, ),   dtype = np.float64, **save_kwargs)
    f.create_dataset('l',      data = l,       shape = (Nmodes,), dtype = np.uint32, scaleoffset=0, **save_kwargs)
    f.create_dataset('m',      data = m,       shape = (Nmodes,), dtype = np.uint32, scaleoffset=0, **save_kwargs)
