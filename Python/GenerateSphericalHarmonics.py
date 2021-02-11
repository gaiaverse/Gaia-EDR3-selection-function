""" Compute the real-valued spherical harmonics for nested Healpix maps. """

import numpy as np
import healpy as hp
import h5py
import sys
import tqdm

nside = int(sys.argv[1])
lmax = int(sys.argv[2])
Npix = hp.nside2npix(nside)

# Form the l's and m's
l = np.concatenate([np.repeat(_l,2*_l+1) for _l in range(lmax+1)]).astype(np.int)
m = np.concatenate([np.arange(-_l,_l+1e-9) for _l in range(lmax+1)]).astype(np.int)
l_hp,m_hp = hp.sphtfunc.Alm.getlm(lmax=lmax)
Nmodes = int((lmax+1)**2)
Nmodes_healpy = int((lmax+1)*(lmax+2)/2)

Ylm = np.zeros((Nmodes,Npix))
alm_hp = np.zeros(Nmodes_healpy)
for i,(_l,_m) in enumerate(zip(tqdm.tqdm(l),m)):
    if _m > 0.0:
        hp_idx = np.where((l_hp==_l) & (m_hp == _m))[0]
        alm_hp[hp_idx] = 1.0
        Ylm[i] = np.power(-1,_m)*hp.pixelfunc.reorder(hp.sphtfunc.alm2map(alm_hp+0.j*alm_hp,nside=nside,verbose=False)/np.sqrt(2), r2n=True)
        alm_hp[hp_idx] = 0.0
    elif _m < 0.0:
        hp_idx = np.where((l_hp==_l) & (m_hp == -_m))[0]
        alm_hp[hp_idx] = 1.0
        Ylm[i] = np.power(-1,-_m+1)*hp.pixelfunc.reorder(hp.sphtfunc.alm2map(0.*alm_hp+1.j*alm_hp,nside=nside,verbose=False)/np.sqrt(2), r2n=True)
        alm_hp[hp_idx] = 0.0
    else:
        hp_idx = np.where((l_hp==_l) & (m_hp == 0))[0]
        alm_hp[hp_idx] = 1.0
        Ylm[i] = hp.pixelfunc.reorder(hp.sphtfunc.alm2map(alm_hp+0.j*alm_hp,nside=nside,verbose=False), r2n=True)
        alm_hp[hp_idx] = 0.0
        
with h5py.File(f'sphericalharmonics_nside{nside}_lmax{lmax}.h5', 'w') as f:
        
    # Create datasets
    f.create_dataset('Ylm', data = Ylm, compression = "lzf", chunks = True, shape = (Nmodes, Npix,), dtype = np.float64, fletcher32 = False, shuffle = True)
    f.create_dataset('l', data = l, compression = "lzf", chunks = True, shape = (Nmodes, ), dtype = np.uint32, fletcher32 = False, shuffle = True, scaleoffset=0)
    f.create_dataset('m', data = m, compression = "lzf", chunks = True, shape = (Nmodes,), dtype = np.uint32, fletcher32 = False, shuffle = True, scaleoffset=0)
        
