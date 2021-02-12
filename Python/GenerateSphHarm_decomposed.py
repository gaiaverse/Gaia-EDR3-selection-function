""" Compute the real-valued spherical harmonics for nested Healpix maps. """

import sys
import numpy as np, healpy as hp, scipy.stats
import h5py, tqdm


def gen_lambda(lmax, l, m, nside, hpx_ring_0):

    alm_hp = np.zeros(int((lmax+1)*(lmax+2)/2))
    _lambda = np.zeros((int((lmax+1)**2),nside))

    for i,(_l,_m) in enumerate(zip(tqdm.tqdm(l),m)):
        if _m > 0.0:
            hp_idx = np.where((l_hp==_l) & (m_hp == _m))[0]
            alm_hp[hp_idx] = 1.0
            _lambda[i] = hp.sphtfunc.alm2map(alm_hp+0.j*alm_hp,nside=nside,verbose=False)[hpx_ring_0]
            alm_hp[hp_idx] = 0.0
        elif _m < 0.0:
            hp_idx = np.where((l_hp==_l) & (m_hp == -_m))[0]
            alm_hp[hp_idx] = 1.0
            _lambda[i] = hp.sphtfunc.alm2map(0.*alm_hp+1.j*alm_hp,nside=nside,verbose=False)[hpx_ring_0]
            alm_hp[hp_idx] = 0.0
        else:
            hp_idx = np.where((l_hp==_l) & (m_hp == 0))[0]
            alm_hp[hp_idx] = 1.0
            _lambda[i] = hp.sphtfunc.alm2map(alm_hp+0.j*alm_hp,nside=nside,verbose=False)[hpx_ring_0]
            alm_hp[hp_idx] = 0.0

    return _lambda

def gen_exponential(lmax, npix, ring_id, ring_x, ring_phi0, N_ring):

    _exponential = np.zeros((2*lmax+1,npix)) + 0.j

    for m in range(-lmax, lmax+1):
        _exponential[m+lmax] = np.exp(1.j*m* (ring_phi0[ring_id] + 2*np.pi*ring_x/N_ring[ring_id]))

    return _exponential


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
hpx_ring_0 = np.argwhere(phi==0)[:,0]

# Number of pixels in the ring of the given pixel
Npix_ring = np.hstack((np.arange(4,nside*4+1,4),
                       np.zeros(nside*2-1).astype(int)+nside*4,
                       np.arange(4,nside*4+1,4)[::-1]))
# Ring ID
ring_id = np.array([], dtype=int)
for i in range(4*nside-1): ring_id=np.hstack((ring_id, np.zeros(Npix_ring[i], dtype=int)+i))
# Minimum phi in the ring of the given pixel
ring_phi0 = scipy.stats.binned_statistic(ring_id, phi, statistic='min', bins=np.arange(nside*4)-0.5).statistic
# Index x of pixel in ring
ring_x = np.array([], dtype=int)
for i in range(4*nside-1): ring_x=np.hstack((ring_x, np.arange(Npix_ring[i])))

_lambda = gen_lambda(lmax, l, m, nside, hpx_ring_0)
_exponential = gen_exponential(lmax, Npix, ring_id, ring_x, ring_phi0, Npix_ring)


with h5py.File('/data/asfe2/Projects/gaia_edr3/sphericalharmonics_decomposed_nside{0}_lmax{1}.h5'.format(nside,lmax), 'w') as f:
    # Create datasets
    f.create_dataset('lambda', data = _lambda, compression = "lzf", chunks = True, shape = (lmax, Npix,), dtype = np.float64, fletcher32 = False, shuffle = True)
    f.create_dataset('exponential', data = _exponential, compression = "lzf", chunks = True, shape = (Nmodes, nside, ), dtype = np.complex128, fletcher32 = False, shuffle = True, scaleoffset=0)
    f.create_dataset('l', data = l, compression = "lzf", chunks = True, shape = (Nmodes, ), dtype = np.uint32, fletcher32 = False, shuffle = True, scaleoffset=0)
    f.create_dataset('m', data = m, compression = "lzf", chunks = True, shape = (Nmodes,), dtype = np.uint32, fletcher32 = False, shuffle = True, scaleoffset=0)
