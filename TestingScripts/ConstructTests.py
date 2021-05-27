import numpy as np
import tqdm
import pandas as pd
import h5py
import healpy as hp
from numba import njit
from scipy import special
import os

NT = 8967691
reset = False

for test in ['zeros','ones','binary','gaps','gaps_parallel','magnitudes','magnitudes_sky']:
    print(test)
    if test == 'zeros':
        N_sources = 10000
        Nm = 1
        healpix_order = 0
        Nl = hp.nside2npix(hp.order2nside(healpix_order))

        xt = np.zeros(NT)
        xml = np.zeros((Nm,Nl))
    elif test == 'ones':
        N_sources = 10000
        Nm = 1
        healpix_order = 0
        Nl = hp.nside2npix(hp.order2nside(healpix_order))

        xt = np.ones(NT)
        xml = np.ones((Nm,Nl))
    elif test == 'binary':
        N_sources = 10000
        Nm = 1
        healpix_order = 0
        Nl = hp.nside2npix(hp.order2nside(healpix_order))

        xt = 5.0*np.ones(NT)
        xt[:int(NT/2)] = -5.0
        xml = 5.0*np.ones((Nm,Nl))
    elif test == 'gaps':
        N_sources = 100000
        Nm = 1
        healpix_order = 0
        Nl = hp.nside2npix(hp.order2nside(healpix_order))
        xml = 5.0*np.ones((Nm,Nl))
        
        xt = 5.0*np.ones(NT)
        t = np.linspace(1666.4384902198801, 2704.3655735533684, NT)
        obmt = 1717.6256+(t + 2455197.5 - 2457023.5 - 0.25)*4
        gaps = pd.read_csv('./edr3_gaps.csv')
        for start,end in zip(gaps['tbeg'].values,gaps['tend'].values):
            xt[(obmt >= start) & (obmt <= end)] = -5.0
    elif test == 'gaps_parallel':
        N_sources = 100000
        Nm = 10
        healpix_order = 0
        Nl = hp.nside2npix(hp.order2nside(healpix_order))
        xml = 5.0*np.ones((Nm,Nl))
        
        xt = 5.0*np.ones(NT)
        t = np.linspace(1666.4384902198801, 2704.3655735533684, NT)
        obmt = 1717.6256+(t + 2455197.5 - 2457023.5 - 0.25)*4
        gaps = pd.read_csv('./edr3_gaps.csv')
        for start,end in zip(gaps['tbeg'].values,gaps['tend'].values):
            xt[(obmt >= start) & (obmt <= end)] = -5.0
    elif test == 'magnitudes':
        N_sources = 100000
        Nm = 10
        healpix_order = 0
        Nl = hp.nside2npix(hp.order2nside(healpix_order))
        xt = 5.0*np.ones(NT)
        xml = np.repeat(np.linspace(-1.0,1.0,Nm)[:,np.newaxis],Nl,axis=1)
    elif test == 'magnitudes_sky':
        N_sources = 100000
        Nm = 10
        healpix_order = 0
        Nl = hp.nside2npix(hp.order2nside(healpix_order))
        xt = 5.0*np.ones(NT)
        xml = np.repeat(np.linspace(-1.0,1.0,Nm)[:,np.newaxis],Nl,axis=1)
        
        healpix_nside = hp.order2nside(healpix_order)
        healpix_npix = hp.nside2npix(healpix_nside)
        l, b = hp.pix2ang(nside=healpix_nside,ipix=np.arange(healpix_npix),lonlat=True,nest=False)
        for m in range(Nm):
            xml[m] *= np.cos(np.deg2rad(b))
    else:
        print(f"{test} not defined.")
        continue

    # Check it exists, if not then create
    directory = f'./TestSets/{test}/'
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        if reset == False:
            print(f"Not resetting {test}.")
            continue

    # Load in simulated times
    with h5py.File('simulated_times.h5', 'r') as f:
        raw_times = {k:f[k][:N_sources] for k in f.keys()}
        data = {i:{'n':int(raw_times['fov_1_n'][i])+int(raw_times['fov_2_n'][i]),'times':np.sort(np.concatenate([raw_times['fov_1_times'][i,:raw_times['fov_1_n'][i]],raw_times['fov_2_times'][i,:raw_times['fov_2_n'][i]]]))} for i in tqdm.tqdm(range(N_sources))}
        del raw_times

    # Random assign magnitudes
    for i in tqdm.tqdm(range(N_sources)):
        data[i]['m'] = np.random.randint(0,Nm)
    
    # Calculate p
    healpix_fov = pd.read_csv(f'./ModelInputs/scanninglaw_to_healpix_{healpix_order}.csv')
    healpix_fov_1 = healpix_fov['fov_1_hpx']
    healpix_fov_2 = healpix_fov['fov_2_hpx']
    for i in tqdm.tqdm(range(N_sources)):
        t = data[i]['times']
        xl = xml[data[i]['m']]
        data[i]['p'] = special.expit(xt[t])*special.expit(xl[healpix_fov_1[t]]+xl[healpix_fov_2[t]])

    # Subsample
    for i in tqdm.tqdm(range(N_sources)):
        data[i]['k'] = np.sum(np.random.binomial(n=1,p=data[i]['p']))

    # Write to file
    source_lines = {m:[] for m in range(Nm)}
    for i in tqdm.tqdm(range(N_sources)):
        if data[i]['k'] >= 5:
            source_lines[data[i]['m']].append(f"{data[i]['k']},{data[i]['n']},"+','.join(map(str, data[i]['times'])))

    for m in range(Nm):

        with open(directory+f'{m}.csv', 'w') as f:

            # Write lines
            f.writelines("{}\n".format(line) for line in source_lines[m])