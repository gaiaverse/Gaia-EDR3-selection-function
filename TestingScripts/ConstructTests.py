
import numpy as np
import tqdm
import pandas as pd
import h5py
import healpy as hp
from scipy import special
import os
import sys

mode = -1
if len(sys.argv) > 1:
	mode = int(sys.argv[1])

largeNumber = 1000.0

class ValidationTestset:

    
    def __init__(self, testsetName, sourceNumber=1, magnitudeBinNumber=1, healpixOrderNumber=0, meanTime = largeNumber, meanMagnitudeSpace = largeNumber, minimumMeasurementNumber = 5):
        self.timestepNumber = 8967691
        
        self.testsetName = testsetName
        self.sourceNumber = int(sourceNumber)
        self.magnitudeBinNumber = int(magnitudeBinNumber)
        self.healpixOrderNumber = int(healpixOrderNumber)
        self.minimumMeasurementNumber = minimumMeasurementNumber
        
        self.healpixNsideNumber = hp.order2nside(self.healpixOrderNumber)
        self.healpixPixelNumber = hp.nside2npix(self.healpixNsideNumber)
        
        self.xTime = meanTime*np.ones(self.timestepNumber)
        self.xMagnitudeSpace = meanMagnitudeSpace*np.ones((self.magnitudeBinNumber,self.healpixPixelNumber))
        
        # Check it exists, if not then create
        self.directoryRoot = '/mnt/extraspace/GaiaSelectionFunction/'
        self.directoryModelInputs = self.directoryRoot + 'ModelInputs/'
        self.directoryTestset = self.directoryRoot + f'Data/TestSets3/{self.testsetName}/'
        if not os.path.exists(self.directoryTestset):
            os.makedirs(self.directoryTestset)

            
    def generateTestset(self):
        
        print(f'Generating {self.testsetName} testset.')
        
        with h5py.File(self.directoryModelInputs+'simulated_times.h5', 'r') as f:
            raw_times = {k:f[k][:self.sourceNumber] for k in f.keys()}
            self.data = {i:{'observations':int(raw_times['fov_1_n'][i])+int(raw_times['fov_2_n'][i]),
                            'times':np.sort(np.concatenate([raw_times['fov_1_times'][i,:raw_times['fov_1_n'][i]],raw_times['fov_2_times'][i,:raw_times['fov_2_n'][i]]]))} for i in tqdm.tqdm(range(self.sourceNumber))}
            del raw_times
        
        # Random assign magnitudes
       
        for i in tqdm.tqdm(range(self.sourceNumber)):
            raw = np.random.rand(1,1)
            transformed = 1.0/(1 + np.exp(-10*(raw - 0.5)))
            extractedG = np.floor(self.magnitudeBinNumber * transformed)
            # ~ self.data[i]['magnitude'] = np.random.randint(0,self.magnitudeBinNumber)
            self.data[i]['magnitude'] = extractedG
        # Calculate probabilities
        healpix_fov = pd.read_csv(self.directoryModelInputs+f'scanninglaw_to_healpix_{self.healpixOrderNumber}.csv')
        healpix_fov_1 = healpix_fov['fov_1_hpx']
        healpix_fov_2 = healpix_fov['fov_2_hpx']
        for i in tqdm.tqdm(range(self.sourceNumber)):
            self.data[i]['probabilities'] = special.expit(self.xTime[self.data[i]['times']])* \
                                        special.expit(self.xMagnitudeSpace[self.data[i]['magnitude'],healpix_fov_1[self.data[i]['times']]] \
                                                     +self.xMagnitudeSpace[self.data[i]['magnitude'],healpix_fov_2[self.data[i]['times']]])
        
        # Subsample
        for i in tqdm.tqdm(range(self.sourceNumber)):
            self.data[i]['measurements'] = np.sum(np.random.binomial(n=1,p=self.data[i]['probabilities']))
        
        # Add noise
        #if 'spikes' in test:
        #    print('Added spikes to '+test)
        #    # Add extra noise
        #    for i in tqdm.tqdm(range(N_sources)):
        #        if np.random.uniform(0,1) < variance_fraction:
        #            data[i]['k'] += np.round(np.random.normal(0,variance_baseline_1+variance_scaling_1*data[i]['n'])).astype(int)
        #        else:
        #            data[i]['k'] += np.round(np.random.normal(0,variance_baseline_2+variance_scaling_2*data[i]['n'])).astype(int)
                           
        # Cull
        skippable = []
        for i in tqdm.tqdm(range(self.sourceNumber)):
            if self.data[i]['measurements'] < self.minimumMeasurementNumber:
                skippable.append(i)
               
       
        # Write to file
        writeLine = {m:[] for m in range(self.magnitudeBinNumber)}
        for i in tqdm.tqdm(range(self.sourceNumber)):
            if i not in skippable:
                 writeLine[self.data[i]['magnitude']].append(f"{self.data[i]['measurements']},{self.data[i]['observations']},"+','.join(map(str, self.data[i]['times'])))

        for m in range(self.magnitudeBinNumber):
            with open(self.directoryTestset+f'{m}.csv', 'w') as f:
                 f.writelines("{}\n".format(line) for line in writeLine[m])
                
        # Write true parameters
        X = np.concatenate([self.xTime,self.xMagnitudeSpace.T.flatten()])
        np.savetxt(self.directoryTestset+'True_TransformedParameters.dat', np.c_[X])
        
        parameters = {}
        missingParameter = -1
        parameters['Nt'] = self.timestepNumber
        parameters['Nm'] = self.magnitudeBinNumber
        parameters['healpix_order'] = self.healpixOrderNumber
        parameters['needlet_order'] = missingParameter
        parameters['Nl'] = self.healpixPixelNumber
        parameters['Ns'] = missingParameter
        parameters['totalRawParams'] = missingParameter
        parameters['totalTransformedParams'] = parameters['Nt'] + parameters['Nm']*parameters['Nl']
        parameters['mu_t'] = missingParameter
        parameters['sigma_t'] = missingParameter
        parameters['l_m'] = missingParameter
        parameters['l_t'] = missingParameter
        line_header,line_values = [], []
        for k,v in parameters.items():
            line_header.append(k)
            line_values.append(v)
        line_header = ','.join(map(str, line_header))+'\n'
        line_values = ','.join(map(str, line_values))
        
        with open(self.directoryTestset+f'Optimiser_Properties.dat', 'w') as f:
            f.write(line_header)
            f.write(line_values)
            
    def applyEdr3Gaps(self):
        timesteps = np.linspace(1666.4384902198801, 2704.3655735533684, self.timestepNumber)
        obmt = 1717.6256+(timesteps + 2455197.5 - 2457023.5 - 0.25)*4
        gaps = pd.read_csv(self.directoryModelInputs+'edr3_gaps.csv')
        for start,end in zip(gaps['tbeg'].values,gaps['tend'].values):
            self.xTime[(obmt >= start) & (obmt <= end)] = -largeNumber
            
    def applyGalacticCrowding(self):
        gaia = -np.log(hp.read_map(self.directoryModelInputs+'densityMap-I_350_gaiaedr3-128.hpx'))
        m = 10/(gaia.max()-gaia.min())
        b = 5 - m*gaia.max()
        gaia = m * gaia + b
        gaia = hp.ud_grade(gaia,nside_out=self.healpixNsideNumber)
        self.xMagnitudeSpace += gaia[np.newaxis,:]
        
    def applyMagnitudeTrend(self,magnitudeArray):
        assert magnitudeArray.size == self.magnitudeBinNumber
        y = -0.075*(magnitudeArray-12)**2/0.5+10
        self.xMagnitudeSpace += y[:,np.newaxis]
        

NMax = 1e4

if mode == 0:
	# Flat recovery
	flatTest = ValidationTestset(testsetName='flat', sourceNumber=NMax, meanTime = 0.0)
	flatTest.generateTestset()

if mode == 1:
	# Gaps recovery
	gapsTest = ValidationTestset(testsetName='gaps', sourceNumber=NMax)
	gapsTest.applyEdr3Gaps()
	gapsTest.generateTestset()


g_bins = np.arange(1.7,23.05,0.2)
g_midbins = 0.5*(g_bins[1:]+g_bins[:-1])

if mode == 2:
	# Galaxy recovery
	galaxyTest = ValidationTestset(testsetName='galaxy', sourceNumber=NMax, magnitudeBinNumber=len(g_midbins), meanMagnitudeSpace=0.0, healpixOrderNumber=7)
	galaxyTest.applyGalacticCrowding()
	galaxyTest.applyMagnitudeTrend(g_midbins)
	galaxyTest.generateTestset()

if mode == 3:
	# Full recovery
	fullTest = ValidationTestset(testsetName='full', sourceNumber=NMax, magnitudeBinNumber=len(g_midbins), meanMagnitudeSpace=0.0, healpixOrderNumber=7)
	fullTest.applyEdr3Gaps()
	fullTest.applyGalacticCrowding()
	fullTest.applyMagnitudeTrend(g_midbins)
	fullTest.generateTestset()
