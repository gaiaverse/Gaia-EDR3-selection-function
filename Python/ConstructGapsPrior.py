import pandas as pd
import numpy as np

# Constants
astrometric_start = 1192.13
astrometric_end = 5230.09

# Load in gaps
known_gaps = pd.read_csv('./edr3_gaps.csv')

# Load in scanning law
jd_time = np.sort(pd.read_csv('./CommandedScanLaw_001.csv',usecols=['jd_time'])['jd_time'].values)
print('There are',jd_time.size,'time points.')

# Compute OBMT
obmt = 1717.6256+(jd_time + 2455197.5 - 2457023.5 - 0.25)*4

# Find gaps
gaps = np.zeros(jd_time.size, dtype = np.int8)
gaps[obmt < astrometric_start] = 1
gaps[obmt > astrometric_end] = 1
for start,end in zip(known_gaps['tbeg'].values,known_gaps['tend'].values):
    in_gap = np.where((obmt > start) & (obmt < end))
    gaps[in_gap] = 1

# Find the consecutive ones
# Taken from https://stackoverflow.com/questions/24885092/finding-the-consecutive-zeros-in-a-numpy-array
def one_runs(a):
    # Pad each end with an extra 0.
    iszero = np.concatenate(([0], a, [0]))
    absdiff = np.abs(np.diff(iszero))
    # Runs start and end where absdiff is 1.
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    return ranges

# Save file
np.savetxt('./ModelInputs/gaps_prior.dat',one_runs(gaps),fmt='%d')