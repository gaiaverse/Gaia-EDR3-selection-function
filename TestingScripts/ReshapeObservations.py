import tqdm
from contextlib import ExitStack
import multiprocessing
import numpy as np
import h5py

N_core = 1
N_maxobs = 256
N_mag = 213+1 # The +1 is to give a bin for the stars without G magnitudes to sit

def mp_worker(args):
    
    mag_idx = args

    with open(f'/mnt/extraspace/GaiaSelectionFunction/NewData/gaiaedr3_selection_function_files/{mag_idx}.csv', 'r') as f:
        
        # Empty lists for storage
        F = []
        K = []
        N = []
        T = []
        
        # Load in data
        for line in f:
            if line.strip():
                line_list = line.strip().split(",")
                F.append(line_list[0])
                K.append(line_list[1])
                N.append(line_list[2])
                T.append(line_list[3:])
        
        # Reformat data
        F = np.array(F).astype(np.float64)
        K = np.array(K).astype(np.uint16)
        N = np.array(N).astype(np.uint16)
        T = np.concatenate(T).astype(np.uint32)
        N_sources = F.size
        N_times = T.size
    
    # Write data
    with h5py.File(f'/mnt/extraspace/GaiaSelectionFunction/NewData/gaiaedr3_selection_function_files_h5/{mag_idx}.h5', 'w') as h: 
        
        h.create_dataset('f', data = F, compression = "lzf", chunks = True, shape = (N_sources,), dtype = np.float64, fletcher32 = False, shuffle = True)
        h.create_dataset('k', data = K, compression = "lzf", chunks = True, shape = (N_sources,), dtype = np.uint16, fletcher32 = False, shuffle = True, scaleoffset=0)
        h.create_dataset('n', data = N, compression = "lzf", chunks = True, shape = (N_sources,), dtype = np.uint16, fletcher32 = False, shuffle = True, scaleoffset=0)
        h.create_dataset('t', data = T, compression = "lzf", chunks = True, shape = (N_times, ), dtype = np.uint32, fletcher32 = False, shuffle = True, scaleoffset=0)
        
    return mag_idx

pool = multiprocessing.Pool(N_core)
for result in tqdm.tqdm(pool.imap(mp_worker, range(N_mag)),total=N_mag):
    continue