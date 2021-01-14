import numpy as np
import h5py
import tqdm

bins = np.array([0.0, 5.0, 16.0, 17.0, 18.0, 18.2, 18.4, 18.6, 18.8, 19.0, 19.1, 19.2, 19.3, 19.4, 19.5, 19.6, 19.7, 19.8, 19.9, 20.0, 20.1, 20.2, 20.3, 20.4, 20.5, 20.6, 20.7, 20.8, 20.9, 21.0, 21.1, 21.2, 21.3, 21.4, 21.5, 25.0])
N_bins = bins.size-1

N_sources = 1811709771
N_chunk = 3*17*6791
N_block = 5231
N_maxobs = 256
        
with h5py.File('./gaiaedr3.h5', 'r') as g, h5py.File('./gaiaedr3_times.h5', 'r') as t, open('./gaiaedr3_selection_function.csv', 'w') as f:
    
    #for block_idx in tqdm.tqdm(range(N_block)):
    for block_idx in tqdm.tqdm(range(10)):
        
        # Load in source_id's and verify that they are the same
        source_id_g = g['source_id'][block_idx*N_chunk:(block_idx+1)*N_chunk]
        source_id_t = t['source_id'][block_idx*N_chunk:(block_idx+1)*N_chunk]
        assert (source_id_g == source_id_t).all()
        
        # Load in data
        phot_g_mean_mag = g['phot_g_mean_mag'][block_idx*N_chunk:(block_idx+1)*N_chunk]
        matched_transits = g['matched_transits'][block_idx*N_chunk:(block_idx+1)*N_chunk]
        fov_1_n = t['fov_1_n'][block_idx*N_chunk:(block_idx+1)*N_chunk]
        fov_2_n = t['fov_2_n'][block_idx*N_chunk:(block_idx+1)*N_chunk]
        fov_1_times = t['fov_1_times'][block_idx*N_chunk:(block_idx+1)*N_chunk]
        fov_2_times = t['fov_2_times'][block_idx*N_chunk:(block_idx+1)*N_chunk]
        
        phot_g_mean_mag_idx = np.digitize(phot_g_mean_mag,bins)-1
        source_times = [np.sort(np.concatenate([fov_1_times[source_idx,:fov_1_n[source_idx]],fov_2_times[source_idx,:fov_2_n[source_idx]]])) for source_idx in range(N_chunk)]
        
        source_lines = [f"{phot_g_mean_mag_idx[source_idx]},{matched_transits[source_idx]},{fov_1_n[source_idx]+fov_2_n[source_idx]},"+','.join(map(str, source_times[source_idx])) for source_idx in range(N_chunk)]
        
        f.writelines("{}\n".format(line) for line in source_lines)
        
        #if block_idx == N_block-1:
        #    f.write('\n')
        
        
        


