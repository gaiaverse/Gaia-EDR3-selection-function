import tqdm
from contextlib import ExitStack

N_sources = 1811709771
N_chunk = 3*17*6791
N_block = 5231
N_maxobs = 256
N_mag = 35+1 # The +1 is to give a bin for the stars without G magnitudes to sit
star_count = {str(mag_idx):0 for mag_idx in range(N_mag)}

with ExitStack() as stack:
    file_box = {str(mag_idx):stack.enter_context(open(f'./gaiaedr3_selection_function_files/{mag_idx}.csv', 'a')) for mag_idx in range(N_mag)}
    
    for block_idx in tqdm.tqdm(range(N_block)):
    #for block_idx in tqdm.tqdm(range(3)):
        with open(f'./store/gaiaedr3_selection_function_{block_idx}.csv', 'r') as g:
            for idx,line in enumerate(g):
                if line.strip():
                    mag_str,rest_of_line = line.split(',', 1)
                    file_box[mag_str].write(rest_of_line)
                    star_count[mag_str] += 1