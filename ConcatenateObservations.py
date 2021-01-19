import tqdm

N_sources = 1811709771
N_chunk = 3*17*6791
N_block = 5231
N_maxobs = 256

with open('./gaiaedr3_selection_function.csv', 'w') as f:
    #for block_idx in tqdm.tqdm(range(N_block)):
    for block_idx in tqdm.tqdm(range(3)):
        with open(f'./store/gaiaedr3_selection_function_{block_idx}.csv', 'r') as g:
            for idx,line in enumerate(g):
                if line.strip():
                    f.write(str(idx+N_chunk*block_idx)+','+line)