import Seeding
import SW_scoring
import Write_BED
import numpy as np
import pandas as pd
import time
from tqdm import tqdm
from pyfaidx import Fasta
import math
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='There should be some descriptions.')
    parser.add_argument('-r', '--ref')
    parser.add_argument('-q', '--query')
    args = parser.parse_args()
    start = time.time()
    ref = Fasta(r'tests/chrY.fa')
    query = Fasta(r'tests/LTR5_Hs.fa')
    ref_name = ref[0].name
    query_name = query[0].name
    ref = ref[ref[0].name][0:500000]  #
    query = query[query[0].name][0:]
    m = len(query)
    n = len(ref)
    k = 0.1
    lam = 0.9  # range 0.8-0.9
    constant = k * m * n
    seed_mismatch = 5   # allowed mismatch during merge nearby seeds
    hamming_results = Seeding.BLAST(query.seq, ref.seq, seed_mismatch)  # seeding and gap free extending
    df = pd.DataFrame(columns=['chrom', 'chrom_start', 'chrom_end'])  # DataFrame to store the final output/result
    i = 0
    gap_penalty = -5  # delete 5 from score if there is a gap
    for index, row in tqdm(hamming_results.iterrows()):
        ref_index = row['r']
        SW_scoring.create_global()
        ref_extract = ref[row['r']:row['r']+row['rl']]
        query_extract = query[row['q']:row['q']+row['l']]
        iter_mat, dir_rec = SW_scoring.create_iterative_matrix(query_extract.seq, ref_extract.seq, gap_penalty)
        index1 = np.argwhere(iter_mat == np.max(iter_mat))[0, 0] - 1
        index2 = np.argwhere(iter_mat == np.max(iter_mat))[0, 1] - 1
        score = np.max(iter_mat)
        if score <= 500:  # set the cutoff to just extend high score pairs
            continue
        SW_scoring.backtracking(iter_mat, dir_rec, index1, index2)
        E = constant * math.exp(-lam * score)
        if E >= 0.1:  # the probability of this extension sequence are found randomly in database
            continue
        df.loc[i] = SW_scoring.get_index_info(ref_name, index2, row['r'])
        i += 1
    Write_BED.df2bed(df)
    end = time.time()
    print("Using time: %fs" % (end - start))
    print('--------------------------------------')


