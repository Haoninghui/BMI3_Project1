import Read_fasta
import Seeding
import SW_scoring
import Score_filter
import Write_BED
import numpy as np
import pandas as pd
import time
from tqdm import tqdm
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='There should be some descriptions.')
    parser.add_argument('-r', '--ref', type=str, help='Input the path of reference sequence')
    parser.add_argument('-q', '--query', type=str, help='Input the path of query sequence')
    parser.add_argument('-m', '--mismatch', default=5, type=int,
                        help='Input the allowed mismatch during merge nearby seeds.')
    parser.add_argument('-g', '--gap', default=-5, type=int, help='Input the gap penalty')
    parser.add_argument('-t', '--threshold', default=100, type=float, help='Input the threshold of SW score.')
    parser.add_argument('-e', '--Escore', default=0.1, type=float, help='Input the threshold of Escore.')
    args = parser.parse_args()
    start = time.time()
    ref_name, ref = Read_fasta.fa2info(args.ref)
    query_name, query = Read_fasta.fa2info(args.query)
    m = len(query)
    n = len(ref)
    seed_mismatch = args.mismatch
    hamming_results = Seeding.BLAST(query.seq, ref.seq, seed_mismatch)  # seeding and gap free extending
    df = pd.DataFrame(columns=['chrom', 'chrom_start', 'chrom_end'])  # DataFrame to store the final output/result
    i = 0
    gap_penalty = args.gap
    for index, row in tqdm(hamming_results.iterrows()):
        ref_index = row['r']
        SW_scoring.create_global()
        ref_extract = ref[row['r']:row['r']+row['rl']]
        query_extract = query[row['q']:row['q']+row['l']]
        iter_mat, dir_rec = SW_scoring.create_iterative_matrix(query_extract.seq, ref_extract.seq, gap_penalty)
        end_index = np.argwhere(iter_mat == np.max(iter_mat))[0, 1] - 1
        score = np.max(iter_mat)
        if Score_filter.SW_score_filter(score, args.threshold):  # set the cutoff to just extend high score pairs
            continue
        SW_scoring.backtracking(iter_mat, dir_rec, np.argwhere(iter_mat == np.max(iter_mat))[0, 0] - 1, end_index)
        if Score_filter.Escore_filter(score, m, n, args.Escore):
            # the probability of this extension sequence are found randomly in database
            continue
        df.loc[i] = SW_scoring.get_index_info(ref_name, end_index, row['r'])
        i += 1
    Write_BED.df2bed(df)
    end = time.time()
    print('The output BED file is ./tests/SampleOutput.BED')
    print("Using time: %fs" % (end - start))
    print('--------------------------------------')
