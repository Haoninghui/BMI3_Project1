import Read_fasta
import Seeding
import SW_scoring
import Score_filter
import Write_BED
import numpy as np
import pandas as pd
import time
import sys
from tqdm import tqdm
import argparse
from pyfaidx import FastaNotFoundError


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='There should be some descriptions. '
                                                 'Please set your working directory to where this .py file is.')
    parser.add_argument('-r', '--ref', default='./tests/chr21.fa', type=str, help='Input the path of reference sequence')
    parser.add_argument('-q', '--query', default='./tests/LTR5_Hs.fa', type=str, help='Input the path of query sequence')
    parser.add_argument('-p', '--path', default='./build', type=str, help='Input the output path')
    parser.add_argument('-o', '--output', default='build', type=str, help='Input the name of output file')
    parser.add_argument('-m', '--mismatch', default=5, type=int,
                        help='Input the allowed mismatch during merge nearby seeds.')
    parser.add_argument('-g', '--gap', default=-5, type=int, help='Input the gap penalty')
    parser.add_argument('-t', '--threshold', default=500, type=float, help='Input the threshold of SW score.')
    parser.add_argument('-e', '--Escore', default=0.1, type=float, help='Input the threshold of Escore.')
    args = parser.parse_args()
    start = time.time()
    try:
        ref_name, ref = Read_fasta.fa2info(args.ref)
        query_name, query = Read_fasta.fa2info(args.query)
    except FastaNotFoundError:
        print('The fasta file does not exist.')
        sys.exit()
    m = len(query)
    n = len(ref)
    seed_mismatch = args.mismatch
    hamming_results = Seeding.seed_gapped(query.seq, ref.seq, seed_mismatch)  # seeding and gap free extending
    df = pd.DataFrame(columns=['chrom', 'chrom_start', 'chrom_end'])  # DataFrame to store the final output/result
    i = 0
    gap_penalty = args.gap
    for index, row in tqdm(hamming_results.iterrows()):
        ref_index = row['r']
        ref_extract = ref[row['r']:row['r']+row['rl']]
        query_extract = query[row['q']:row['q']+row['l']]
        iter_mat, dir_rec = SW_scoring.create_iterative_matrix(query_extract.seq, ref_extract.seq, gap_penalty)
        end_index = np.argwhere(iter_mat == np.max(iter_mat))[0, 1] - 1
        score = np.max(iter_mat)
        if Score_filter.SW_score_filter(score, args.threshold):  # set the cutoff to just extend high score pairs
            continue
        try:
            SW_scoring.backtracking(iter_mat, dir_rec, np.argwhere(iter_mat == np.max(iter_mat))[0, 0] - 1, end_index)
        except RecursionError or MemoryError:
            continue
        del iter_mat, dir_rec
        if Score_filter.Escore_filter(score, m, n, args.Escore):
            # the probability of this extension sequence are found randomly in database
            continue
        df.loc[i] = SW_scoring.get_index_info(ref_name, end_index, row['r'])
        i += 1
    Write_BED.makefile(args.path, args.output, df)
    end = time.time()
    print("Using time: %fs" % (end - start))
    print('--------------------------------------')
