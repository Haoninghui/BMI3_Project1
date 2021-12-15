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
    parser = argparse.ArgumentParser(description='''Please enter ONE FASTA file of the consensus sequence as a query 
    sequence, and ONE FASTA file with exactly ONE chromosome sequence. \n
    The output file will be in BED format.''')
    parser.add_argument('-r', '--ref', default='./tests/chr21.fa', type=str,
                        help='input the path of the reference FASTA file')
    parser.add_argument('-q', '--query', default='./tests/LTR5_Hs.fa', type=str,
                        help='input the path of the query FASTA file')
    parser.add_argument('-p', '--path', default='./build', type=str, help='define the path of an output BED file')
    parser.add_argument('-o', '--output', default='build', type=str, help='define the name of the output BED file')
    parser.add_argument('-m', '--mismatch', default=5, type=int,
                        help='input the number of mismatches allowed during merging nearby seeds, default is 5')
    parser.add_argument('-g', '--gap', default=-5, type=int, help='input the value of gap penalty, default is -5')
    parser.add_argument('-t', '--threshold', default=500, type=float,
                        help='input the threshold of Smith–Waterman score, default is 500')
    parser.add_argument('-e', '--Escore', default=0.1, type=float, help='input the threshold E-score, default is 0.1')
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
