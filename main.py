import BLAST_new
import SW_scoring
import Write_BED
import numpy as np
import pandas as pd
import time
from tqdm import tqdm
from pyfaidx import Fasta


if __name__ == '__main__':
    start = time.time()
    ref = Fasta(r'tests/chrY.fa')
    query = Fasta(r'tests/families.fa')
    ref = ref['chrY'][0:]
    query = query['DF0000558.4'][0:]
    #query = {'DF0000558.4':'GGGCTAGCTACGTCCCAGGTGGCCGGTACGATCGGGCTAGCATGCCGTACG'}
    #ref = {'chrY': 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCGGGGGTAGGGATCAAACGTCCCAGGTGGGCTAGCTACGTCCCGGGATCCCTACATGCCGTACGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'}
    hamming_results = BLAST_new.BLAST(query.seq, ref.seq)
    df = pd.DataFrame(columns=['chrom', 'chrom_start', 'chrom_end'])
    for index, row in tqdm(hamming_results.iterrows()):
        ref_index = row['r']
        SW_scoring.create_global()
        '''
        ref_extract = SW_scoring.index2seq(row['r'], row['r']+row['rl']-1, ref.seq)
        query_extract = SW_scoring.index2seq(row['q'], row['q']+row['l']-1, query.seq)
        '''
        ref_extract = ref[row['r']:row['r']+row['rl']]
        query_extract = query[row['q']:row['q']+row['l']]
        iter_mat, dir_rec = SW_scoring.create_iterative_matrix(query_extract.seq, ref_extract.seq, -5)
        index1 = np.argwhere(iter_mat == np.max(iter_mat))[0, 0] - 1
        index2 = np.argwhere(iter_mat == np.max(iter_mat))[0, 1] - 1
        score = np.max(iter_mat)
        if score <= 500:
            continue
        SW_scoring.backtracking(iter_mat, dir_rec, index1, index2)
        df.loc[index] = SW_scoring.get_index_info('chrY', index2, row['r'])
    Write_BED.df2bed(df)
    end = time.time()
    print("Using time: %fs" % (end - start))
    print('--------------------------------------')

