# Smith-Waterman Algorithm

import numpy as np
from pyfaidx import Fasta


def create_global():
    """
    Create some global variables.
    :return: NULL
    """
    global finalList
    global start_index
    global mat
    finalList = []
    start_index = []
    mat = np.array([[10, -3, -1, -4, 10],
                    [-3, 9, -5, 0, -3],
                    [-1, -5, 7, -3, -1],
                    [-4, 0, -3, 8, -4],
                    [10, -3, -1, -4, 10]])


def index2seq(start, end, seq):
    """
    Transform indexes to sequence.
    :param start: the starting position of a subsequence in a certain sequence
    :param end: the end position of a subsequence in a certain sequence
    :param seq: the whole-length sequence
    :return: the subsequence
    """
    return seq[start:end+1]


def sym2no(sym: str) -> int:
    """
    Transform symbols to numbers.
    :param sym: symbol, ordering as ACGT, N = A
    :return: transformed symbol
    """
    trans = str.maketrans('ACTGNacgtn', '0123001230')
    return int(sym.translate(trans))


def get_score(sym1: str, sym2: str) -> int:
    """
    Compare the similarity of two bases via scoring matrix.
    :param sym1: first symbol to be compared
    :param sym2: second symbol yo be compared
    :return: the score of the similarity of 2 symbols
    """
    global mat
    return mat[sym2no(sym1)][sym2no(sym2)]


def create_iterative_matrix(seq1: str, seq2: str, gap=-5):
    """
    Create the iterative matrix, thus getting the highest score.
    :param seq1: the first sequence.
    :param seq2: the second sequence.
    :param gap: Gap penalty.
    :return: A iterative matrix with score of each pair of bases comparison,
    and a record of possible directions to reach the highest score.
    """
    length1 = len(seq1)
    length2 = len(seq2)
    iter_mat = np.zeros((length1+1, length2+1))
    dir_rec = np.zeros((length1+1, length2+1, 3))
    for i in range(1, length1+1):
        for j in range(1, length2+1):
            neighbors = [iter_mat[i, j-1], iter_mat[i-1, j-1], iter_mat[i-1, j]]
            neighboring_scores = [gap, get_score(seq1[i-1], seq2[j-1]), gap]
            temp_score = np.add(neighbors, neighboring_scores)
            if max(temp_score) > 0:
                iter_mat[i, j] = max(temp_score)
                possible_dirs = [y for y, x in enumerate(temp_score) if x == max(temp_score)]
                for k in range(len(possible_dirs)):
                    dir_rec[i][j][possible_dirs[k]] = 1
    return iter_mat, dir_rec


def backtracking(iter_mat, dir_rec, index1, index2):
    """
    Find the route from end position to start position, in the iterative matrix.
    :param iter_mat: Iterative matrix for scoring and recording the scores.
    :param dir_rec: Record of possible directions from left/up-left/up to current position.
    :param index1: Initialize with end index of sequence1.
    :param index2: Initialize with end index of sequence2.
    :return: NULL
    """
    if dir_rec[index1][index2][0] == 1:
        if all(dir_rec[index1][index2-1] == [0, 0, 0]):
            start_index.append(index1)
            start_index.append(index2-1)
            return
        else:
            finalList.append(0)
            backtracking(iter_mat, dir_rec, index1, index2-1)
            finalList.pop()
    if dir_rec[index1][index2][1] == 1:
        if all(dir_rec[index1-1][index2-1] == [0, 0, 0]):
            start_index.append(index1-1)
            start_index.append(index2-1)
            return
        else:
            finalList.append(1)
            backtracking(iter_mat, dir_rec, index1-1, index2-1)
            finalList.pop()
    if dir_rec[index1][index2][2] == 1:
        if all(dir_rec[index1-1][index2] == [0, 0, 0]):
            start_index.append(index1-1)
            start_index.append(index2)
            return
        else:
            finalList.append(2)
            backtracking(iter_mat, dir_rec, index1-1, index2)
            finalList.pop()


def get_index_info(chrom: str, chrom_end: int, r: int):
    """
    Output a dataframe with information like BED.
    :param chrom: The name of the chromosome from which sequence1 came.
    :param chrom_end: The end position of sequence1 in the chromosome.
    :param r: Index of reference genome/chromosome.
    :return: A list with information including chrom, chrom_start, chrom_end, and score.
    """
    global start_index
    df = [chrom, start_index[1]+r, chrom_end+r]
    return df
