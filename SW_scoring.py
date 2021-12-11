# Smith-Waterman Algorithm

import numpy as np
import pandas as pd


def index2seq(start: int, end: int, seq: str):
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
    :param sym: symbol, ordering as ACGT
    :return: transformed symbol
    """
    trans = str.maketrans('ACTGacgt', '01230123')
    return int(sym.translate(trans))


def get_score(sym1: str, sym2: str) -> int:
    """
    Compare the similarity of two bases via scoring matrix.
    :param sym1: first symbol to be compared
    :param sym2: second symbol yo be compared
    :return: the score of the similarity of 2 symbols
    """
    mat = pd.DataFrame([[4, -7, -5, -7],
                        [-7, 4, -7, -5],
                        [-5, -7, 4, -7],
                        [-7, -5, -7, 4]])
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


def route_rec():
    """
    Record the route during backtracking.
    :return: NULL
    """
    global finalList
    global start_index
    finalList = []
    start_index = []


def backtracking(iter_mat, dir_rec, seq1, seq2, index1, index2):
    """
    Find the route from end position to start position, in the iterative matrix.
    :param iter_mat: Iterative matrix for scoring and recording the scores.
    :param dir_rec: Record of possible directions from left/up-left/up to current position.
    :param seq1: The first sequence.
    :param seq2: The second sequence.
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
            backtracking(iter_mat, dir_rec, seq1, seq2, index1, index2-1)
            finalList.pop()
    if dir_rec[index1][index2][1] == 1:
        if all(dir_rec[index1-1][index2-1] == [0, 0, 0]):
            start_index.append(index1-1)
            start_index.append(index2-1)
            return
        else:
            finalList.append(1)
            backtracking(iter_mat, dir_rec, seq1, seq2, index1-1, index2-1)
            finalList.pop()
    if dir_rec[index1][index2][2] == 1:
        if all(dir_rec[index1-1][index2] == [0, 0, 0]):
            start_index.append(index1-1)
            start_index.append(index2)
            return
        else:
            finalList.append(2)
            backtracking(iter_mat, dir_rec, seq1, seq2, index1-1, index2)
            finalList.pop()


def get_index_info(chrom1: str, chrom2: str, start1: int, start2: int, end1: int, end2: int):
    """
    Output a dataframe with information like BED.
    :param chrom1: The name of the chromosome from which sequence1 came.
    :param chrom2: The name of the chromosome from which sequence2 came.
    :param start1: The start position of sequence1 in the chromosome.
    :param start2: The start position of sequence1 in the chromosome.
    :param end1: The end position of sequence1 in the chromosome.
    :param end2: The end position of sequence2 in the chromosome.
    :return: Dataframe with information including chrom, chrom_start, and chrom_end.
    """
    df = pd.DataFrame([[chrom1, start1, end1], [chrom2, start2, end2]],
                      columns=['chrom', 'chrom_start', 'chrom_end'])
    return df


if __name__ == '__main__':
    route_rec()
    seq1 = 'ACGTA'
    seq2 = 'GGGGGGGACGT'
    gap = -5
    iter_mat, dir_rec = create_iterative_matrix(seq1, seq2, gap)
    index1 = np.argwhere(iter_mat == np.max(iter_mat))[0, 0]-1
    index2 = np.argwhere(iter_mat == np.max(iter_mat))[0, 1]-1
    backtracking(iter_mat, dir_rec, seq1, seq2, index1, index2)
    print(get_index_info('chrX', 'chrY', start_index[0], start_index[1], index1, index2))
