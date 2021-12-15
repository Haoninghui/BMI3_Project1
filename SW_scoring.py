# Smith-Waterman Algorithm, and an output function of dataframe

import numpy as np

# define some global variables
finalList = []  # temporarily storing the path in backtracking.
start_index = []
mat = np.array([[10, -3, -1, -4, 10],
                [-3, 9, -5, 0, -3],
                [-1, -5, 7, -3, -1],
                [-4, 0, -3, 8, -4],
                [10, -3, -1, -4, 10]])  # scoring matrix for single-bases comparison


def index2seq(start: int, end: int, seq):
    """
    Extract a sub-sequence from a sequence.
    :param start: The starting position of a subsequence in a certain sequence
    :param end: The end position of a subsequence in a certain sequence
    :param seq: The whole-length sequence, a FastaRecord object.
    :return: the subsequence
    """
    return seq[start:end+1]


def sym2no(sym: str) -> int:
    """
    Transform symbols to numbers, easier for numpy matrix to manipulate.
    :param sym: Symbol, ordering as ACGTN, N = A.
    :return: Transformed symbol.
    """
    trans = str.maketrans('ACTGNacgtn', '0123001230')
    return int(sym.translate(trans))


def get_score(sym1: str, sym2: str) -> int:
    """
    Compare the similarity of two bases via scoring matrix.
    :param sym1: First symbol to be compared.
    :param sym2: Second symbol to be compared.
    :return: The score of the similarity of 2 symbols.
    """
    global mat
    return mat[sym2no(sym1)][sym2no(sym2)]


def create_iterative_matrix(seq1: str, seq2: str, gap=-5):
    """
    Create an iterative matrix, thus getting the highest score.
    :param seq1: The first sequence.
    :param seq2: The second sequence (ref).
    :param gap: Gap penalty.
    :return: An iterative matrix with score of each pair of bases comparison,
    and a record of possible directions to reach the highest score.
    """
    global finalList
    global start_index
    finalList = []
    start_index = []
    length1 = len(seq1)
    length2 = len(seq2)
    iter_mat = np.zeros((length1+1, length2+1))
    dir_rec = np.zeros((length1+1, length2+1, 3))
    for i in range(1, length1+1):
        for j in range(1, length2+1):
            neighbors = [iter_mat[i, j-1], iter_mat[i-1, j-1], iter_mat[i-1, j]]
            neighboring_scores = [gap, get_score(seq1[i-1], seq2[j-1]), gap]
            # Apart from the first row & column, the score in each cell is calculated by the neighboring cells:
            # left, up-left, and up. Gap penalty is also considered.
            temp_score = np.add(neighbors, neighboring_scores)
            # The score of the current cell may be come from multiple directions,
            # judge and assign the values.
            if max(temp_score) > 0:
                iter_mat[i, j] = max(temp_score)
                possible_dirs = [y for y, x in enumerate(temp_score) if x == max(temp_score)]
                for k in range(len(possible_dirs)):
                    dir_rec[i][j][possible_dirs[k]] = 1
    return iter_mat, dir_rec


def backtracking(iter_mat, dir_rec, index1: int, index2: int):
    """
    Find the route from end position to start position w.r.t. the iterative matrix.
    :param iter_mat: Iterative matrix for scoring and recording the scores.
    :param dir_rec: Record of possible directions from left|up-left|up to current cell.
    :param index1: Initialize with end index of sequence1.
    :param index2: Initialize with end index of sequence2 (ref).
    :return: NULL
    """
    if dir_rec[index1][index2][1] == 1:  # First, we like the the up-left direction representing a match.
        # If the position points to the up-left, judge whether its neighbors are all zero (end backtracking).
        if np.all(dir_rec[index1-1][index2-1] == [0, 0, 0]):
            start_index.append(index1-1)
            start_index.append(index2-1)
            return  # End backtracking.
        else:
            finalList.append(1)  # Up-left.
            backtracking(iter_mat, dir_rec, index1-1, index2-1)
            finalList.pop()  # Take the last element.
    # ibid.
    elif dir_rec[index1][index2][0] == 1:
        if np.all(dir_rec[index1][index2-1] == [0, 0, 0]):
            start_index.append(index1)
            start_index.append(index2-1)
            return
        else:
            finalList.append(0)  # Left.
            backtracking(iter_mat, dir_rec, index1, index2-1)
            finalList.pop()
    # ibid.
    else:
        if np.all(dir_rec[index1-1][index2] == [0, 0, 0]):
            start_index.append(index1-1)
            start_index.append(index2)
            return
        else:
            finalList.append(2)  # Up.
            backtracking(iter_mat, dir_rec, index1-1, index2)
            finalList.pop()


def get_index_info(chrom: str, chrom_end: int, r: int):
    """
    Export a dataframe with information like BED.
    :param chrom: The name of the chromosome from which sequence came.
    :param chrom_end: The end position of sequence1 in the chromosome.
    :param r: Index of reference genome/chromosome.
    :return: A list with information including chrom, chrom_start, chrom_end, and score.
    """
    global start_index
    df = [chrom, start_index[1]+r, chrom_end+r]
    return df
