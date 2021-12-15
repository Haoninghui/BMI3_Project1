import pandas as pd


#  Generate reverse and complementary seq
def Rev_Complementary(DNA):
    """
    Generate the reverse and complementary seq: (-) strand
    :para DNA: query seq or ref genome seq
    :return: the complementary string from 5' to 3'
    """
    com = {'A': 'T', 'a': 'T',
           'T': 'A', 't': 'A',
           'C': 'G', 'c': 'G',
           'G': 'C', 'g': 'C'}
    dna = []
    for n in DNA:
        dna.append(com[n])
    return dna


# transform base pair to number to increase the speed
def SymToNum(seq):
    """
    Transform nucleotide into number
    :para word: nucleotide seq
    :return: number list of input seq
    """
    trans = {'a': 0, 'A': 0, 'c': 1, 'C': 1, 'g': 2, 'G': 2, 't': 3, 'T': 3}
    num = []
    for n in seq:
        num.append(trans[n])
    return num


# Simply score two nucleotides using hamming distance
# After matched and merged seeds and before the SW algorithm extending
def hammingScore(q, r):
    """
    Gap free extension
    :para q,r: base pair from query seq and ref seq
    :return: the hamming distance score
    """
    if q == r:
        return 2
    else:
        return -1


# Generate seeds from query seq and store in hash table
def Seed(query):
    """
    Using sliding window to generate seeds from input seq
    :para query: input seq
    :retrun: a hash table (dict) stores 11-mer seeds and their start index
    """
    seed = dict()
    ends = len(query) - 11 + 1
    query = query.upper()
    for i in range(ends):
        q = query[i:i+11]  # slicing window
        if q in seed:  # write in python3
            seed[q].append(i)  # a seed occurs more than one time
        else:
            seed[q] = [i]
    return seed


# Transform seq into sum-value in 4-scale
# one optimization to increase the matched speed when comparing the seeds from query and ref
def SeqToNum(seq, length):
    """
    transfer sequence into a sum of SymToNum in 4-scale
    :return: unique sum of seq in 4-scale
    """
    summ = 0
    seq_num = SymToNum(seq)
    for i, n in enumerate(seq_num):
        summ += n*(4**(length-i))
    return summ


# Merge the overlap and nearby seed into one seed
def merge_seed(match, seed_mismatch):
    # match is a dataframe with start index of query and ref and length
    """
    Merge the overlapped seeds (end index > another start index) and nearby seeds (end to end) together
    :para match: Start and end index of query and ref genome and length of matched seeds (DataFrame)
    :return: merge_seeds with start and end index of query and ref genome and their length (DataFrame)
    """
    match = match.explode('q')  # split multiple index in query into single row
    match = match.reset_index(drop=True)  # reset the row index
    match = match.explode('r')  # split multiple index in ref into single row
    # sorted by query index
    match.sort_values('q', axis=0, ascending=True, inplace=True, kind='quicksort', na_position='last')
    match = match.reset_index(drop=True)  # reset the row index
    for i in range(len(match) - 1):  # previous seed
        j = i + 1
        while j < len(match):
            # distance on query and ref <= seed length + seed mismatch
            if match['q'][j]-match['q'][i] == match['r'][j]-match['r'][i] &\
                    (match['q'][j]-match['q'][i]) <= match['l'][i] + seed_mismatch:
                match['l'][i] = match['q'][j] + match['l'][j] - match['q'][i]
                match.drop(j, inplace=True)
                match = match.reset_index(drop=True)  # reset the row index
            else:
                j += 1
    return match


# Core BLAST function
# further can add seed thershold and extend thershold
def BLAST(query, ref, seed_mismatch=5):
    """
    Core BLAST function, start from query seq and ref genome
    :para query: query seq that already read-in in the main function
    :para ref: ref genome seq that already read-in in the main function
    :return: merged seed DataFrame
    """
    # matched seeds put in the dataframe
    match_seed = pd.DataFrame(columns=['q', 'r', 'l'])
    query_seed = Seed(query)
    ref_seed = Seed(ref)
    # further improvement:define seed selection thershold
    # further improvement: transfer seed seq to sum of 4-scale number
    # compare seeds from query dict and ref genome dict
    for i in query_seed:  # for every seed in query dict
        if i in ref_seed:  # write in python3
            match_seed = match_seed.append({'q': query_seed[i], 'r': ref_seed[i], 'l': 11}, ignore_index=True)
            # q: start index on query, r: start index on ref, l: length of seed on query
            # index on query and ref of matched seed store in a DF
    # Merge the overlapped and nearby seeds
    merge_seeds = merge_seed(match_seed, seed_mismatch)  # for further extend
    # release the memory, only leave merge_seed（Dataframe）
    del ref_seed, query_seed
    # First extend using simple hamming distance, gap free extension
    # generate a new column to store the length of matched seq on ref genome, rl
    merge_seeds['rl'] = 0
    # now, l is query length(ii), rl is ref length(iir)
    for row in range(len(merge_seeds)):  # for every merged seed in dataframe (1 seed/line)
        sumscore = 1 * merge_seeds['l'][row]  # Initial score
        q_start = merge_seeds['q'][row]  # start index in query seq
        r_start = merge_seeds['r'][row]  # start index in ref seq
        l = merge_seeds['l'][row]  # length of the seed
        rl = merge_seeds['l'][row]  # length of the seed in ref
        # extend on both left and right sides
        while sumscore > 0:
            if q_start > 0 and q_start + l <= len(query) - 1:  # normal case
                sumscore += hammingScore(query[q_start - 1], ref[r_start - 1]) + hammingScore(query[q_start + l],
                                                                                              ref[r_start + rl])
                q_start -= 1
                r_start -= 1
                l += 2  # 1 nucleotide on both sides(lefy and right)
                rl += 2
            elif q_start <= 0 and q_start + l <= len(query) - 1:  # out of start of query but within the end
                sumscore += hammingScore(query[q_start + l], ref[r_start + rl]) - 1
                q_start = 0
                r_start -= 1
                l += 1
                rl += 2
            elif q_start > 0 and q_start + l > len(query) - 1:  # out of end of query but within the start
                sumscore += hammingScore(query[q_start - 1], ref[r_start - 1]) - 1
                q_start -= 1
                r_start -= 1
                l += 1
                rl += 2
            elif q_start <= 0 and q_start + l > len(query) - 1:  # out of query length on both sides
                q_start = 0
                r_start -= int(sumscore/2)
                rl += sumscore
                sumscore = 0
        merge_seeds.iloc[row, 0] = q_start
        merge_seeds.iloc[row, 1] = r_start
        merge_seeds.iloc[row, 2] = l
        merge_seeds.iloc[row, 3] = rl
    # Set gap free extension threshold, only 1 in 50 seq pass
    # merge_seeds = merge_seeds[merge_seeds['l'] > query_len]
    return merge_seeds
