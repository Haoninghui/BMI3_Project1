# Final output, writing a BED file.

def write_bed(chrom, chrom_start, chrom_end, score):
    """
    Writing a BED file.
    :param chrom: The name of the chromosome (e.g. chr1, chrX) or scaffold.
    :param chrom_start: The starting position of the feature in the chromosome or scaffold.
    The number of the first base is regarded as 0.
    :param chrom_end: The end position of the feature in the chromosome or scaffold.
    :param score: The score made by SW algorithm.
    :return: NULL
    """
    f = open(r'./tests/sampleOutput.bed', 'a')
    f.write('\n{}\t{}\t{}\t{}'.format(chrom, chrom_start, chrom_end, score))
    f.close()
    # If all the records are saving to the BED file, the first line need to be delete.


def df2bed(df):
    """
    Writing BED file from dataframe.
    :param df: a dataframe containing the name of the chromosome (chrom),
    the starting/end position of the feature in the chromosome (chrom_start, chrom_end),
    and SW score.
    :return: NULL
    """
    f = open(r'./tests/sampleOutput.bed', 'w')
    for i in range(len(df)):
        chrom = df.chrom[i]
        chrom_start = df.chrom_start[i]
        chrom_end = df.chrom_end[i]
        score = df.score[i]
        if i == 0:
            f.write('{}\t{}\t{}\t{}'.format(chrom, chrom_start, chrom_end, score))
        else:
            f.write('\n{}\t{}\t{}\t{}'.format(chrom, chrom_start, chrom_end, score))
    f.close()


if __name__ == '__main__':
    import pandas as pd
    df = pd.DataFrame([['chr1', 0, 1000, 500], ['chrX', 9, 1500, 400], ['chr10', 1, 1001, 501]],
                      columns=['chrom', 'chrom_start', 'chrom_end', 'score'])
    df2bed(df)
