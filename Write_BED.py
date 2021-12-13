# Final output, writing a BED file.


def df2bed(df):
    """
    Writing BED file from dataframe.
    :param df: a dataframe containing the name of the chromosome (chrom),
    and the starting/end position of the feature in the chromosome (chrom_start, chrom_end).
    :return: NULL
    """
    f = open(r'./tests/sampleOutput.bed', 'w')
    for i in range(len(df)):
        chrom = df.chrom[i]
        chrom_start = df.chrom_start[i]
        chrom_end = df.chrom_end[i]
        if i == 0:
            f.write('{}\t{}\t{}'.format(chrom, chrom_start, chrom_end))
        else:
            f.write('\n{}\t{}\t{}'.format(chrom, chrom_start, chrom_end))
    f.close()
