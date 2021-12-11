# Final output, writing a BED file.

def write_bed(chrom, chrom_start, chrom_end):
    """
    Writing a BED file.
    :param chrom: The name of the chromosome (e.g. chr1, chrX) or scaffold.
    :param chrom_start: The starting position of the feature in the chromosome or scaffold.
    The number of the first base is regarded as 0.
    :param chrom_end: The end position of the feature in the chromosome or scaffold.
    :return: NULL
    """
    f = open(r'./tests/sampleOutput.bed', 'a')
    f.write('\n{}\t{}\t{}'.format(chrom, chrom_start, chrom_end))
    f.close()
    # If all the records are saving to the BED file, the first line need to be delete.


def df2bed(df):
    """
    Writing BED file from dataframe.
    :param df: a dataframe containing the name of the chromosome (chrom),
    the starting/end position of the feature in the chromosome (chrom_start, chrom_end).
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
