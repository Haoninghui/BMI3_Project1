# Reading .fasta file into a dictionary object, containing the chromosome index and the correspond nucleotide sequence.

def fasta2dict(filename):
    """
    :param filename: the repository of the .fasta file
    :return: a dictionary which values are the nucleotide sequences
    """
    f = open(filename, 'r')
    dic = {}
    for line in f:
        if line.startswith('>'):
            index = line.replace('>', '').split()[0]
            dic[index] = ''
        else:
            dic[index] += line.replace('\n', '').strip()
    f.close()
    return dic
