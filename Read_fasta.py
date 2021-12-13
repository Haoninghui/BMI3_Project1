# Reading .fasta file into a dictionary object, containing the chromosome index and the correspond nucleotide sequence.
from pyfaidx import Fasta


def fasta2dict(filename):
    """
    :param filename: the repository of the .fasta file
    :return: a dictionary which values are the nucleotide sequences
    """
    f = open(filename, 'r')
    dic = {}
    chr_name = ""
    chr_seq = ""
    while True:
        line = f.readline()
        if not line:
            break
        if line.startswith('>'):
            if chr_name != "":
                dic[chr_name] = chr_seq
            chr_name = line.replace('>', '').split()[0]
            chr_seq = ''
        else:
            chr_seq += line.rstrip()
    dic[chr_name] = chr_seq
    f.close()
    return dic


if __name__ == '__main__':
    ref = Fasta(r'tests/families.fa')
    print(ref.keys())
    print(len(ref['DF0000558.4'][0:].seq))
    print(len(fasta2dict(r'tests/families.fa')['DF0000558.4']))
