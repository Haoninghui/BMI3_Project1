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


def fa2info(path):
    fa_rec = Fasta(path)
    name = fa_rec[0].name
    seq = fa_rec[fa_rec[0].name][0:]
    return name, seq
