# Reading .fasta file into a dictionary object, containing the chromosome index and the correspond nucleotide sequence.
from pyfaidx import Fasta


def fasta2dict(filename):
    """
    :param filename: The repository of the .fasta file
    :return: A dictionary which values are the nucleotide sequences
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
    """
    Extract name and sequence from FASTA file.
    :param path: The path of the input FASTA file.
    :return: The name of the sequence, e.g. chr1; and the sequence itself, but still as FastaRecord object.
    """
    fa_rec = Fasta(path)
    name = fa_rec[0].name
    seq = fa_rec[fa_rec[0].name][0:]
    return name, seq
