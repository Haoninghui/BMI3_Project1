import Read_fasta

if __name__ == '__main__':
    path = r'tests/mock_sample.txt'
    sample = Read_fasta.fasta2dict(path)
    print(sample)
