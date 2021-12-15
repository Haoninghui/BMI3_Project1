# Final output, writing a BED file.
import os


def makefile(path, filename: str, df):
    if os.path.exists(path):
        if os.path.isdir(path):
            f = open(path+'/'+filename+'.bed', 'w+')
            for i in range(len(df)):
                chrom = df.chrom[i]
                chrom_start = df.chrom_start[i]
                chrom_end = df.chrom_end[i]
                if i == 0:
                    f.write('{}\t{}\t{}'.format(chrom, chrom_start, chrom_end))
                else:
                    f.write('\n{}\t{}\t{}'.format(chrom, chrom_start, chrom_end))
            f.close()
            print('The output file is '+path+'/'+filename+'.bed')
        else:
            print('Please enter a directory name.')
    else:
        print('The directory does not exist.')
