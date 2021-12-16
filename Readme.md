LTR Predictor
=====
LTR predictor is a pure Python3-based BLAST-like endogenous retrovirus (ERV) alignment tool, helping you identify the chromosomal coordinates of the input consensus sequence. The LTR predictor can be run in Linux or Windows command line.

You can enter **ONE FASTA** file of the consensus sequence as a query sequence, and **ONE FASTA** file with **exactly ONE** chromosome sequence. (Please note that it may take ~2 hours to get the alignment result of human chromosome 1.)

The output file will be in **BED format** with the **NAME** of the chromosome, **START & END POSITIONS** of the feature in chromosomal coordinates.


![Alignment result of human chromosome 1](https://github.com/Haoninghui/BMI3_Project1/blob/master/docs/chr1.svg)


### Documentation
* Operating systems: Linux, Windows
* Requirements: Python (3.6+)
* Supported Format: FASTA (.fa)


Quickstart Guide
------
### Development
Load `main.py` in your favorite Python development platform.

### Execute
```shell
python main.py
```

#### Basic usage
##### Optional Arguments:
-h, --help  
show this help message and exit

-r REF, --ref REF   
input the path of the reference FASTA file, default is ./tests/chr21.fa

-q QUERY, --query QUERY 
input the path of the query FASTA file, default is ./tests/LTR5_Hs.fa

-p PATH, --path PATH  
define the path of an output BED file, default is ./build

-o OUTPUT, --output OUTPUT  
define the name of the output BED file, default is build

-m MISMATCH, --mismatch MISMATCH    
input the number of mismatches allowed during merging nearby seeds, default is 5

-g GAP, --gap GAP   
input the value of gap penalty, default is -5

-t THRESHOLD, --threshold THRESHOLD 
input the threshold of Smith¨CWaterman score, default is 500

-e ESCORE, --Escore ESCORE  
input the threshold E-score, default is 0.1

##### Visualization
The output BED file can be imported into a genomics data viewer to visualize, for instance, [IGV](http://www.igv.org), as shown above.

Acknowledgments
-----
LTR predictor uses the following third-party libraries:
- numpy for large matrices storage and processing 
- pandas for dataframe storage and processing
- argparse for creating command-line interfaces
- tqdm for creating progress bar
- pyfaidx for reading FASTA files
