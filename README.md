# RNA-seq pipeline

Pipeline for taking raw FASTQ files, aligning them to a genome of choice, doing differential expression and splicing analysis.
These bash shell scripts were written to work on Compute Canada servers with slurm scheduler.
They were also designed to keep the scripts in the home/project space and the large input and output files in the scratch space for respecting storage space considerations and for faster read/write operations.

The scripts should be run from the scripts folder.

## 1. fastqc
Before
The script 01_fastq.sh will take a folder or multiple fastq.gz files as inputs and
