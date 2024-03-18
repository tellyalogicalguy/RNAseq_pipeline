# RNA-seq pipeline

Pipeline for taking raw FASTQ files, aligning them to a genome of choice, doing differential expression and splicing analysis.

These bash shell scripts were written to work on *Compute Canada* HPC servers with the **SLURM** scheduler.  
They were also designed to keep the scripts in the `HOME`/`PROJECT` space and the large input and output files in the `SCRATCH` space for not running into storage limitaions and for faster read/write operations.  
This is done by mirroring the folder structure in both the `HOME` and `SCRATCH` spaces, with only symlinks pointing to the large files being stored in the `HOME` space. 

## Setting up
The project folders and scripts can be set-up by running the `setup_initial_folders_hnRNPL_RNAseq.sh` script.  
This script will take the name of a cell line/type as an argument and set up several folders necessary for the rest of the pipeline.
Since this pipeline was set-up for analysing publicly available RNA-seq data of different cell types where a particular protein of interest (hnRNPL) had been depleted, the default folder structure currently is `./exampleCellLine/rna_seq/hnrnpl/` in which the `data`, `scripts` and `results` folders will be made.  
Here is an example output when running `setup_initial_folders_hnRNPL_RNAseq.sh exampleCellLine`, which contains useful instructions to finish set-up.
```
Setting up folders in /home/subrampg/binf_analyses/hnrnpl_project/ (local)
/home/subrampg/binf_analyses/hnrnpl_project/exampleCellLine
└── rna_seq
    └── hnrnpl
        ├── data
        │   ├── data_info.txt
        │   ├── initial_raw_fastq_setup_commands.txt
        │   ├── raw_fastq
        │   ├── sra_download.sh
        │   └── sra_list.txt
        └── scripts
            ├── 01_fastq.sh
            ├── 12_trimAlign.sh
            ├── 12_trimAlign_sbatch.sh
            ├── 12_trimAlign_sbatch_SE.sh
            ├── 34_star2ndPassBW_sbatch.sh
            ├── 34_star2ndPassBW_sbatch_SE.sh
            ├── 56_fcDESeq_bwQC.sh
            ├── 56_fcDESeq_bwQC_sbatch.sh
            ├── 56_fcDESeq_bwQC_sbatch_unstranded.sh
            ├── 7_rMATS_sbatch.sh
            ├── deepTools_requirements.txt
            ├── featureCounts_DESeq.R
            ├── getSirt1rMATSevent.R
            ├── multiqc_requirements.txt
            ├── rMATS2sashimi_requirements.txt
            ├── seff_outputs
            └── slurm_outputs

7 directories, 19 files
Setting up folders in /lustre07/scratch/subrampg/binf_analyses_data/hnrnpl_project/ (scratch)
/lustre07/scratch/subrampg/binf_analyses_data/hnrnpl_project/exampleCellLine
└── rna_seq
    └── hnrnpl
        └── data
            └── raw_fastq

4 directories, 0 files
Transfer raw FASTQ files with rsync from local computer to:
subrampg@narval.computecanada.ca:/lustre07/scratch/subrampg/binf_analyses_data/hnrnpl_project//exampleCellLine/rna_seq/hnrnpl/data/raw_fastq/

Then do commands in ./exampleCellLine/rna_seq/hnrnpl/data/initial_raw_fastq_setup_commands.txt, i.e.:

After folders are set up and files were transferred with rsync:

------------------------------------------------------------------------------
!!!MODIFY "data_info.txt" IN "data" folder to reflect file and sample names!!!
------------------------------------------------------------------------------

Go to the scratch "raw_fastq" folder, run:
for file in $(ls *.gz); do md5sum $file >> md5sums.txt; done
to do md5 checksums.

Go to the home "raw_fastq" folder, run to symlink files:
ln -s ~/scratch/binf_analyses_data/"??PREFIX??"/"??FACTOR??"/"??scratch_data_folder??/raw_fastq/* .
```

The default folder structures can be changed by altering the `setup_initial_folders_hnRNPL_RNAseq.sh` folder.

## Downloading files from SRA (optional)
You can chose to download files from the SRA.  
Update the `sra_list.txt` in the `./exampleCellLine/rna_seq/hnrnpl/data/` folder and then running `./sra_download.sh`. This script will downloaded files from SRA to the `./exampleCellLine/rna_seq/hnrnpl/data/raw_fastq/` folder in the SCRATCH space and symlink the raw files to the mirrored folder in HOME space.  
Alternatively, the raw files can be uploaded to the SCRATCH space directly. In this case, please symlink the files by `ln -s`.  

Do `md5sum` to check file integrity.

##

This completes the setup.

Once the files are setup, scripts can be run from the `./scripts` folder.

#### fastqc
Before running the entire pipeline, the quality of the raw files can be checked using fastqc.  
The script `01_fastq.sh` will take multiple `*.fastq.gz` files or a folder containing said files as inputs.  
The output will be in `../results/fastqc/` folder.

## Running the pipeline
The pipeline can be starting by running `./12_trimAlign.sh`.  
This will take you through a series of prompts 
