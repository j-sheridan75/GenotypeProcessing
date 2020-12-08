# GenotypeProcessing

Basic workflow for calling variants from genotype-by-sequencing (GBS) in plants (https://doi.org/10.1371/journal.pone.0019379).

GBSv2.sh is an example of shell scripting and cannot be run as it requires software dependencies as well as large compressed FASTQ files.  To run GBSv2.sh, the 40 GB supporting FASTQ file is also required.  The last step of GBSv2.sh, VCF2geno.py can be run outside of GBSv2.sh.

GBSv2.sh utilizes the below software dependencies:
Bowtie2
SAMTools
TASSELv5standalone ->https://bitbucket.org/tasseladmin/tassel-5-source/wiki/https://bitbucket.org/tasseladmin/tassel-5-source/wiki/Tassel5GBSv2Pipeline
perl
python3
