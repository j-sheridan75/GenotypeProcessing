#!/bin/bash

#shell script to run the TASSEL GBSv2 pipeline and filter vcf output file

## set variables
SEQUENCE=".data/fastq/"
FASTA="./data/genomes/test_genome.faa" 
GENOME="./data/indexes/test_genome" 
DATABASE="./data/intermediate/test_DB.db"
KEY="./data/intermediate/KEYFILE_test.txt" 
TAGS="./data/intermediate/testTags.fastq"
SAM="./data/intermediate/testTagsALL.sam"
SAMhq= "./data/intermediate/testTags_hq.sam"
SAMlq= "./data/intermediate/testTags_lq.sam"
SAMunal= "./data/intermediate/testTags_UNAL.sam"
STATS="./data/intermediate/testTags_stats.txt"
QSCORE="./data/intermediate/testTags_qual.txt"
H5="./data/intermediate/test_H5.h5"
VCF="./data/intermediate/test_VCF.vcf"
H5CLOSED="./data/intermediate/test_H5CLOSED.h5"
OUT="./data/GBSv2pipeline.out"


## SNP Calling on Genotype-by-sequencing data using TASSEL GBSv2 pipeline
## TASSEL5 GBSSeqToTagDBPlugin: RUN Tags to DB  
perl /home/solutions/tassel-5-standalone/run_pipeline.pl -Xms10G -Xmx10G -GBSSeqToTagDBPlugin -e ApeKI -i ${SEQUENCE} -db ${DATABASE} -k ${KEY} -mnQS 20 -mxKmerNum 300000000 >> ${OUT}

## TASSEL5 TagExportToFastqPlugin: EXPORT Tags to FASTQ format
perl run_pipeline.pl -Xms10G -Xmx10G -TagExportToFastqPlugin -db ${DATABASE} -o ${TAGS} -c 20 >> ${OUT}

## BOWTIE2 to index genome if not already done
bowtie2-build-l --threads 8 ${FASTA} ${GENOME}

## BOWTIE2 to align TAGS to genome indexes
bowtie2 -p 20 --very-sensitive-local -x ${GENOME} -U ${TAGS} -S ${SAM}

## SAMTOOLS to separate  (alignments of tags mapping to one location on genome),  (alignments of tags mapping to multiple locations on genome), and UNALIGNED tags)
samtools view -F 4 ${SAM} | awk '$0~"XS:"' > ${SAMlq}

samtools view -F 4 ${SAM} | grep -v "XS:" > ${SAMhq}

samtools view -f 4 ${SAM} > ${SAMunal}

##TASSEL5 SAMToGBSdbPlugin: Stores position information from aligned tag in SAM file
perl /home/solutions/tassel-5-standalone/run_pipeline.pl -Xms10G -Xmx10G -SAMToGBSdbPlugin -i ${SAM} -db ${DATABASE} -aProp 0.96 -aLen 20 >> ${OUT}
 
##  TASSEL5 DiscoverySNPCallerPluginV2: Identify SNPs from aligned TAGS
perl /home/solutions/tassel-5-standalone/run_pipeline.pl -Xms10G -Xmx10G -DiscoverySNPCallerPluginV2 -db ${DATABASE} -deleteOldData true -mnLCov 0.50 -mnMAF 0.02 >> ${OUT}

##  TASSEL5 SNPQualityProfilerPlugin: Scores SNPs for coverage, depth, and stats
perl /home/solutions/tassel-5-standalone/run_pipeline.pl -Xms10G -Xmx10G -SNPQualityProfilerPlugin -db ${DATABASE} -statFile ${STATS} >> ${OUT}

##  TASSEL5 UpdateSNPPositionQUalityPlugin: optional to store quality scores for SNP pisitions (can be created from STAT file)
perl /home/solutions/tassel-5-standalone/run_pipeline.pl -Xms10G -Xmx10G -UpdateSNPPositionQualityPlugin -db ${DATABASE} -qsFile ${QSCORE} >> ${OUT}
 
##  TASSEL5 ProductionSNPCallerPluginV2: adds genotype calls to database
perl /home/solutions/tassel-5-standalone/run_pipeline.pl -Xms10G -Xmx10G -ProductionSNPCallerPluginV2 -db ${DATABASE} -e ApeKI -i ${SEQUENCE} -k ${KEY} -minPosQS 0 -ko true -o ${H5} >> ${OUT}

##  TASSEL5 BuildUnfinishedHDF5GenotypesPlugin: Closes h5 file
perl /home/solutions/tassel-5-standalone/run_pipeline.pl -Xms10G -Xmx10G -BuildUnfinishedHDF5GenotypesPlugin -i ${H5} -o ${H5CLOSED}

##  TASSEL5 Export as vcf
perl /home/solutions/tassel-5-standalone/run_pipeline.pl -h5 ${H5CLOSED} -export ${VCF} -exportType VCF >> ${OUT}

## python script to parse VCF and filter individuals and markers
time python3 ./VCF2Geno.py -i ${VCF} -o ./TESTout.vcf