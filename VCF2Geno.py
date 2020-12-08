#!/usr/bin/python -tt

"""
Script utilizes pandas dataframes to parse a VCF file and filter individuals and markers based on user provided filters.
Logic checks are used instead of Unit Tests as the logfile is sufficient to check logic the corresponding method.
"""

import sys
import getopt
import os
from collections import Counter
import pandas as pd


__author__ = "Jaime Sheridan"
__copyright__ = "Copyright 2020"
__credits__ = ["Jaime Sheridan"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Jaime Sheridan"
__email__ = "sheridanjl75@gmail.com"
__status__ = "Beta-Testing"


def importVCF(vcfFile):
    """
    Reads in file in VCF format (Variant Call Format) and creates a dataframe.
    """
    with open(vcfFile, 'r') as inputText:
        for line in inputText:
            if line.startswith('##'):
                continue
            elif line.startswith('#CHROM'):
                vcf_header = line.strip().split('\t')
                vcf_header[0] = vcf_header[0][1:]
                vcfDF = pd.read_csv(inputText, 
                                    delimiter='\t',
                                    comment='#',
                                    names=vcf_header)
        #print(vcf_header)                                
    return(vcfDF)

def processDataFrame(vcfDF, FilterStep=0, outputFile=''):
    """
    Pre-processes the vcf dataframe in multiple passes
    First pass it re-creates a marker ID column, removes in-variant marker sites, and removes multi-allele SNPs followed by frequency calculations for rows and columns.
    Genotype calls are re-mapped to continuous "0"= homozygous reference allele, "1"= heterozygous, and "2" = homozygous alternate allele.  This continuous genotype format is used in downstream genomic modeling workflows.
    Subsequent passes do frequency calculations for rows and columns and performs genotype re-mapping.
    """     
    if FilterStep==0:

        #drop column 'ID'
        vcfDF.drop(['ID'],axis=1, inplace=True)     

        #remove non-variable sites
        vcfDF[vcfDF.ALT != '.']

        #remove multi-alleleic sites
        vcfDF[vcfDF.ALT.str.len() < 2]
        
        #add new column 'ID' which is concatenation of CHROM and POS
        #moving last column to correct position in dataframe
        vcfDF.insert(loc=2, column="ID", value= vcfDF.CHROM.astype(str).str.cat(vcfDF.POS.astype(str),sep='_').replace('\n',''))
    if FilterStep=='final':
        for column in vcfDF.iloc[:,9:]:
            vcfDF[column]=vcfDF[column].str.split(':').str.get(0).map({'0/0': '0', '1/0': '1', '0/1': '1',
               '1/1': '2', './.': '.'})
        columns2drop=[0,1,3,4,5,6,7,8]
        vcfDF.drop(vcfDF.columns[columns2drop],axis=1, inplace=True)
        vcfDF.set_index('ID', inplace=True)
        vcfDF.to_csv(outputFile)
        return None
    individualDict={}
    #individualDict is a dictionary of key being individual column ids from vcfDF and value being a dictionary of key=genotypes and value=frequency
    #explain why continuous
    for column in vcfDF.iloc[:,9:]:
        individualDict[column]= vcfDF[column].str.split(':').str.get(0).map({'0/0': '0', '1/0': '1', '0/1': '1',
               '1/1': '2', './.': '.'}).value_counts(normalize=True).to_dict()

    #markerDict is a dictionary of key being marker "ID" column from vdfDF and value being a dictionary of key=genotypes and value=frequency
    markerDict={}
    for i, row in vcfDF.iterrows():     
        markerDict[i]= row[9:].str.split(':').str.get(0).map({'0/0': '0', '1/0': '1', '0/1': '1',
               '1/1': '2', './.': '.'}).value_counts(normalize=True).to_dict()

    return individualDict, markerDict

def filterIndividuals(individualDict, newVCFdf, filters, log_file):
    """
    Filters individuals (columns) given a missing data allowance BUT is deprecated for filterMissing method.
    """ 
    fail_counter=0
    log_file.write("Failed Individuals\n")
    for individualID, frequencyDict in individualDict.items():
        if frequencyDict['.'] > 0.30:
            newVCFdf.drop([individualID],axis=1)
            fail_counter+=1
            individualMissingStats="{}\t{}\n".format(individualID, frequencyDict['.'])
            log_file.write(individualMissingStats)

    log_file.write("\nFailed Individuals Percent: {:.2f}\n".format(fail_counter/len(individualDict)*100))    
    print("Failed Individuals Percent: {:.2f}\n".format(fail_counter/len(individualDict)*100))
    log_file.flush()
    return None


def filterMissingMarkers(markerDict, newVCFdf, filters, log_file):
    """
    Filters markers (rows) given a missing data allowance BUT is deprecated for filterMissing method.
    """ 
    fail_counter=0
    #iterates through rows of frequencyDict and markerDict and filters out markers based on filters
    for i, frequencyDict in markerDict.items():
        #missing marker data allowance
        if frequencyDict['.'] > 0.3:
            newVCFdf.drop([i],axis=0)
            fail_counter+=1
    
    log_file.write("\nFailed Markers Percent: {:.2f}\n".format(fail_counter/len(markerDict)*100))    
    print("Failed Markers Percent: {:.2f}\n".format(fail_counter/len(markerDict)*100))
    log_file.flush()
    log_file.close()
    return None        


def filterMissing(vcfDict, newVCFdf, filters, log_file, filterType):
    """
    Filters rows (markers) or columns (individuals) for missing data based on given filters.
    """ 
    #logic check
    print("Pre-filter: {}".format(newVCFdf.shape))
    
    axis_variable=1
    if filterType=='markers':
        axis_variable=0
    fail_counter=0
    log_file.write("Failed {}\n".format(filterType))
    for i, frequencyDict in vcfDict.items():
        missingFreq=frequencyDict.get('.')
        if type(missingFreq)==float and missingFreq > filters:
            newVCFdf.drop([i],axis=axis_variable, inplace=True)
            fail_counter+=1
            if filterType=='individuals':
                individualMissingStats="{}\t{}\n".format(i, frequencyDict['.'])
                log_file.write(individualMissingStats)
        else:
            log_file.write("No missing {} data found for {}\n".format(filterType, i))
    log_file.write("\nFailed {} Percent: {:.2f}\n".format(filterType, fail_counter/len(vcfDict)*100))    
    print("\nFailed {} Percent: {:.2f}\n".format(filterType, fail_counter/len(vcfDict)*100))
    individualDict, markerDict=processDataFrame(newVCFdf, FilterStep=1)

    #logic check
    print("Post-filter: {}".format(newVCFdf.shape))

    log_file.flush()
    return individualDict, markerDict    

def filterHeterozygousMarkers(markerDict, newVCFdf, filters, log_file):
    """
    Filters markers based on heterozygocity.  Highly inbred crops should have heterozygocity filters allowing for less than 40% heterozygocity in a marker.
    Highly heterozygous crops such as hybrids my have more lax heterozygocity filters allowing for up to 60% heterozygocity across a marker.
    """ 
    #logic check
    print("Pre-filter: {}".format(newVCFdf.shape))

    fail_counter=0
    #iterates through rows of frequencyDict and markerDict and filters out markers based on filters
    for i, frequencyDict in markerDict.items():
        #heterozygous marker data allowance
        HetFreq=frequencyDict.get('1')
        if type(HetFreq)==float and HetFreq > filters:
            newVCFdf.drop([i],axis=0, inplace=True)
            fail_counter+=1
    individualDict, markerDict= processDataFrame(newVCFdf, FilterStep=1)
    log_file.write("\nFailed Heterozygous Markers Percent: {:.2f} of {} markers\n".format(fail_counter/len(markerDict)*100,len(markerDict)))    
    print("\nFailed Heterzygous Markers Percent: {:.2f} of {} markers\n".format(fail_counter/len(markerDict)*100, len(markerDict)))
    
    #logic check
    print("Post-filter: {}".format(newVCFdf.shape))

    return individualDict, markerDict     

def filterRareMarkers(markerDict, newVCFdf, filters, log_file):
    """
    Filters markers based on minor allele frequency (MAF).  Rare variants in crops that are highly inbred should be removed as they are most likely sequencing errors.
    A common filter would require the minor allele to appear in at least 5% of the population.  Diversity panels commonly would allow the minor allele frequency to be as low as 2%.
    """
    #logic check
    print("Pre-filter: {}".format(newVCFdf.shape))

    fail_counter=0
    #iterates through rows of frequencyDict and markerDict and filters out markers based on filters
    for i, frequencyDict in markerDict.items():
        #MAF marker data allowance
        Freq1=frequencyDict.get('1')
        Freq2=frequencyDict.get('2')
        MAF=0
        if type(Freq1)==float and type(Freq2)==float:
            MAF=Freq1+Freq2
        elif type(Freq1) == float and type(Freq2) != float:
            MAF=Freq1
        elif type(Freq1) != float and type(Freq2) == float:
            MAF=Freq2
        if MAF < filters:
            newVCFdf.drop([i],axis=0, inplace=True)
            fail_counter+=1
    individualDict, markerDict= processDataFrame(newVCFdf, FilterStep=1)
    log_file.write("\nFailed Rare Markers Percent: {:.2f} of {} markers\n".format(fail_counter/len(markerDict)*100, len(markerDict)))    
    print("\nFailed Rare Markers Percent: {:.2f} of {} markers\n".format(fail_counter/len(markerDict)*100, len(markerDict)))
    
    #logic check
    print("Post-filter: {}".format(newVCFdf.shape))

    return individualDict, markerDict     
        
def masterFilter(markerDict, individualDict, newVCFdf, filters, log_file, outputFile):
    """
    Performs ordered filtering of missing data in individuals, missing data in markers, followed by optional heterozygocity and rare variants/minor allele frequency.
    """ 
    individualDict, markerDict = filterMissing(individualDict, newVCFdf, 0.50, log_file, filterType='individual')
    individualDict, markerDict =filterMissing(markerDict, newVCFdf, filters[1], log_file, filterType='markers')
    individualDict, markerDict =filterHeterozygousMarkers(markerDict, newVCFdf, filters[2], log_file)
    individualDict, markerDict =filterRareMarkers(markerDict, newVCFdf, filters[3], log_file)
    individualDict, markerDict =filterMissing(individualDict, newVCFdf, filters[0], log_file, filterType='individual')
    processDataFrame(newVCFdf, 'final', outputFile)
    #logic check
    print("Post-filter: {}".format(newVCFdf.shape))

    return None

def main(argv):
    inputFile = ''
    outputFile = ''
    #create list of filters in the order of filter input (missing data for individual, missing data for marker, heterozygocity, minor allele frequency)
    filters = []
    prompts=["Enter in missing data maximum for an individual: ", "Missing data maximum for a marker: ", "Heterozygocity maximum (ie. inbred crop populations should have a maximum heterozygocity of 0.40): ", "Minor allele frequency (ie. inbred crop populations should have a minor allele frequency of 0.05): "]
    try:
        opts, args = getopt.getopt(
            argv, "hi:o:f:t", ["inFile=", "outFile=", "filters=", "filterType="])
    except getopt.GetoptError:
        print('VCF2geno.py -i <inFile> -o <outFile>\n\n \
		example usage: \n\n\
		python VCF2geno.py 	-i TASSELoutputFile.vcf \
							-o new.vcf \n\n \
        Optional Arguments:\n \
        -f "0.2,0.02,0.02,0.3,30.0"\n')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('vcf2geno.py -i <inFile> -o <outFile>\n')
            sys.exit()
        elif opt in ("-i", "--inFile"):
            inputFile = arg
        elif opt in ("-o", "--outFile"):
            outputFile = arg
        elif opt in ("-f", "--filters"):
            filters = arg
        elif opt in ("-t", "--filterType"):
            filters = arg
    for prompt in prompts:
        filters.append(float(input(prompt)))
    log_filePATH = '/'.join(inputFile.split('/')[:-1])
    log_file=open(log_filePATH + '/vcfLog.log', 'w')
    log_file.write("Missing individual filter: {}\nMissing marker filter: {}\nMaximum Heterozygocity filter: {}\nMinimum Minor Allele Frequency: {}\n\n".format(filters[0], filters[1], filters[2], filters[3]))

    vcfDF=importVCF(inputFile)
    individualDict, markerDict =processDataFrame(vcfDF)
    #filterIndividuals(individualDict, newVCFdf, filters, log_file)
    #filterMissingMarkers(markerDict, newVCFdf, filters, log_file)
    #filterMissing(vcfDict, newVCFdf, filters, log_file, filterType)
    #filterHeterozygousMarkers(markerDict, newVCFdf, filters, log_file)
    #filterRareMarkers(markerDict, newVCFdf, filters, log_file)
    masterFilter(markerDict, individualDict, vcfDF, filters, log_file, outputFile)

    log_file.close()

if __name__ == "__main__":
    main(sys.argv[1:])
