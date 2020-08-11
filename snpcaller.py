import os
import argparse
import glob
from Bio import SeqIO
from itertools import combinations
from functools import partial
import multiprocessing
import subprocess


# Calculate read coverage of positions in the genome
def calcCoverage(refSeq, inBam, baseName):
    """Uses bcftools mpileup to read coverage of positions in genome

    Params
    ------
    refSeq: String
         File location of the reference multifasta file (already indexed)
    inBam: String
         File location of the input bam file
    baseName: String
         Basename of the isolate 

    Returns
        -------
        bcfFile: String
             File location of the output bcf file

    """
    bamFile = baseName+'.bam'
    bcfFile = baseName+'.bcf'
    
    try:
        r1 = subprocess.check_output(['bcftools', 'mpileup', '-O','b','-o', bcfFile,'-f', refSeq, bamFile])
        # options: -O b to generate bcf format output, -f flags path to reference genome
    except subprocess.CalledProcessError as error:
        errorMsg = ('Status : Fail',error.returncode,error.output.strip('\n'))
    
    return(bcfFile)


# Filter and report the SNP variants in variant calling format
def callVariants(inBcf, baseName):

    """Uses bcftools call to detect SNPs

    Params
    ------
    refSeq: String
         File location of the reference multifasta file (already indexed)
    inBcf: String
         File location of the input bcf file
    baseName: String
         Basename of the isolate 

    Returns
    -------
    vcfFile: String
         File location of the output vcf file

    """
    bcfFile = baseName+'.bcf'
    vcfFile = baseName+'.vcf'
    
    try:
        r1 = subprocess.check_output(['bcftools', 'call', '-ploidy','1', '-m','-v','-o', vcfFile, bcfFile])
        # options: -m for multiallelic & rare-variant calling, -v to output variant sites only, ploidy is 1 for haploid genome 
    except subprocess.CalledProcessError as error:
        errorMsg = ('Status : Fail',error.returncode,error.output.strip('\n'))
    
    return(vcfFile)
    

# Assess the alignment - visualization (maybe in a separate script?)

    # Using bcftools view



parser = argparse.ArgumentParser(description='Script takes sam files and performs SNP calling against a reference sequence using bcftools mpileup')
parser.add_argument('-d', 'directory',type=str,help='Enter the path to the folder containing the sam files to be used')
parser.add_argument('reference',type=str,help='Enter the path to the reference isolate')
args = parser.parse_args()


samFiles = sorted(glob.glob(args.directory+'*.sam'))
logFile = 'snpcaller.log'

bwaIndexFile = bwaIndex(args.reference)
sortedReferenceFile = sortReference(args.reference)
seqDictFile = sequenceDictionary(args.reference,args.picardPath)



if __name__ == '__main__':
    # Run mpileup here

