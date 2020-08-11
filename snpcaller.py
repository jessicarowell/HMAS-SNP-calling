import os
import argparse
import glob
from Bio import SeqIO
from itertools import combinations
from functools import partial
import multiprocessing
import subprocess


# Calculate read coverage of positions in the genome





# Detect the SNPs
def callVariants(refSeq,fastq,baseName):

    """Uses bcftools mpileup to call SNVs

    Params
    ------
    refSeq: String
         File location of the reference multifasta file
    bamFile: String
         File location of the input bam file
    baseName: String
         Basename of the isolate being aligned against the reference

    Returns
    -------
    samFile: String
         File location of the created sam file

    """

    samFile = baseName+'.sam'
    
    try:
        mpileup_cmd = 'bwa mem -R "@RG\\tID:FLOWCELL1.LANE1\\tPL:ILLUMINA\\tLB:test\\tSM:PA01" '+refSeq+' '+fastq+' > '+samFile
        if os.system(bwa_cmd) != 0:
            raise Exception('Aint work yo')
    except:
        print('borked')

    return(samFile)

def callVariants(refSeq,deDupBam,baseName):


    try:
        r1 = subprocess.check_output(['gatk','--java-options','-Xmx4g','HaplotypeCaller','-ploidy','1','-R',refSeq,'-I',deDupBam,'-O',baseName+'.vcf'])
    except subprocess.CalledProcessError as error:
        errorMsg = ('Status : Fail',error.returncode,error.output.strip('\n'))
    return(baseName+'.vcf')


# Filter and report the SNP variants in variant calling format









# Assess the alignment - visualization (maybe in a separate script?)


parser = argparse.ArgumentParser(description='Script takes sam files and performs SNP calling against a reference sequence using bcftools mpileup')
parser.add_argument('-d', 'directory',type=str,help='Enter the path to the folder containing the sam files to be used')
parser.add_argument('reference',type=str,help='Enter the path to the reference isolate')
args = parser.parse_args()


samFiles = sorted(glob.glob(args.directory+'*.sam'))
logFile = 'snpcaller.log'

bwaIndexFile = bwaIndex(args.reference)
sortedReferenceFile = sortReference(args.reference)
seqDictFile = sequenceDictionary(args.reference,args.picardPath)

# use refFileCheck from Sean's script to check the ref file?
refPass = refFileCheck(bwaIndexFile,sortedReferenceFile,seqDictFile)

if refPass == False:
    os.exit('Problem encountered with the reference file')


if __name__ == '__main__':
    # Run mpileup here

