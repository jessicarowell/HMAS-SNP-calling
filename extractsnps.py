#!/usr/bin/env python
# coding: utf-8

import os
import argparse
import glob
from itertools import combinations
from functools import partial
import multiprocessing
import subprocess

class Pairwise:
    """Calculated pairwise differnce counts between two isolate multifasta files

    Attributes:
        isolate1: String representing the name of the first isolate in the comparison
        snpDiff: Integer representing the number of snp differnces calculated by GATK HaplotypeCaller

    """
    def __init__(self,isolate1,snpDiff):
        self.isolate1 = isolate1
        self.snpDiff = snpDiff 

    """
    Returns a new Pairwise object
    """

def vcfOutputParser(logList,baseName): 
    """Parses vcf output file

    Params
    ------
    logList: String
        A list of the lines in the log file created by vcf filtering command
    baseName: String
        Basename of the isolate 

    Returns
    -------
    snpCount: String
        The number of SNPs after filtering

    """   
    snps = 0
    #baseName = baseName.split('/')[-1]
    
    for line in logList:
        if line.startswith('After filtering'):
            if 'possible' in line:
                snps = int(line.split(' ')[3])

    snpCount = Pairwise(baseName,snps)
    return(snpCount)

def extractSNPs(vcfFile, baseName):
    """Performs vcf filtering

    Params
    ------
    vcfFile: String
        Input vcf file
 
    Returns
    -------
    snpCount: String
        The number of SNPs after filtering

    """

    r1 = subprocess.Popen(('vcftools','--vcf',vcfFile,'--freq'),stderr=subprocess.PIPE)
    r2 =(r1.stderr.read().decode('ascii'))
    r2 = r2.split('\n')
    snpCount = vcfOutputParser(r2,baseName)
    return(snpCount)

def listener_process(q):
    """Sets the logging file for the listener process

    Params
    ------
    q: queue
         queue for logging events

    Returns
    -------

    """
    with open(logFile,'w') as f:
        while 1:
            msg = q.get()
            if msg == 'kill':
                break
            f.write(str(msg)+'\n')
            f.flush()

def bcftoolsParallelFunctions(vcfFile, q):
    """Run SNP calling functions 

    Params
    ------
    vcfFile: String
        File path of the input vcf file 

    Returns
    -------
    snpCount: String
    	Count of SNP differences between isolate and reference

    """
    baseName = os.path.splitext(os.path.basename(vcfFile))[0]
    snpCount = extractSNPs(vcfFile, baseName)
    return(snpCount)

def ParallelRunner(parallelFunctions, inputFiles):
    """ Run functions on multiple files in parallel

    Params
    ------
    myFunction: String
        Name of the function to run
    inputFiles: String
        list of input files to run in parallel

    Returns
    -------
    
    """

    manager = multiprocessing.Manager()
    q = manager.Queue()
    pool = multiprocessing.Pool(processes=args.numcores)
    watcher = pool.apply_async(listener_process,(q,))
    parallelStatic = partial(parallelFunctions,q=q)
    result_list = pool.map(parallelStatic,inputFiles)
    q.put('kill')
    pool.close()
    pool.join()
    return(result_list)

parser = argparse.ArgumentParser(description='Script takes bam files and performs SNP calling against a reference sequence using bcftools mpileup')
parser.add_argument('-d', '--directory',type=str,help='Enter the path to the folder containing the bam files to be used')
parser.add_argument('-n', '--numcores',type=int,help='Enter the number of cores you would like to use during processing')
parser.add_argument('-o', '--output',type=str,help='Enter a name for the output file')
args = parser.parse_args()


# Starting from the vcf files
vcfFiles = sorted(glob.glob(args.directory+'/'+'*deDup.vcf'))
logFile = 'snpcaller.log'



# Use parallel processing from Sean's functions to process the samFiles
if __name__ == '__main__':
    # Run mpileup here and tabulate the SNP counts
    results = ParallelRunner(bcftoolsParallelFunctions, vcfFiles)


# Output isolate SNP count to a file
outFile = open(args.output, 'w')
for i in results:
    print(i.isolate1+'\t'+str(i.snpDiff), file=outFile)
outFile.close()


