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

def calcCoverage(refSeq, bamFile, baseName):
    """Uses bcftools mpileup to read coverage of positions in genome

    Params
    ------
    refSeq: String
        File location of the reference multifasta file (already indexed)
    bamFile: String
        Name of the input bam file
    baseName: String
        Basename of the isolate 

    Returns
        -------
        bcfFile: String
             File location of the output bcf file

    """
    inBam = os.path.dirname(bamFile)
    bcfFile = inBam+'/'+baseName+'.bcf'
    
    try:
        r1 = subprocess.check_output(['bcftools', 'mpileup', '-O','b','-o', bcfFile,'-f', refSeq, bamFile])
        # options: -O b to generate bcf format output, -f flags path to reference genome
    except subprocess.CalledProcessError as error:
        errorMsg = ('Status : Fail',error.returncode,error.output.strip('\n'))
    
    return(bcfFile)


def callVariants(bcfFile, baseName):
    """Uses bcftools call to detect SNPs

    Params
    ------
    bcfFile: String
        Input bcf file
    baseName: String
        Basename of the isolate 

    Returns
    -------
    vcfFile: String
        File location of the output vcf file

    """
    inBcf = os.path.dirname(bcfFile)
    vcfFile = inBcf+'/'+baseName+'.vcf'
    
    try:
        r1 = subprocess.check_output(['bcftools', 'call', '--ploidy','1', '-m','-v','-o', vcfFile, bcfFile])
        # options: -m for multiallelic & rare-variant calling, -v to output variant sites only, ploidy is 1 for haploid genome 
    except subprocess.CalledProcessError as error:
        errorMsg = ('Status : Fail',error.returncode,error.output.strip('\n'))
    
    return(vcfFile)
    


def vcfOutputParser(logList, baseName):
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


def bcftoolsParallelFunctions(bamFile, q):
    """Run SNP calling functions 

    Params
    ------
    bamFile: String
        File path of the input bam file 

    Returns
    -------
    snpCount: String
        Count of SNP differences between isolate and reference

    """
    baseName = os.path.basename(bamFile)
    bcf = calcCoverage(args.reference, bamFile, baseName)
    vcf = callVariants(bcf, baseName)
    snpCount = extractSNPs(vcf, baseName)
    
    return(snpCount)



def ParallelRunner(parallelFunctions, inputFiles):
    """ Run functions on multiple files in parallel

    Params
    ------
    parallelFunctions: String
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
parser.add_argument('-r', '--reference',type=str,help='Enter the path to the reference isolate')
parser.add_argument('-n', '--numcores',type=int,help='Enter the number of cores you would like to use during processing')
args = parser.parse_args()


# Starting from the de-duplicated bam files
bamFiles = sorted(glob.glob(args.directory+'/'+'*deDup.bam'))
logFile = 'snpcaller.log'



# Use parallel processing from Sean's functions to process the samFiles
if __name__ == '__main__':
    # Run mpileup here
    results = ParallelRunner(bcftoolsParallelFunctions, bamFiles)


# Output isolate SNP count to a file (once I have written a function to tabulate the SNPs called)
