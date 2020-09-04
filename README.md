# Summary of current directory

Contains info on the SNP Calling pipeline, for downstream of analysis after QC of HMAS data

## Location of Data

Current Path: `/scicomp/home/ick4/data/snpcalling`  (referred to herein as `$HOME/$PROJECT`)
Path for data used in analysis: `$HOME/$PROJECT`

Raw data located here: `/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/Salmonella/SnpCalling`



## Description of Contents

`$HOME/$PROJECT/snpcaller.py`
Script for calling SNPs from HMAS data.  
At the moment this script follows from Sean's SNP calling (GATK) script and takes in sorted, de-duplicated .sam files

`$HOME/$PROJECT/SNP_Calling_Flowchart.pdf`
A simple flowchart to track the logic of the SNP calling
	
Created by: Jessica Rowell
Last modified: 8/10/2020

		
## Data manipulations done in SNP Calling tests

After running snpcaller.py on the 266 files in DesignSet and LyveSet, some of them stalled.
I created a table to show bam, bcf, and vcf files (columns) for each sample (row)
There are blanks where files did not get generated.

4 of the 5 samples that didn't finish were excessively large files (>1G) and 1 is corrupt.

### Trying to figure out which files mpileup could not complete

```
import os, glob
import numpy as np
import pandas as pd

def check_if_exists(x, orig_list, new_list):
	na = ""
	if any(x in s for s in orig_list):
		orig_x = [i for i in orig_list if x in i]
		return(orig_x)
	else:
		return(na)


bams = sorted(glob.glob('./*deDup.bam'))
bcfs = sorted(glob.glob('./*deDup.bcf'))
vcfs = sorted(glob.glob('./*deDup.vcf'))

bams = [os.path.basename(x) for x in bams]
bcfs = [os.path.basename(x) for x in bcfs]
vcfs = [os.path.basename(x) for x in vcfs]

basenames = [os.path.splitext(os.path.basename(x))[0] for x in bams]

base_bams = [os.path.splitext(os.path.basename(x))[0] for x in bams] # Same as basenames
base_bcfs = [os.path.splitext(os.path.basename(x))[0] for x in bcfs]
base_vcfs = [os.path.splitext(os.path.basename(x))[0] for x in vcfs]



df = pd.DataFrame(columns = np.arange(3))
for f in basenames:
	bam = check_if_exists(f, bams, new_list)
	bcf = check_if_exists(f, bcfs, new_list)
	vcf = check_if_exists(f, vcfs, new_list)
	pd_list = pd.Series([bam,bcf,vcf], index = df.columns)
	df = df.append(pd_list, ignore_index = True)

df.to_csv('snp_caller_struggles.txt',sep='\t', header = False, index = False, index_label = None)
```

Struggled to figure out how to do this with pandas, so after 5 minutes of Googling I just used sed

```
sed -i "s/'//g"  snp_caller_struggles.txt
sed -i "s/\[//g"  snp_caller_struggles.txt
sed -i "s/\[//g"  snp_caller_struggles.txt
```

Files that samtools mpileup choked on:
1.8G	Sal_JEG_EC20121176deDup.bam # ran this one by hand; took a full day
1.2G	Sal_JEG_EC20121179deDup.bam
1.3G	Sal_JEG_EC20121180deDup.bam
856K	Sal_JEG_OLF-SE4-0317-8deDup.bam
3.4M	Sal_JR3_62-z36_RKS2983deDup.bam # bcf was produced but not the vcf (tried again and it was still corrupted 'No BGZF EOR marker')

Commands to run by hand:
bcftools mpileup -O b -o Sal_JEG_EC20121176deDup.bcf -f /scicomp/home/ick4/data/snpcalling/test/ReferenceSequenceFiles/GCA_000022165.1_ASM2216v1_genomic.fna Sal_JEG_EC20121176deDup.bam

bcftools call --ploidy 1 -m -v -o Sal_JEG_EC20121176deDup.vcf Sal_JEG_EC20121176deDup.bcf

## END 
