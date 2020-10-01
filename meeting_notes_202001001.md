# Notes on variant calling and potential next steps

What has been done since last meeting:
1. Technical call with Lee Katz, Taylor Griswold, Sean Lucking, Jessica Rowell (September 29, 2020)
2. Started preliminary run of CFSAN SNP Pipeline on 266 isolates
3. Wrote a script that arranges files for CFSAN SNP Pipeline (github.com/jessicarowell/Helpers)  

Note-taker: Jessica

## Background

I scheduled this call to discuss the comparison of GATK and Samtools variant calling pipeline results.

I have never done variant calling prior to this, and my first foray into it was writing a pipeline.  Initially, Sean and I had conducted a literature review to find out what tools people use to do variant calling and whether there were existing pipelines to do this. Since I had never done variant calling, I focused on the tools used and came up with the plan to compare GATK and Samtools on a set of amplicon data that matches what we expect to feed the pipeline. Sean found some pipelines, but none of them appeared to do what we needed.  The results of our search are summarized in an Excel spreadsheet (path is in the OneNote notebook).

We both wrote pipelines (Sean: GATK, Jess: Samtools) and ran them on 266 salmonella isolates. These 266 isolates are fastq files that Sean made from fasta files (using made-up quality data). The fasta files comprise synthetic reads pulled out via in silico PCR using the T3PO primers (93% homology across genomes).

The main questions and discussion points for this call were:

1. How variant calling results validation works

2. Why do the very large files produce such different results for GATK vs. Samtools?

3. If you have two pipelines that both implement the same tools, what's the validation process for those pipelines? Don't you expect the same results from both?

## An alternate pipeline

Lee directed us to CFSAN SNP Pipeline, which looks like it might do what we need. They compared Lyve-SET to this pipeline in their paper, and the two papers on it are referenced there (Davis 2015 and Pettengill 2014).

Validation of a variant calling method takes an extensive amount of work, so we want to avoid re-inventing the wheel if possible.

I looked into CFSAN further, and I am running it on the 266 isolates now. Someone/something interrupted the program last night at 1:33 AM, so I'm now looking into either running it on Aspen (emailed SciComp about this) or running it on Monolith3 over the weekend. Monolith3 is experiencing someweird issues.

CFSAN SNP Pipeline has really nice documentation for installation and usage.  It requires 8 programs, all of which we have modules for, and all the latest versions work with the exception of GATK (v3.7 must be used instead of v4).

## Some thoughts about our initial results

Taylor said that some salmonella serotypes (e.g. ones that carry plasmids or ones with a larger range of SNPs) have a history of "behaving badly."  I asked her whether there is any documentation regarding serotypes and she said no, it's all "common knowledge among them" so we would just have to rely on noticing strange things, thinking of the right questions, and then asking them.  In my opinion, this is problematic.

She noticed that among the 8 outliers were 5 non-subspecies type I serotypes and said that there is a lot of variation among non-subspecies type I serotypes.  I asked about the 3 enteritiditis isolates, and she said those are a special case of subspecies type I that has a lot of mobility and can be highly clonal. So one thing we are certainly dealing with is differences in mapping by salmonella subspecies.


We should also compare the default options between GATK and Samtools.  Personally, I suspect this the cause of the 3x difference in number of variants called lies here.  Maybe Samtools lacks a filtering step that GATK implements at the beginning, so it's calling from all the bases rather than only the ones that mapped (for example).

In my opinion, we're dealing with two major variables here: 
1. GATK vs. Samtools performance
2. Salmonella serotype variation in mapping ability

## Potential next steps

We could run both pipelines again, using the original reference from which the synthetic reads were generated. Our null hypothesis in this case would be 0 SNPs called, so this could be a form of "sanity check."

We can run CFSAN SNP Pipeline on all 266 isolates and compare the results to ours. If CFSAN does what we want, we should just use it for our downstream analysis. I'm already doing a preliminary run and so far it seems very straightforward and easy to use, and contains a lot of informative logging as well.

I'd like to quickly make scatterplots for Samtools vs. GATK and color-code by subspecies I (enteritidis) vs. subspecies I (others) vs. non-subspecies I just to see if anything jumps out. 
