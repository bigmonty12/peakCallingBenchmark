# Benchmarking ATAC-seq peak calling
ATAC-seq (Assay for Transposase-Accessible Chromatin using Sequencing) is becoming increasing popular as a method to investigate the accessibility of DNA in different biological states. Our lab has been using ATAC-seq for the last couple of years to see how Drosophila chromatin changes in response to alcohol exposure (Yes, we get flies drunk! It's pretty interesting!). As our lab's bioinformatician, I have developed an in-house pipeline for our analysis needs. The pipeline is based off suggestions from ENCODE, Harvard FAS Informatics, and other scientific bodies/individuals. While most scientists agree about the majority of the steps needed to analyze ATAC-seq data, there is a huge question mark about how to identify accessible regions, also known as 'peaks.' Until very recently, the community simply tried to adapt tools originally intended for CHiP-seq data (MACS2, f-seq, ZINBA). This past year brought HMMRATAC, an ATAC-seq specific peak caller, as well as the unpublished Genrich, a new general peak-caller with an ATAC-seq option. The goal of this project is to benchmark the performance of common ATAC-seq peak-calling methods.
## Peak Callers
Based on published papers and message boards, MACS2 seems to be the most popular tool to analyze ATAC-seq peaks. Though unpublished, Genrich appears to be gaining traction among users. HMMRATAC is also increasingly mentioned, often with caveats about run time. This section is dedicated to these three methods (as well as multiple MACS2 options).
### MACS2
MACS2, or Model-based Analysis for ChIP-Seq, is the most well-known and -used tool for peak calling. As evidenced by its name, it was originally developed for ChIP-seq. It uses the Poisson distribution as the null basis for detecting genome biases and enrichment. A sliding window technique is employed to find more accessible regions. It was originally created by [Tao Liu]().
One problem with MACS2 is that there is no consensus on how to use it with ATAC-seq. There are three commonly used options. First, the default is to use MACS2 in 'BAM' mode and extend each read fragment in both directions. This mode throws away half of properly paired data and has issues with false positive peaks due to fragments being shifted and extended. The second common method is to use 'BAMPE' mode. This method utilizes both reads of paired-end data and infers full fragments instead of cut sites. This method ignores any extension or read shift parameters. A third manner in which to use MACS2 for ATAC-seq is to use 'BED' mode. In this case, a paired-end BAM file is converted to a BED file, and then the extension and shift parameters can be used. It is, in a sense, the combination of the typical 'BAM' and 'BAMPE' methods.
### Genrich
Genrich was recently developed by [John Gaspar](), previous of Harvard FAS Informatics. The default method is for ChIP-seq peak-calling but there is an ATAC-seq specific mode which can be enabled by the `-j` parameter. To account for the manner in which the Tn5 enzyme inserts itself into the genome, Genrich centers a 100 bp interval on the end of each correlating fragment. It uses the log-normal distribution to calculate p-values for each base of the genome.
Some advantages over MACS include how much faster it runs, as well as the capability of calling peaks for multiple replicates collectively. Unfortunately, Genrich is yet to be published, creating some hesitancy among researchers.
### HMMRATAC
HMMRATAC is the only peak-caller designed specifically for ATAC-seq. It was developed by Eddie Tarbell and Tao Liu (of MACS2 fame). HMMRATAC uses a three-state Hidden Markov Model (HMM) to segment the genome into open chromatin regions, nucleosomal regions, and background regions. It employs a semi-supervised machine learning approach to train its HMM states.
Though HMMRATAC has much potential as a peak-caller, it is quite computationally intensive. It requires a lot of memory and often takes hours to run. HMMRATAC also creates a gappedPeak file rather than narrowPeak files like the tools above. While not inherently wrong, it does not easily lend itself to traditional peak understanding or differential peak finding.
## Materials and Methods
### Materials
The human GM12878 cell line ATAC-seq paired-end data used was publicly available and downloaded under five SRA accession numbers. The first three SRA accession numbers, SRR891269-SRR891271, were used in the HMMRATAC paper [] while the other two SRA samples, SRR5427884-SRR5427885, were used by Gaspar in his preprint, Improving Peak-Calling with MACS2.
### Preprocessing
I created a [Nextflow script]() to automate all of the preprocessing steps. To summarize, fastq pairs are first evaluated with [FastQC](). [TrimGalore]() is used to remove adapters and low quality reads from the fastq pairs. The reads are then aligned against the hg38 genome with [Bowtie2]() and sorted with [Samtools](). [Picard]() is used to remove duplicates. Note that for [Genrich](), the bam file must be sorted by read names which is also accomplished with [Samtools]() but with the parameter `-n`. Also, replicates were merged with [Samtools]() and then preprocessed using the aforementioned steps.
### 'Gold Standard' Lists
The 'gold standard' dataset was downloaded from the supplementary data from the HMMRATAC paper, which was originally downloaded from the UCSC genome browser track wgEncodeBroadHmmGm12878HMM. This dataset, called 'active regions', is made by merging two chromatin states 'Active Promoters' and 'Strong Enhancers' generated by [chromHMM]() in GM12878 cells.
The 'real negative' dataset was also downloaded from the supplementary data from the HMMRATAC paper. This set came from the 'Heterochromatin' state annotation from [chromHMM](). As both BED files were originally created with the hg19 genome, I used the USCS [LiftOver]() tool to convert the coordinates to hg38.
## Results

## Conclusion

### Ideas
Find summaries of per position reads
```
bedtools coverage -a A.narrowPeak -b A.bam -d > A.per-base-coverage.txt
awk '{if ($12 <= 1) { print $1, $2, $3, $11, $12 } }' A.per-base-coverage.txt > A.false-positives.txt
```

Get total base pairs in bed file (make sure to `bedtools merge` first so base pairs aren't counted more than one time.
```
cat $bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
```

Find overlapping base pairs (between called peak and real negative/positive) (make sure to use `bedtools merge` on both files so base pairs aren't counted more than once)
```
bedtools intersect -a standard.bed -b peaks.bed -wo | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$(NF) }END{print SUM}'
```

