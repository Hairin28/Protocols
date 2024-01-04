- [RNA-seq Data Processing](#rna-seq-data-processing)
  - [1. Download and process SRA files](#1-download-and-process-sra-files)
    - [1.1 Download sra files below](#11-download-sra-files-below)
    - [1.2 Process using sratoolkit](#12-process-using-sratoolkit)
  - [2. fastqc](#2-fastqc)
  - [3. Cut adaptor](#3-cut-adaptor)
  - [4. Trim to 40bp](#4-trim-to-40bp)
  - [5. Map reads to reference genome](#5-map-reads-to-reference-genome)
    - [5.1 BWA mapping](#51-bwa-mapping)
  - [6. Remove blacklist](#6-remove-blacklist)
  - [7. Convert to bigwig](#7-convert-to-bigwig)
  - [8. Featurecounts](#8-featurecounts)
  - [9. further quantified analysis](#9-further-quantified-analysis)
- [ChIP-seq Data Processing](#chip-seq-data-processing)
  - [1. download file - chip SRR4892826 SRR4892828 \& input SRR4892865 SRR4892866](#1-download-file---chip-srr4892826-srr4892828--input-srr4892865-srr4892866)
  - [2. fastqc](#2-fastqc-1)
  - [2. cut adaptor](#2-cut-adaptor)
  - [3. fastqc](#3-fastqc)
  - [4. trim to 30bp](#4-trim-to-30bp)
  - [4. bwa mapping](#4-bwa-mapping)
  - [5. samtobam, sorted \& index](#5-samtobam-sorted--index)
  - [6. remove duplicates](#6-remove-duplicates)
  - [7. remove blacklist](#7-remove-blacklist)
  - [8. convert to bigwig](#8-convert-to-bigwig)
  - [9. merge replicates](#9-merge-replicates)
  - [10. convert merged files to bigwig](#10-convert-merged-files-to-bigwig)
  - [11. peakcall, RiP calculation and so on](#11-peakcall-rip-calculation-and-so-on)


# RNA-seq Data Processing

- Use "RNA-seq data of ezh2 cKO male PGC (Huang et al.,2021)" as an example

## 1. Download and process SRA files

### 1.1 Download sra files below

|method| |stage| | | |url|
|:----|:----|:----|:----|:----|:----|:----|
|RNA-seq|Ezh2 CKO male PGC rep1|E13.5|male|ezh2_cKO|GSM4196632|https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR10560079|
|RNA-seq|Ezh2 CKO male PGC rep2|E13.5|male|ezh2_cKO|GSM4196633|https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR10560081|
|RNA-seq|Ezh2 CKO male PGC rep3|E13.5|male|ezh2_cKO|GSM4196634|https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR10560083|
|RNA-seq|Ezh2 Ctrl male PGC rep1|E13.5|male|Ezh2_ctrl|GSM4196638|https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR10560091|
|RNA-seq|Ezh2 Ctrl male PGC rep2|E13.5|male|Ezh2_ctrl|GSM4196639|https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR10560093|
|RNA-seq|Ezh2 Ctrl male PGC rep3|E13.5|male|Ezh2_ctrl|GSM4196640|https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR10560095|

### 1.2 Process using sratoolkit

```bash
# take one file SRR10560079 as example
$ module load sratoolkit
$ prefetch SRR10560079
$ fastq-dump SRR10560079 --split-files # --split-files for paired-end files
```

## 2. fastqc

```bash
$ module load FastQC
$ fastqc SRR10560079_1.fastq  SRR10560079_2.fastq -o ./ 
# fastqc results check
# check Per base sequence quality > 28
# check Sequence Length Distribution and Adapter Content
```

## 3. Cut adaptor

```bash
$ module load cutadapt 
# AGATCGGAAGAG is an Illumina Universe adaptor, change it for specific samples

# paired-end mode for SRR10560079 file
$ cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -m 20 -o SRR10560079_1_cutadapt.fastq -p SRR10560079_2_cutadapt.fastq SRR10560079_1.fastq SRR10560079_2.fastq -j 8

# single-end mode for SRR10560079 file
$ cutadapt -b AGATCGGAAGAG -m 20 -j 8 -o XXXX_cutadapt.fastq XXXX.fastq 
```

## 4. Trim to 40bp

```bash
# specifically, SRR10560079 adaptors start from 40bp position of read, trim the reads to 40bp
$ module load seqkit
$ cat SRR10560079_1_cutadapt.fastq | seqkit subseq -r 1:40 > SRR10560079_1_cutadapt_40bp.fastq
$ cat SRR10560079_2_cutadapt.fastq | seqkit subseq -r 1:40 > SRR10560079_2_cutadapt_40bp.fastq
```

## 5. Map reads to reference genome  

- several tools for mapping, hisat2, bowtie2, BWA et al.
  
### 5.1 BWA mapping

```bash
# first download and index mm10 genome fa file, using it as reference in the following steps
$ wget "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz"
$ gunzip mm10.fa.gz

$ module load bwa
$ bwa index mm10.fa -p ./mm10

# for paired-end reads, input read1.fastq read2.fastq together
$ bwa mem -t 8 -M /home/Li/Log/test_data/mm10 SRR10560079_1_cutadapt_40bp.fastq SRR10560079_2_cutadapt_40bp.fastq > SRR10560079_cutadapt_40bp.sam 2> SRR10560079_cutadapt_40bp.log
# for single-end reads
$ bwa mem -t 8 -M /home/Li/Log/test_data/mm10 SRR10560079_cutadapt_40bp.fastq > SRR10560079_cutadapt_40bp.sam 2> SRR10560079_cutadapt_40bp.log
```

## 6. Remove blacklist

```bash
# download blacklist file from https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz
# move it to the working filter and gunzip it with $ gunzip mm10-blacklist.v2.bed.gz

$ module load samtools
$ samtools view -Sb SRR10560079_cutadapt_40bp.sam | samtools sort > SRR10560079_cutadapt_40bp_sorted.bam
$ bedtools intersect -v -abam SRR10560079_cutadapt_40bp_sorted.bam -b mm10-blacklist.v2.bed > SRR10560079_cutadapt_40bp_sorted_rmblack.bam
```

## 7. Convert to bigwig

```bash
$ samtools index SRR10560079_cutadapt_40bp_sorted_rmblack.bam
$ bamCoverage --bam SRR10560079_cutadapt_40bp_sorted_rmblack.bam -o SRR10560079_cutadapt_40bp_sorted_rmblack.bw -of bigwig --binSize=10 --minMappingQuality 10
```

## 8. Featurecounts

```bash
# download gtf file from https://nov2020.archive.ensembl.org/Mus_musculus/Info/Index
$ featureCounts -p -t exon -g gene_id -a Mus_musculus.GRCm38.102.gtf -o SRR10560079_featureCounts.txt SRR10560079_cutadapt_40bp_sorted_rmblack.bam # -t exon or gene
```

## 9. further quantified analysis

```
e.g. PCA, DESeq2 log2FC, GO term analysis, TPM & FPKM calculation and so on
```


<br/>


***
# ChIP-seq Data Processing

- Use SSC H3K27me3 ChIP-seq dataset GSE89502 as an example

## 1. download file - chip SRR4892826 SRR4892828 & input SRR4892865 SRR4892866

```bash
$ prefetch SRR4892826 SRR4892828|fastq-dump SRR4892826 SRR4892828
$ prefetch SRR4892865 SRR4892866|fastq-dump SRR4892865 SRR4892866
```

## 2. fastqc

```bash
$ module load FastQC
$ fastqc /home/Li/DataBase/SSC/fastq/SRR4892826.fastq -o /home/Li/DataBase/SSC/fastqc
$ fastqc /home/Li/DataBase/SSC/fastq/SRR4892828.fastq -o /home/Li/DataBase/SSC/fastqc
$ fastqc /home/Li/DataBase/SSC/fastq/SRR4892865.fastq -o /home/Li/DataBase/SSC/fastqc
$ fastqc /home/Li/DataBase/SSC/fastq/SRR4892866.fastq -o /home/Li/DataBase/SSC/fastqc
```

## 2. cut adaptor

```bash
# cut two adaptors CTGTCTCTTATACACATCT & AGATCGGAAGAG, single-end mode
$ module load cutadapt 
$ cutadapt -b CTGTCTCTTATACACATCT -b AGATCGGAAGAG -m 20 -j 8 -o /home/Li/DataBase/SSC/cutadapt/SRR4892826_cutadapt.fastq /home/Li/DataBase/SSC/fastq/SRR4892826.fastq
$ cutadapt -b CTGTCTCTTATACACATCT -b AGATCGGAAGAG -m 20 -j 8 -o /home/Li/DataBase/SSC/cutadapt/SRR4892828_cutadapt.fastq /home/Li/DataBase/SSC/fastq/SRR4892828.fastq
$ cutadapt -b CTGTCTCTTATACACATCT -b AGATCGGAAGAG -m 20 -j 8 -o /home/Li/DataBase/SSC/cutadapt/SRR4892865_cutadapt.fastq /home/Li/DataBase/SSC/fastq/SRR4892865.fastq
$ cutadapt -b CTGTCTCTTATACACATCT -b AGATCGGAAGAG -m 20 -j 8 -o /home/Li/DataBase/SSC/cutadapt/SRR4892866_cutadapt.fastq /home/Li/DataBase/SSC/fastq/SRR4892866.fastq
```

## 3. fastqc

```bash
$ module load FastQC
$ fastqc /home/Li/DataBase/SSC/cutadapt/SRR4892826_cutadapt.fastq -o /home/Li/DataBase/SSC/fastqc
$ fastqc /home/Li/DataBase/SSC/cutadapt/SRR4892828_cutadapt.fastq -o /home/Li/DataBase/SSC/fastqc
$ fastqc /home/Li/DataBase/SSC/cutadapt/SRR4892865_cutadapt.fastq -o /home/Li/DataBase/SSC/fastqc
$ fastqc /home/Li/DataBase/SSC/cutadapt/SRR4892866_cutadapt.fastq -o /home/Li/DataBase/SSC/fastqc
# fastqc results check
# check Per base sequence quality > 28
# check Sequence Length Distribution and Adapter Content
```

## 4. trim to 30bp

```bash
$ module load seqkit
$ cat /home/Li/DataBase/SSC/cutadapt/SRR4892826_cutadapt.fastq | seqkit subseq -r 1:30 > /home/Li/DataBase/SSC/trimming/SRR4892826_cutadapt_30bp.fastq
$ cat /home/Li/DataBase/SSC/cutadapt/SRR4892828_cutadapt.fastq | seqkit subseq -r 1:30 > /home/Li/DataBase/SSC/trimming/SRR4892828_cutadapt_30bp.fastq
$ cat /home/Li/DataBase/SSC/cutadapt/SRR4892865_cutadapt.fastq | seqkit subseq -r 1:30 > /home/Li/DataBase/SSC/trimming/SRR4892865_cutadapt_30bp.fastq
$ cat /home/Li/DataBase/SSC/cutadapt/SRR4892866_cutadapt.fastq | seqkit subseq -r 1:30 > /home/Li/DataBase/SSC/trimming/SRR4892866_cutadapt_30bp.fastq
```

## 4. bwa mapping

```bash
# use the indexed mm10 fa files, which are in /home/Li/DataBase/BWA/mm10/mm10
$ module load bwa
$ bwa mem -t 8 -M /home/Li/DataBase/BWA/mm10/mm10 /home/Li/DataBase/SSC/trimming/SRR4892826_cutadapt_30bp.fastq > /home/Li/DataBase/SSC/bwa/SRR4892826_cutadapt_30bp.sam 2> /home/Li/DataBase/SSC/bwa/SRR4892826_cutadapt_30bp.log
$ bwa mem -t 8 -M /home/Li/DataBase/BWA/mm10/mm10 /home/Li/DataBase/SSC/trimming/SRR4892828_cutadapt_30bp.fastq > /home/Li/DataBase/SSC/bwa/SRR4892828_cutadapt_30bp.sam 2> /home/Li/DataBase/SSC/bwa/SRR4892828_cutadapt_30bp.log
$ bwa mem -t 8 -M /home/Li/DataBase/BWA/mm10/mm10 /home/Li/DataBase/SSC/trimming/SRR4892865_cutadapt_30bp.fastq > /home/Li/DataBase/SSC/bwa/SRR4892865_cutadapt_30bp.sam 2> /home/Li/DataBase/mESC/bwa/SRR4892865_cutadapt_30bp.log
$ bwa mem -t 8 -M /home/Li/DataBase/BWA/mm10/mm10 /home/Li/DataBase/SSC/trimming/SRR4892866_cutadapt_30bp.fastq > /home/Li/DataBase/SSC/bwa/SRR4892866_cutadapt_30bp.sam 2> /home/Li/DataBase/mESC/bwa/SRR4892866_cutadapt_30bp.log
```

## 5. samtobam, sorted & index

```bash
$ module load samtools

$ samtools view -@ 8 -Sb /home/Li/DataBase/SSC/bwa/SRR4892826_cutadapt_30bp.sam | samtools sort -@ 8 > /home/Li/DataBase/SSC/bwa/SRR4892826_cutadapt_30bp_sorted.bam
$ samtools index -@ 8 /home/Li/DataBase/SSC/bwa/SRR4892826_cutadapt_30bp_sorted.bam

$ samtools view -@ 8 -Sb /home/Li/DataBase/SSC/bwa/SRR4892828_cutadapt_30bp.sam | samtools sort -@ 8 > /home/Li/DataBase/SSC/bwa/SRR4892828_cutadapt_30bp_sorted.bam
$ samtools index -@ 8 /home/Li/DataBase/SSC/bwa/SRR4892828_cutadapt_30bp_sorted.bam

$ samtools view -@ 8 -Sb /home/Li/DataBase/SSC/bwa/SRR4892865_cutadapt_30bp.sam | samtools sort -@ 8 > /home/Li/DataBase/SSC/bwa/SRR4892865_cutadapt_30bp_sorted.bam
$ samtools index -@ 8 /home/Li/DataBase/SSC/bwa/SRR4892865_cutadapt_30bp_sorted.bam

$ samtools view -@ 8 -Sb /home/Li/DataBase/SSC/bwa/SRR4892866_cutadapt_30bp.sam | samtools sort -@ 8 > /home/Li/DataBase/SSC/bwa/SRR4892866_cutadapt_30bp_sorted.bam
$ samtools index -@ 8 /home/Li/DataBase/SSC/bwa/SRR4892866_cutadapt_30bp_sorted.bam
```

## 6. remove duplicates

```bash
$ module load picard

$ picard MarkDuplicates ASSUME_SORTED=true REMOVE_DUPLICATES=true \
            I=/home/Li/DataBase/SSC/bwa/SRR4892826_cutadapt_30bp_sorted.bam \
            O=/home/Li/DataBase/SSC/picard/SRR4892826_cutadapt_30bp_sorted_picard.bam \
            M=/home/Li/DataBase/SSC/picard/SRR4892826_cutadapt_30bp_sorted_picard.log
$ samtools index -@ 8 /home/Li/DataBase/SSC/picard/SRR4892826_cutadapt_30bp_sorted_picard.bam

$ picard MarkDuplicates ASSUME_SORTED=true REMOVE_DUPLICATES=true \
            I=/home/Li/DataBase/SSC/bwa/SRR4892828_cutadapt_30bp_sorted.bam \
            O=/home/Li/DataBase/SSC/picard/SRR4892828_cutadapt_30bp_sorted_picard.bam \
            M=/home/Li/DataBase/SSC/picard/SRR4892828_cutadapt_30bp_sorted_picard.log
$ samtools index -@ 8 /home/Li/DataBase/SSC/picard/SRR4892828_cutadapt_30bp_sorted_picard.bam

$ picard MarkDuplicates ASSUME_SORTED=true REMOVE_DUPLICATES=true \
            I=/home/Li/DataBase/SSC/bwa/SRR4892865_cutadapt_30bp_sorted.bam \
            O=/home/Li/DataBase/SSC/picard/SRR4892865_cutadapt_30bp_sorted_picard.bam \
            M=/home/Li/DataBase/SSC/picard/SRR4892865_cutadapt_30bp_sorted_picard.log
$ samtools index -@ 8 /home/Li/DataBase/SSC/picard/SRR4892865_cutadapt_30bp_sorted_picard.bam

$ picard MarkDuplicates ASSUME_SORTED=true REMOVE_DUPLICATES=true \
            I=/home/Li/DataBase/SSC/bwa/SRR4892866_cutadapt_30bp_sorted.bam \
            O=/home/Li/DataBase/SSC/picard/SRR4892866_cutadapt_30bp_sorted_picard.bam \
            M=/home/Li/DataBase/SSC/picard/SRR4892866_cutadapt_30bp_sorted_picard.log
$ samtools index -@ 8 /home/Li/DataBase/SSC/picard/SRR4892866_cutadapt_30bp_sorted_picard.bam
```

## 7. remove blacklist

```bash
$ bedtools intersect -v \
        -abam /home/Li/DataBase/SSC/picard/SRR4892826_cutadapt_30bp_sorted_picard.bam \
        -b /home/Li/DataBase/Mus_musculus/UCSC/mm10/mm10-blacklist.v2.bed \
        > /home/Li/DataBase/SSC/picard/SRR4892826_cutadapt_30bp_sorted_picard_rmblack.bam
$ samtools index -@ 8 /home/Li/DataBase/SSC/picard/SRR4892826_cutadapt_30bp_sorted_picard_rmblack.bam

$ bedtools intersect -v \
        -abam /home/Li/DataBase/SSC/picard/SRR4892828_cutadapt_30bp_sorted_picard.bam \
        -b /home/Li/DataBase/Mus_musculus/UCSC/mm10/mm10-blacklist.v2.bed \
        > /home/Li/DataBase/SSC/picard/SRR4892828_cutadapt_30bp_sorted_picard_rmblack.bam
$ samtools index -@ 8 /home/Li/DataBase/SSC/picard/SRR4892828_cutadapt_30bp_sorted_picard_rmblack.bam

$ bedtools intersect -v \
        -abam /home/Li/DataBase/SSC/picard/SRR4892865_cutadapt_30bp_sorted_picard.bam \
        -b /home/Li/DataBase/Mus_musculus/UCSC/mm10/mm10-blacklist.v2.bed \
        > /home/Li/DataBase/SSC/picard/SRR4892865_cutadapt_30bp_sorted_picard_rmblack.bam
$ samtools index -@ 8 /home/Li/DataBase/SSC/picard/SRR4892865_cutadapt_30bp_sorted_picard_rmblack.bam

$ bedtools intersect -v \
        -abam /home/Li/DataBase/SSC/picard/SRR4892866_cutadapt_30bp_sorted_picard.bam \
        -b /home/Li/DataBase/Mus_musculus/UCSC/mm10/mm10-blacklist.v2.bed \
        > /home/Li/DataBase/SSC/picard/SRR4892866_cutadapt_30bp_sorted_picard_rmblack.bam
$ samtools index -@ 8 /home/Li/DataBase/SSC/picard/SRR4892866_cutadapt_30bp_sorted_picard_rmblack.bam
```

## 8. convert to bigwig

```bash
$ bamCoverage -b /home/Li/DataBase/SSC/rmblack_bam/replicates/SRR4892826_cutadapt_30bp_sorted_picard_rmblack.bam -o /home/Li/DataBase/SSC/bigwig/SRR4892826_cutadapt_30bp_sorted_picard_rmblack.bw --binSize 50 --extendReads 50 -p 4

$ bamCoverage -b /home/Li/DataBase/SSC/rmblack_bam/replicates/SRR4892828_cutadapt_30bp_sorted_picard_rmblack.bam -o /home/Li/DataBase/SSC/bigwig/SRR4892828_cutadapt_30bp_sorted_picard_rmblack.bw --binSize 50 --extendReads 50 -p 4

$ bamCoverage -b /home/Li/DataBase/SSC/rmblack_bam/replicates/SRR4892865_cutadapt_30bp_sorted_picard_rmblack.bam -o /home/Li/DataBase/SSC/bigwig/SRR4892865_cutadapt_30bp_sorted_picard_rmblack.bw --binSize 50 --extendReads 50 -p 4

$ bamCoverage -b /home/Li/DataBase/SSC/rmblack_bam/replicates/SRR4892866_cutadapt_30bp_sorted_picard_rmblack.bam -o /home/Li/DataBase/SSC/bigwig/SRR4892866_cutadapt_30bp_sorted_picard_rmblack.bw --binSize 50 --extendReads 50 -p 4

# CPM normalization, not necessary for normal procession
$ bamCoverage -b /home/Li/DataBase/SSC/rmblack_bam/replicates/SRR4892826_cutadapt_30bp_sorted_picard_rmblack.bam -o /home/Li/DataBase/SSC/bigwig/SRR4892826_cutadapt_30bp_sorted_picard_rmblack.CPM.bw --binSize 50 --extendReads 50 -p 4 --normalizeUsing CPM

$ bamCoverage -b /home/Li/DataBase/SSC/rmblack_bam/replicates/SRR4892828_cutadapt_30bp_sorted_picard_rmblack.bam -o /home/Li/DataBase/SSC/bigwig/SRR4892828_cutadapt_30bp_sorted_picard_rmblack.CPM.bw --binSize 50 --extendReads 50 -p 4 --normalizeUsing CPM

$ bamCoverage -b /home/Li/DataBase/SSC/rmblack_bam/replicates/SRR4892865_cutadapt_30bp_sorted_picard_rmblack.bam -o /home/Li/DataBase/SSC/bigwig/SRR4892865_cutadapt_30bp_sorted_picard_rmblack.CPM.bw --binSize 50 --extendReads 50 -p 4 --normalizeUsing CPM

$ bamCoverage -b /home/Li/DataBase/SSC/rmblack_bam/replicates/SRR4892866_cutadapt_30bp_sorted_picard_rmblack.bam -o /home/Li/DataBase/SSC/bigwig/SRR4892866_cutadapt_30bp_sorted_picard_rmblack.CPM.bw --binSize 50 --extendReads 50 -p 4 --normalizeUsing CPM
```

## 9. merge replicates

```bash
$ samtools merge /home/Li/DataBase/SSC/rmblack_bam/merged_files/SSC_K27me3_merged.bam /home/Li/DataBase/SSC/rmblack_bam/replicates/SRR4892826_cutadapt_30bp_sorted_picard_rmblack.bam /home/Li/DataBase/SSC/rmblack_bam/replicates/SRR4892828_cutadapt_30bp_sorted_picard_rmblack.bam
$ samtools index /home/Li/DataBase/SSC/rmblack_bam/merged_files/SSC_K27me3_merged.bam

$ samtools merge /home/Li/DataBase/SSC/rmblack_bam/merged_files/SSC_WCE_merged.bam /home/Li/DataBase/SSC/rmblack_bam/replicates/SRR4892865_cutadapt_30bp_sorted_picard_rmblack.bam /home/Li/DataBase/SSC/rmblack_bam/replicates/SRR4892866_cutadapt_30bp_sorted_picard_rmblack.bam
$ samtools index /home/Li/DataBase/SSC/rmblack_bam/merged_files/SSC_WCE_merged.bam
```

## 10. convert merged files to bigwig

```bash
$ bamCoverage -b /home/Li/DataBase/SSC/rmblack_bam/merged_files/SSC_K27me3_merged.bam -o /home/Li/DataBase/SSC/bigwig/SSC_K27me3_merged.bw --binSize 50 --extendReads 50 -p 4

$ bamCoverage -b /home/Li/DataBase/SSC/rmblack_bam/merged_files/SSC_WCE_merged.bam -o /home/Li/DataBase/SSC/bigwig/SSC_WCE_merged.bw --binSize 50 --extendReads 50 -p 4
```

```bash
# CPM normalization
$ bamCoverage -b /home/Li/DataBase/SSC/rmblack_bam/merged_files/SSC_K27me3_merged.bam -o /home/Li/DataBase/SSC/bigwig/SSC_K27me3_merged.CPM.bw --binSize 50 --extendReads 50 -p 4 --normalizeUsing CPM

$ bamCoverage -b /home/Li/DataBase/SSC/rmblack_bam/merged_files/SSC_WCE_merged.bam -o /home/Li/DataBase/SSC/bigwig/SSC_WCE_merged.CPM.bw --binSize 50 --extendReads 50 -p 4 --normalizeUsing CPM
```

## 11. peakcall, RiP calculation and so on
