# Mapping

Author: Sarah Buddenborg, skb[at]sanger.ac.uk

## Contents
- Project setup
- Prepare reference genome
- Raw sequence data
- Metadata
- Trimming of the raw reads
- Merge bam file from multiple lanes
- Mapping

## Project setup

```bash
# Set working directory
cd /lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq
WORKING_DIR=/lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq

mkdir 00_SCRIPTS  01_REF  02_RAW  03_TRIMMING  04_MAPPING  05_ANALYSIS

```

## Prepare reference genome

```bash
# Make STAR genome index for Mapping
cd ${WORKING_DIR}/01_REF
/lustre/scratch118/infgen/team333/skb/software/STAR-2.7.9a/source/STAR \
--runMode genomeGenerate \
--runThreadN 6 \
--genomeDir /lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/01_REF \
--genomeFastaFiles /nfs/users/nfs_s/skb/ENA/SM_V9_ENA.fa \
--sjdbGTFfile /nfs/users/nfs_s/skb/SM_V9_16Mar.gtf \
--genomeSAindexNbases 10 \
--limitGenomeGenerateRAM 90000000000

```


## Raw sequencing data
- Sequencing was performed at the Wellcome Sanger Institute, Hinxton UK
- Raw sequencing data was symlink'd to the lustre environment into the directory 02_RAW

```bash
cd ${WORKING_DIR}/02_RAW
pf data -t study -i 5837 --symlink

```
---

## Metadata
- Googlesheet: https://docs.google.com/spreadsheets/d/1w-RX5himr6UNS5NTCBAQsxYM-q2qdkcL/edit#gid=1245866365

---
## Trimming of raw reads

```bash
cd ${WORKING_DIR}/03_TRIMMING
module load trimgalore/v0.4.4

```
where "samples.txt" is a one-column list containing the prefix of fastq file names
and where "trim.sh" is:

```bash
#!/bin/bash
for i in `cut -f1 /lustre/scratch118/infgen/team133/skb/RNAseq/samples.txt`; do trim_galore \
--quality 20 \
--phred33 \
--retain_unpaired \
--fastqc \
--fastqc_args "--outdir /lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/03_TRIMMING" \
--output_dir /lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/03_TRIMMING \
--paired $i'_1.fastq.gz' $i'_2.fastq.gz' \
; done

bsub -q long -n 4 -R "span[hosts=1] select[mem>50000] rusage[mem=50000]" -M 50000 -J trim -o trim.o -e trim.e ./trim.sh

```

## Mapping
- Paired-end (PE) mapping of reads

```bash
cd ${WORKING_DIR}/04_MAPPING
for i in `cut -f1 ./samples.txt`; do bsub -q normal -n 4 -R "span[hosts=1] select[mem>100000] rusage[mem=100000]" -M 100000 -J $i -o $i'.o' -e $i'.e' \
/lustre/scratch118/infgen/team333/skb/software/STAR-2.7.9a/source/STAR \
--runThreadN 4 \
--genomeDir /lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/01_REF \
--readFilesIn /lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/03_TRIMMING/$i'_1_val_1.fq' /lustre/scratch118/infgen/team133/skb/RNAseq/trimmed/$i'_2_val_2.fq' \
--alignIntronMin 10 \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 90000000000 \
--outFileNamePrefix $i'_'; done

```

###FIX: some samples did not map because they require more RAM than originally requested above
where "samples2.txt" is a one-column list containing the prefix of fastq file names that require more memory for mapping

```bash
# use hugemem queue
for i in `cut -f1 ./samples2.txt`; do bsub -q normal -n 4 -R "span[hosts=1] select[mem>100000] rusage[mem=100000]" -M 100000 -J $i -o $i'.o' -e $i'.e' \
/lustre/scratch118/infgen/team333/skb/software/STAR-2.7.9a/source/STAR \
--runThreadN 4 \
--genomeDir /lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/01_REF \
--readFilesIn /lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/03_TRIMMING/$i'_1_val_1.fq' /lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/03_TRIMMING/$i'_2_val_2.fq' \
--alignIntronMin 10 \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix $i'_'; done

```

## Merge bam files from the same sample run over two lanes
- The same pooled libraries for all samples were run over two Novaseq lanes so these bam files need to be merged for the same samples

where "samples.txt" is a three-column list containing the prefix of the two bam file names and the new merged filename
and where "merge.sh" is:

```bash
#!/bin/bash
samtools merge Eggs_R1_Aligned.sortedByCoord.merged.bam 35207_2#3_Aligned.sortedByCoord.out.bam 35101_2#3_Aligned.sortedByCoord.out.bam;

cd /lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/04_MAPPING/merge
module load samtools/1.9
bsub -q normal -n 12 -R "span[hosts=1] select[mem>100000] rusage[mem=100000]" -M 100000 -J merge -o merge.o -e merge.e ./merge.sh


```

```bash
cd ${WORKING_DIR}/04_MAPPING
module load samtools/1.9
