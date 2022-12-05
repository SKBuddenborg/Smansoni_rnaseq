# Mapping

Author: Sarah Kay Buddenborg, skb[at]sanger.ac.uk

## Contents
- Project setup
- Prepare reference genome
- Raw sequence data
- Metadata
- Trimming of the raw reads
- Merge bam file from multiple lanes
- Mapping

## Project setup
Set working directory:

```bash
cd /lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq

WORKING_DIR=/lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq

mkdir 00_SCRIPTS  01_REF  02_RAW  03_TRIMMING  04_MAPPING  05_ANALYSIS
```

Collect pathfind suppelementary info for samples to match lane ID with sample name
where "sanger_sample_id.txt" is a single column list of the internal Sanger Sample ID numbers.

```bash
cd ${WORKING_DIR}
module load pathfind

pf supplementary --id ./sanger_sample_id.txt --type file --file-id-type samples > samples_info.txt
```



## Prepare reference genome
Convert GFF3 to GTF:
```bash
cd ${WORKING_DIR}/01_REF/v10
module load cufflinks/2.2.1--py36_2
gffread SM_V10.gff3 -T -o SM_V10.gtf
```

where "genomeGen.sh" is:

```bash
#!/bin/bash
/lustre/scratch118/infgen/team333/skb/software/STAR-2.7.9a/source/STAR \
--runMode genomeGenerate \
--runThreadN 6 \
--genomeDir ${WORKING_DIR}/01_REF/v10 \
--genomeFastaFiles ${WORKING_DIR}/01_REF/v10/SM_V10.fa \
--sjdbGTFfile ${WORKING_DIR}/01_REF/v10/SM_V10.gtf \
--genomeSAindexNbases 10 \
--limitGenomeGenerateRAM 90000000000
```

```bash
bsub -q normal -n 6 -R "span[hosts=1] select[mem>100000] rusage[mem=100000]" -M 100000 -J genomeGen -o genomeGen.o -e genomeGen.e ./genomeGen.sh
```
---

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
for i in `cut -f1 ${WORKING_DIR}/samples.txt`; do trim_galore \
--quality 20 \
--phred33 \
--retain_unpaired \
--fastqc \
--fastqc_args "--outdir ${WORKING_DIR}/03_TRIMMING" \
--output_dir ${WORKING_DIR}/03_TRIMMING \
--paired $i'_1.fastq.gz' $i'_2.fastq.gz' \
; done
```

```bash
bsub -q long -n 4 -R "span[hosts=1] select[mem>50000] rusage[mem=50000]" -M 50000 -J trim -o trim.o -e trim.e ./trim.sh
```

## Mapping
- Paired-end (PE) mapping of reads using STAR v2.7.9a
```bash
cd ${WORKING_DIR}/04_MAPPING
```

```bash
for i in `cut -f1 ${WORKING_DIR}/samples.txt`; do bsub -q normal -n 4 -R "span[hosts=1] select[mem>100000] rusage[mem=100000]" -M 100000 -J $i -o $i'.o' -e $i'.e' \
/lustre/scratch118/infgen/team333/skb/software/STAR-2.7.9a/source/STAR \
--runThreadN 4 \
--genomeDir ${WORKING_DIR}/01_REF/v10 \
--readFilesIn ${WORKING_DIR}/03_TRIMMING/$i'_1_val_1.fq' ${WORKING_DIR}/03_TRIMMING/$i'_2_val_2.fq' \
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
--genomeDir ${WORKING_DIR}/01_REF \
--readFilesIn ${WORKING_DIR}/03_TRIMMING/$i'_1_val_1.fq' ${WORKING_DIR}/03_TRIMMING/$i'_2_val_2.fq' \
--alignIntronMin 10 \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix $i'_'; done
```

## Merge bam files from the same sample run over two lanes
- The same pooled libraries for all samples were run over two Novaseq lanes so these bam files need to be merged for the same samples

```bash
cd ${WORKING_DIR}/04_MAPPING/merge
module load samtools/1.9
```
where "samples_merge.txt" is a three-column list containing the prefix of the two bam file names and the new merged filename
and where "merge.sh" is:

```bash
#!/bin/bash

```

```bash
bsub -q normal -n 12 -R "span[hosts=1] select[mem>100000] rusage[mem=100000]" -M 100000 -J merge -o merge.o -e merge.e ./merge.sh
```

## Read counts to matrix for genes
- Use HTSeq to output a read count matrix for coding gene features
- Where "samples.txt" is a one-column list containing the prefix of the merged bam filenames
- Remove superfluous lines from output file using sed

```bash
module load cufflinks/2.2.1--py36_2
gffread /lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/01_REF/SM_V9_31Aug21.gff -T -o /lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/01_REF/SM_V9_31Aug21.gtf

module load htseq/0.11.2--py27h637b7d7_1
for i in `cut -f1 /lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/05_ANALYSIS/samples.txt`; do bsub -q normal -n 4 -R "span[hosts=1] select[mem>100000] rusage[mem=100000]" -M 100000 -J $i'htseq' -o $i'htseq.o' -e $i'htseq.e' \
htseq-count -f bam -r pos -s no -m intersection-nonempty --nonunique=all -t exon -i gene_id \
/lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/04_MAPPING/merged_bams/$i'_Aligned.sortedByCoord.merged.bam' \
/lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/01_REF/SM_V9_31Aug21.gtf; done

for i in `cut -f1 /lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/05_ANALYSIS/samples.txt`; do cat $i'htseq.o' | egrep "^Smp" $i'htseq.o' > $i'_fixhtseq.o'; done 
```

## Read counts to matrix for lncRNAs
- Use HTSeq to output a read count matrix for lncRNA features
- Where "samples.txt" is a one-column list containing the prefix of the merged bam filenames
- Remove superfluous lines from output file using sed

```bash
module load cufflinks/2.2.1--py36_2
gffread /lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/01_REF/lncRNA_April2021.NameChanged.addedType.gff3 -T -o /lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/01_REF/lncRNA_April2021.NameChanged.addedType.gtf

cd /lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/05_ANALYSIS/lncRNA_htseq
module load htseq/0.11.2--py27h637b7d7_1
for i in `cut -f1 /lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/05_ANALYSIS/samples.txt`; do bsub -q normal -n 4 -R "span[hosts=1] select[mem>100000] rusage[mem=100000]" -M 100000 -J $i'htseq' -o $i'htseq.o' -e $i'htseq.e' \
htseq-count -f bam -r pos -s no -m intersection-nonempty --nonunique=all -t exon -i gene_id \
/lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/04_MAPPING/merged_bams/$i'_Aligned.sortedByCoord.merged.bam' \
/lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/01_REF/lncRNA_April2021.NameChanged.addedType.gtf; done

for i in `cut -f1 /lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/05_ANALYSIS/samples.txt`; do cat $i'htseq.o' | egrep "^Smp" $i'htseq.o' > $i'_fixhtseq.o'; done 
```

# !!!!!!!!! X_Eggs_R1 and X_32d_Sporo_R1 have 0 counts for all genes
bsub -q normal -n 6 -R "span[hosts=1] select[mem>100000] rusage[mem=100000]" -M 100000 -J 32d_Sporo_htseq -o 32d_Sporo_htseq.o -e 32d_Sporo_htseq.e \
htseq-count -f bam -r pos -s no -m intersection-nonempty --nonunique=all -t exon -i gene_id \
/lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/04_MAPPING/merged_bams/32d_Sporocysts_R1_Aligned.sortedByCoord.merged.bam \
/lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/01_REF/SM_V9_31Aug21.gtf

bsub -q normal -n 6 -R "span[hosts=1] select[mem>100000] rusage[mem=100000]" -M 100000 -J Eggs_R1_htseq -o Eggs_R1_htseq.o -e Eggs_R1_htseq.e \
htseq-count -f bam -r pos -s no -m intersection-nonempty --nonunique=all -t exon -i gene_id \
/lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/04_MAPPING/merged_bams/Eggs_R1_Aligned.sortedByCoord.merged.bam \
/lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/01_REF/SM_V9_31Aug21.gtf



egrep "^Smp" Eggs_R1_htseq.o > Eggs_R1_fixhtseq.o
egrep "^Smp" 32d_Sporocysts_R1_htseq.o > 32d_Sporocysts_R1_fixhtseq.o

