# Mapping

Author: Sarah Kay Buddenborg, skb[at]sanger.ac.uk

## Contents
- Paired-end (PE) mapping of reads to genome using STAR v2.7.9a
- Merge genome mapped bam files from multiple lanes

## STAR mapping
Set working directory:
```bash
cd ${WORKING_DIR}/04_MAPPING
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
/data/pam/team333/skb/scratch/software/STAR-2.7.9a/source/STAR \
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

```bash
for i in `cut -f1 ${WORKING_DIR}/samples.txt`; do bsub -q normal -n 4 -R "span[hosts=1] select[mem>100000] rusage[mem=100000]" -M 100000 -J $i -o $i'.o' -e $i'.e' \
/data/pam/team333/skb/scratch/software/STAR-2.7.9a/source/STAR \
--runThreadN 4 \
--genomeDir ${WORKING_DIR}/01_REF/v10 \
--readFilesIn ${WORKING_DIR}/03_TRIMMING/$i'_1_val_1.fq' ${WORKING_DIR}/03_TRIMMING/$i'_2_val_2.fq' \
--alignIntronMin 10 \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 90000000000 \
--outFileNamePrefix $i'_'; done
```

FIX: some samples did not map because they require more memory for bam sorting or longer runtime queues than originally requested above.

where "samples_hugemem_.txt" and "samples_long.txt" are one-column lists containing the prefix of fastq file names that require more memory for bam sorting or longer runtime queue for mapping, respectively. 

```bash
# hugemem queue
for i in `cut -f1 ${WORKING_DIR}/04_MAPPING/samples_hugemem.txt`; do bsub -q hugemem -n 1 -R "span[hosts=1] select[mem>200000] rusage[mem=200000]" -M 200000 -J $i -o $i'.o' -e $i'.e' /data/pam/team333/skb/scratch/software/STAR-2.7.9a/source/STAR --runThreadN 6 --genomeDir ${WORKING_DIR}/01_REF/v10 --readFilesIn ${WORKING_DIR}/03_TRIMMING/$i'_1_val_1.fq' ${WORKING_DIR}/03_TRIMMING/$i'_2_val_2.fq' --alignIntronMin 10 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 190000000000  --outFileNamePrefix $i'_'; done
```

```bash
# long queue
for i in `cut -f1 ${WORKING_DIR}/04_MAPPING/samples_long.txt`; do bsub -q long -n 6 -R "span[hosts=1] select[mem>100000] rusage[mem=100000]" -M 100000 -J $i -o $i'.o' -e $i'.e' /data/pam/team333/skb/scratch/software/STAR-2.7.9a/source/STAR --runThreadN 6 --genomeDir ${WORKING_DIR}/01_REF/v10 --readFilesIn ${WORKING_DIR}/03_TRIMMING/$i'_1_val_1.fq' ${WORKING_DIR}/03_TRIMMING/$i'_2_val_2.fq' --alignIntronMin 10 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 90000000000  --outFileNamePrefix $i'_'; done
```

For samples not run for unknown reason:
```bash
for i in `cut -f1 ${WORKING_DIR}/04_MAPPING/samples_notrun.txt`; do bsub -q normal -n 6 -R "span[hosts=1] select[mem>100000] rusage[mem=100000]" -M 100000 -J $i -o $i'.o' -e $i'.e' \
/data/pam/team333/skb/scratch/software/STAR-2.7.9a/source/STAR \
--runThreadN 6 \
--genomeDir ${WORKING_DIR}/01_REF/v10 \
--readFilesIn ${WORKING_DIR}/03_TRIMMING/$i'_1_val_1.fq' ${WORKING_DIR}/03_TRIMMING/$i'_2_val_2.fq' \
--alignIntronMin 10 \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 90000000000 \
--outFileNamePrefix $i'_'; done
```

For samples "could not open read files":
```bash
for i in `cut -f1 ${WORKING_DIR}/04_MAPPING/samples_noinput.txt`; do bsub -q normal -n 6 -R "span[hosts=1] select[mem>100000] rusage[mem=100000]" -M 100000 -J $i -o $i'.o' -e $i'.e' \
/data/pam/team333/skb/scratch/software/STAR-2.7.9a/source/STAR \
--runThreadN 6 \
--genomeDir ${WORKING_DIR}/01_REF/v10 \
--readFilesIn ${WORKING_DIR}/03_TRIMMING/$i'_1_val_1.fq' ${WORKING_DIR}/03_TRIMMING/$i'_2_val_2.fq' \
--alignIntronMin 10 \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 90000000000 \
--outFileNamePrefix $i'_'; done
```

Output counts of total reads uniquely mapped:
```bash
grep 'Uniquely mapped reads number' *Log.final.out > mapped_counts.txt```
```

Output counts of total reads filtered as input into STAR:
```bash
grep 'Number of input reads' *Log.final.out > filtered_counts.txt
```

## Merge bam files from the same sample run over two lanes
- The same pooled libraries for all samples were run over two Novaseq lanes so these bam files need to be merged for the same samples

```bash
cd ${WORKING_DIR}/04_MAPPING
mkdir /merged
cd ${WORKING_DIR}/04_MAPPING/merged
module load samtools/1.9
```
where "samples_merge.txt" is a tab-delimited file containing the output filename, bam filename from lane 1, and bam filename from lane 2

```bash
while read -r SAMPLE FILE1 FILE2; do bsub -q normal -n 2 -R "span[hosts=1] select[mem>100000] rusage[mem=100000]" -M 100000 -J ${SAMPLE} -o ${SAMPLE}.o -e ${SAMPLE}.e samtools merge ${SAMPLE} ${WORKING_DIR}/04_MAPPING/${FILE1} ${WORKING_DIR}/04_MAPPING/${FILE2}; done < samples_merge.txt
```

## Kallisto mapping to the transcriptome
Create genome index file for Kallisto mapping
```bash
cd ${WORKING_DIR}/01_REF/v10 
module load kallisto/0.46.2--h4f7b962_1
nano kallisto_index.sh
kallisto index --make-unique -I index.fa.idx index.fa

bsub -q normal -n 6 -R "span[hosts=1] select[mem>50000] rusage[mem=50000]" -M 50000 -J index -o index.o -e index.e kallisto index --make-unique -i SM_V10.fa.idx SM_V10.fa
```
Kallisto pseudomapping
```bash
for i in `cut -f1 ${WORKING_DIR}/04_MAPPING/kallisto/samples_hash.txt`; \
do bsub -q normal -n 10 -R "span[hosts=1] select[mem>20000] rusage[mem=20000]" -M 20000 -J $i -o $i'.o' -e $i'.e' \
kallisto quant \
-i ${WORKING_DIR}/01_REF/v10/SM_V10.fa.idx \
-o ${WORKING_DIR}/04_MAPPING/kallisto/$i \
-t 10 \
-g ${WORKING_DIR}/01_REF/v10/SM_V10/SM_V10.gtf \
-c ${WORKING_DIR}/01_REF/v10/chrNameLength.txt \
-b 100 \
${WORKING_DIR}/03_TRIMMING/35207_2#$i'_1_val_1'.fq \
${WORKING_DIR}/03_TRIMMING/35207_2#$i'_2_val_2'.fq \
${WORKING_DIR}/03_TRIMMING/35101_2#$i'_1_val_1'.fq \
${WORKING_DIR}/03_TRIMMING/35101_2#$i'_2_val_2'.fq \
;done

```
## TEST RUN
```bash
head -4000 ${WORKING_DIR}/03_TRIMMING/35101_2#62_1_val_1.fq > ${WORKING_DIR}/04_MAPPING/kallisto/test.fq
head -4000 ${WORKING_DIR}/03_TRIMMING/35101_2#62_2_val_2.fq > ${WORKING_DIR}/04_MAPPING/kallisto/test_2.fq
```

```bash
bsub -q normal -n 1 -R "span[hosts=1] select[mem>100000] rusage[mem=100000]" -M 100000 -J test -o test.o -e test.e \
kallisto quant \
-i ${WORKING_DIR}/01_REF/v10/SM_V10.fa.idx \
-o ${WORKING_DIR}/04_MAPPING/kallisto \
-t 1 \
-g ${WORKING_DIR}/01_REF/v10/SM_V10.gtf \
-c ${WORKING_DIR}/01_REF/v10/chrNameLength.txt \
-b 10 \
${WORKING_DIR}/04_MAPPING/kallisto/test.fq \
${WORKING_DIR}/04_MAPPING/kallisto/test_2.fq
```

## STAR bam files to read counts
- Use StringTie to generate a read count matrices
- Where "samples.txt" is a one-column list containing the prefix of the merged bam filenames
- Remove superfluous lines from output file using sed

```bash
module load stringtie/2.1.4--h7e0af3c_0
for i in `cut -f1 ${WORKING_DIR}/04_MAPPING/merged/samples_merge.txt`; do bsub -q yesterday -n 1 -R "span[hosts=1] select[mem>50000] rusage[mem=50000]" -M 50000 -J $i -o $i'.o' -e $i'.e' \
stringtie ${WORKING_DIR}/04_MAPPING/merged/$i \
-G ${WORKING_DIR}/01_REF/v10/SM_V10.gtf \
-p 1 \
-A ${WORKING_DIR}/05_ANALYSIS/$i'_STAR_gene_abund.out'\
; done