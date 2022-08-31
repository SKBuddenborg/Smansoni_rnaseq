## Read counts to matrix
- Use HTSeq to output a read count matrix for features
- Where "samples.txt" is a one-column list containing the prefix of the merged bam filenames

```bash
module load cufflinks/2.2.1--py36_2
gffread /lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/01_REF/SM_V9_31Aug21.gff -T -o /lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/01_REF/SM_V9_31Aug21.gtf

module load htseq/0.11.2--py27h637b7d7_1
for i in `cut -f1 /lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/05_ANALYSIS/samples.txt`; do bsub -q normal -n 4 -R "span[hosts=1] select[mem>100000] rusage[mem=100000]" -M 100000 -J $i'htseq' -o $i'htseq.o' -e $i'htseq.e' \
htseq-count -f bam -r pos -s yes -m intersection-nonempty --nonunique=all --stranded=no -t exon -i gene_id -o $i'_readcounts' \
/lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/04_MAPPING/merged_bams/$i'_Aligned.sortedByCoord.merged.bam' \
/lustre/scratch118/infgen/team333/skb/Smansoni_rnaseq/01_REF/SM_V9_31Aug21.gtf; done

```
