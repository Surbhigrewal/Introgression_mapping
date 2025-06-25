#!/bin/bash

# ------------------------
# Introgression Mapping Pipeline - Main Bash Workflow
# ------------------------

# Set the sample name and parent identifiers
prefix=Introgression_line_name          # e.g., Mut1
wheat_parent_1=Wheat_parent_1_name      # e.g., CS_mut
wheat_parent_2=Wheat_parent_2_name      # e.g., Paragon_mut

# ------------------------
# Step 1: Build combined reference index of wheat and wild relative genome fasta sequence
# ------------------------
hisat2-build -p 80 combinedref_CSMut.fasta combinedref_CSMut

# ------------------------
# Step 2: Adapter trimming using Trimmomatic
# ------------------------
trimmomatic PE -threads 96 \
  ${prefix}/${prefix}_1.fastq.gz ${prefix}/${prefix}_2.fastq.gz \
  ${prefix}/${prefix}_1.trimmed.fastq.gz ${prefix}/${prefix}_unpaired_1.trimmed.fastq.gz \
  ${prefix}/${prefix}_2.trimmed.fastq.gz ${prefix}/${prefix}_unpaired_2.trimmed.fastq.gz \
  ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads \
  SLIDINGWINDOW:4:15 MINLEN:40 CROP:150

# ------------------------
# Step 3: Align trimmed reads to the combined reference using HISAT2
# ------------------------
hisat2 -p 96 -x combinedref_CSMut \
  -1 ${prefix}/${prefix}_1.trimmed.fastq.gz \
  -2 ${prefix}/${prefix}_2.trimmed.fastq.gz \
  -S ${prefix}/${prefix}.sam \
  --no-spliced-alignment --no-unal &> ${prefix}/${prefix}.log

# ------------------------
# Step 4: Filter for primary alignments, map quality ≥10, then sort
# ------------------------
samtools view -h ${prefix}/${prefix}.sam | \
  awk '$1 ~ /^@/ || $2 == 65 || $2 == 129 || $2 == 67 || $2 == 131 || $2 == 113 || $2 == 177 || $2 == 81 || $2 == 161 || $2 == 163 || $2 == 83 || $2 == 97 || $2 == 145 || $2 == 99 || $2 == 147 || $2 == 137 || $2 == 73 {print $0}' | \
  samtools view -u -q 10 - | samtools sort -o ${prefix}/${prefix}_filt_srt.bam -
samtools index -c ${prefix}/${prefix}_filt_srt.bam

# ------------------------
# Step 5: Remove PCR duplicates and reassign tags with Picard
# ------------------------
picard -Xmx128g MarkDuplicates \
  -I ${prefix}/${prefix}_filt_srt.bam \
  -M ${prefix}/${prefix}-metrics.txt \
  -REMOVE_DUPLICATES TRUE \
  -ASSUME_SORTED TRUE \
  -MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000 \
  -VALIDATION_STRINGENCY LENIENT \
  -O /dev/stdout | \
  picard -Xmx16g SetNmMdAndUqTags \
  -I /dev/stdin \
  -O ${prefix}/${prefix}_filt_srt_rmdup.bam \
  -R combinedref_CSMut.fasta
samtools index -c ${prefix}/${prefix}_filt_srt_rmdup.bam

# ------------------------
# Step 6: Create 1 Mb genomic windows from reference index
# ------------------------
samtools faidx combinedref_CSMut.fasta
awk '{OFS="\t"; print $1,$2}' combinedref_CSMut.fasta.fai > CSMut_genome_file.txt
bedtools makewindows -w 1000000 -g CSMut_genome_file.txt > CSMut_1Mb_windows.bed

# ------------------------
# Step 7: Count reads per 1 Mb window
# ------------------------
hts_nim_tools count-reads CSMut_1Mb_windows.bed ${prefix}/${prefix}_filt_srt_rmdup.bam | \
  sort -k1,1 -k2,2n | \
  awk '{OFS="\t"; print}' > ${prefix}/${prefix}_cov_1Mb_windows.tsv

# ------------------------
# Step 8: Compute coverage deviation against parental controls
# ------------------------
# ⚠️ Make sure you have already created:
# ${wheat_parent_1}_cov_1Mb_windows.tsv
# ${wheat_parent_2}_cov_1Mb_windows.tsv
python3 cov_deviation_new_high.py \
  ${prefix}/${prefix}_cov_1Mb_windows.tsv \
  ${wheat_parent_1}_cov_1Mb_windows.tsv \
  ${wheat_parent_2}_cov_1Mb_windows.tsv \
  1Mb ${prefix}

# ------------------------
# Step 9: Plot genome-wide coverage deviation and highlight introgressions
# ------------------------
Rscript Chr.plot.R ${prefix}
