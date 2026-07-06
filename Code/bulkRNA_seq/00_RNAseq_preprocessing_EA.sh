#!/bin/bash
#
# CREATED USING THE BIOHPC PORTAL on Sat Oct 04 2025 13:27:26 GMT-0500 (Central Daylight Time)
#
# This file is batch script used to run commands on the BioHPC cluster.
# The script is submitted to the cluster using the SLURM `sbatch` command.
# Lines starting with # are comments, and will not be run.
# Lines starting with #SBATCH specify options for the scheduler.
# Lines that do not start with # or #SBATCH are commands that will run.

# Name for the job that will be visible in the job queue and accounting tools.
#SBATCH --job-name EA_RNA_seq_1

# Name of the SLURM partition that this job should run on.
#SBATCH -p 512GB       # partition (queue)
# Number of nodes required to run this job
#SBATCH -N 1

# Memory (RAM) requirement/limit in MB.
#SBATCH --mem 501760      # Memory Requirement (MB)

# Time limit for the job in the format Days-H:M:S
# A job that reaches its time limit will be cancelled.
# Specify an accurate time limit for efficient scheduling so your job runs promptly.
#SBATCH -t 1-2:0:00

# The standard output and errors from commands will be written to these files.
# %j in the filename will be replace with the job number when it is submitted.
#SBATCH -o job_%j.out
#SBATCH -e job_%j.err

# Send an email when the job status changes, to the specfied address.
#SBATCH --mail-type ALL
#SBATCH --mail-user emilija.aleksandrovic@utsouthwestern.edu

# COMMAND GROUP 1: directories
RAW_DIR=/project/pathology/SiZhang_lab/shared/Active_Projects/EA_WM_RNA-seq/01.RawData

# List of sample folders (adjust to your folder names)
SAMPLES=(
Reg_C2I_1 Reg_C2I_2 Reg_C2I_3 Reg_C2P_1 Reg_C2P_2 Reg_C2P_3 Reg_NI_1 Reg_NI_2 Reg_NI_3 Reg_NP_1 Reg_NP_2 Reg_NP_3 ULA_C2I_1 ULA_C2I_2 ULA_C2I_3 ULA_C2P_1 ULA_C2P_2 ULA_C2P_3 ULA_NI_1 ULA_NI_2 ULA_NI_3 ULA_NP_1 ULA_NP_2 ULA_NP_3
# ... add all your sample folders here
)

PROCESSED_DIR=/project/pathology/SiZhang_lab/shared/Active_Projects/EA_RNA_seq_processed
mkdir -p "$PROCESSED_DIR"
cd "$PROCESSED_DIR"

# COMMAND GROUP 2
source $(conda info --base)/etc/profile.d/conda.sh
conda activate rna_seq_env

# COMMAND GROUP 3
echo "==== ENVIRONMENT CHECK ===="
which trimmomatic
which hisat2
which samtools
which featureCounts
which python
conda info --envs | grep "*"
echo "==========================="

# COMMAND GROUP 4: Samples
SAMPLES=(
Reg_C2I_3 Reg_C2I_2 Reg_C2I_1
Reg_C2P_3 Reg_C2P_2 Reg_C2P_1
Reg_NI_3  Reg_NI_2  Reg_NI_1
Reg_NP_3  Reg_NP_2  Reg_NP_1
ULA_C2I_3 ULA_C2I_2 ULA_C2I_1
ULA_C2P_3 ULA_C2P_2 ULA_C2P_1
ULA_NI_3  ULA_NI_2  ULA_NI_1
ULA_NP_3  ULA_NP_2  ULA_NP_1
)

# COMMAND GROUP 5: Trimmomatic
echo "==== RUNNING TRIMMOMATIC ===="
for sample in "${SAMPLES[@]}"; do
    R1="$RAW_DIR/$sample/${sample}_1.fq.gz"
    R2="$RAW_DIR/$sample/${sample}_2.fq.gz"
    trimmomatic PE -threads 8 -phred33 \
        "$R1" "$R2" \
        "${sample}_R1_paired.fastq.gz" "${sample}_R1_unpaired.fastq.gz" \
        "${sample}_R2_paired.fastq.gz" "${sample}_R2_unpaired.fastq.gz" \
        ILLUMINACLIP:/home2/s438978/RNA_seq/references/novogene_adapters.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 &
done
wait
echo "==== TRIMMOMATIC DONE ===="


# COMMAND GROUP 6: HISAT2 Alignment in Parallel
GENOME_INDEX=/project/pathology/SiZhang_lab/shared/UTSW_CITE-seq_RAW/Ref_Ms/mm10-3.0.0_FPs/fasta/genome_index
echo "==== RUNNING HISAT2 ===="
for sample in "${SAMPLES[@]}"; do
    echo "Aligning sample: $sample"
    hisat2 -p 8 -x $GENOME_INDEX \
        -1 ${sample}_R1_paired.fastq.gz \
        -2 ${sample}_R2_paired.fastq.gz \
        -S ${sample}_paired_aligned.sam &
done
wait
echo "==== HISAT2 DONE ===="


# COMMAND GROUP 7: Convert, Sort, and Index Alignment Files

echo "==== RUNNING SAMTOOLS ===="
for sample in "${SAMPLES[@]}"; do
    echo "Processing SAMtools for sample: $sample"
    samtools view -@ 8 -bS ${sample}_paired_aligned.sam | \
        samtools sort -@ 8 -o ${sample}_sorted.bam && \
        samtools index ${sample}_sorted.bam &
done
wait
echo "==== SAMTOOLS DONE ===="


# COMMAND GROUP 8: Quantification with featureCounts

GTF_FILE=/project/pathology/SiZhang_lab/shared/UTSW_CITE-seq_RAW/Ref_Ms/mm10-3.0.0_FPs/genes/genes.gtf
echo "==== RUNNING featureCounts ===="
for sample in "${SAMPLES[@]}"; do
    echo "Quantifying sample: $sample"
    featureCounts -T 8 -p -a $GTF_FILE \
        -o ${sample}_counts.txt \
        ${sample}_sorted.bam &
done
wait
echo "==== featureCounts DONE ===="

echo "==== ALL STEPS COMPLETED ===="

# END OF SCRIPT