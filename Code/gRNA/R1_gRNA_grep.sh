# load modules
module load python/3.8.x-anaconda bbmap/38.46

# to activate Conda env
conda activate ISS

# install fastp using : conda install -c bioconda fastp

# write path to input (where your files are) and output (where you want to output files) directories. 
datadir=/project/pathology/SiZhang_lab/shared/Active_Projects/EA_ova-gRNA/Fastq
output=/project/pathology/SiZhang_lab/shared/Active_Projects/EA_ova-gRNA/Fastq/gRNA

# 1. run fastq QC and export html reports

fastp -i $datadir/OvCa-1_S1_L001_R1_001.fastq.gz -I $datadir/OvCa-1_S1_L001_R2_001.fastq.gz -o $datadir/out.OvCa-1_S1_R1.fastq.gz -O $datadir/out.OvCa-1_S1_R2.fastq.gz -h OvCa-1_S1_report.html
fastp -i $datadir/OvCa-2_S2_L001_R1_001.fastq.gz -I $datadir/OvCa-2_S2_L001_R2_001.fastq.gz -o $datadir/out.OvCa-2_S2_R1.fastq.gz -O $datadir/out.OvCa-2_S2_R2.fastq.gz -h OvCa-2_S2_report.html

fastp -i $datadir/OvCa-1_S1_L001_R1_001.fastq.gz -o $datadir/out.OvCa-1_S1_R1_only.fastq.gz -h $datadir/OvCa-1_S1_report_R1.html
fastp -i $datadir/OvCa-1_S1_L001_R2_001.fastq.gz -o $datadir/out.OvCa-1_S1_R2_only.fastq.gz -h $datadir/OvCa-1_S1_report_R2.html

# 2. select reads with > Q25 score 
bbduk.sh in1=$datadir/OvCa-1_S1_L001_R1_001.fastq.gz \
         out1=$output/filtered_S1_R1.fastq.gz \
         maq=25

bbduk.sh in1=$datadir/OvCa-2_S2_L001_R1_001.fastq.gz \
         out1=$output/filtered_S2_R1.fastq.gz \
         maq=25

# 3. Extract gRNA from R1
zcat $output/filtered_S1_R1.fastq.gz | awk 'NR%4==2 {
    gRNA_index = index($0, "TGCTGTTTCCAGCATAGCTCTGAAAC");
    if (gRNA_index > 0) {
        gRNA = substr($0, gRNA_index + length("TGCTGTTTCCAGCATAGCTCTGAAAC"), 21);
        print gRNA
    } else {
        print "NA"
    }
}' > $output/extracted_S1_gRNA.csv 

zcat $output/filtered_S2_R1.fastq.gz | awk 'NR%4==2 {
    gRNA_index = index($0, "TGCTGTTTCCAGCATAGCTCTGAAAC");
    if (gRNA_index > 0) {
        gRNA = substr($0, gRNA_index + length("TGCTGTTTCCAGCATAGCTCTGAAAC"), 21);
        print gRNA
    } else {
        print "NA"
    }
}' > $output/extracted_S2_gRNA.csv 

# filter extracted gRNA

awk '!/NA/ && !/^(A+|T+|C+|G+)$/ && length($0) == 21' $output/extracted_S1_gRNA.csv > $output/extracted_S1_gRNA_filtered.csv
awk '!/NA/ && !/^(A+|T+|C+|G+)$/ && length($0) == 21' $output/extracted_S2_gRNA.csv > $output/extracted_S2_gRNA_filtered.csv

# Calculate Uni gRNA, calculate freq and plot a histgram
cat $output/extracted_S1_gRNA_filtered.csv | tr -d '\r' | sed 's/ *$//' | sort | uniq -c | sort -k1,1nr > $output/R1_filtered_Uni_extracted_S1_gRNA_filtered.csv_frequencies.txt
cat $output/extracted_S2_gRNA_filtered.csv | tr -d '\r' | sed 's/ *$//' | sort | uniq -c | sort -k1,1nr > $output/R1_filtered_Uni_extracted_S2_gRNA_filtered.csv_frequencies.txt

# plot in python using "plot_gRNA_freq.py"

# calculate coverage in python using "Matching_analysis.py"

