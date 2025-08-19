#!/bin/bash
module load python/3.10.x-anaconda
conda activate ISS # most packages is installed under ISS environment. 

# use cat command to combined fastqs: cat LARRY_RPI1_UTSW02_L2_CKDL230026701-1A_22CGKWLT3_S1_L004_R1_001.fastq.gz LARRY_L2_RPI2_CKDL240007529-1A_22GVFFLT3_S1_L004_R1_001.fastq.gz > L2_combined_R1.fastq.gz

# write path to input (where your files are) and output (where you want to output files) directories
datadir=/project/pathology/SiZhang_lab/shared/Active_Projects/Ova_LARRY/RPI/Z/L1
output=/project/pathology/SiZhang_lab/shared/Active_Projects/Ova_LARRY/RPI/Z/L1

# 1. Extract cellBC, UMI, Lib-ID, LARRY
zcat L1_combined_S1_L003_R1_001.fastq.gz | awk 'NR%4==2 {
    cell_barcode = substr($0, 1, 16);
    UMI = substr($0, 17, 12);
    lib_id_index = index($0, "TTGCTAGGACCGGCCTTAAAGC");
    if (lib_id_index > 0) {
        lib_id = substr($0, lib_id_index + length("TTGCTAGGACCGGCCTTAAAGC"), 6);
    } else {
        lib_id = "NA";
    }
    LARRY_index = index($0, "CCACGGTGGCGATATCGGATCCAGACAT");
    if (LARRY_index > 0) {
        LARRY = substr($0, LARRY_index + length("CCACGGTGGCGATATCGGATCCAGACAT"), 40);
    } else {
        LARRY = "NA";
    }
    print cell_barcode "\t" UMI "\t" lib_id "\t" LARRY
}' > extracted_cellBC_UMI_LibIDs_LARRY.txt

# 2. Filter out NA and long strech of A/C/C/G
awk '!/NA/ && length($4) == 40 && !($1 ~ /(A{8}|T{8}|G{8}|C{8})/ || $2 ~ /(A{8}|T{8}|G{8}|C{8})/ || $3 ~ /(A{8}|T{8}|G{8}|C{8})/ || $4 ~ /(A{8}|T{8}|G{8}|C{8})/)' extracted_cellBC_UMI_LibIDs_LARRY.txt > extracted_cellBC_UMI_LibIDs_LARRY_filtered.txt
# >> extracted_cellBC_UMI_LibIDs_LARRY_filtered.txt will be used for 03_correctBC_parallel_v3.py

# 3. Generate LARRY whitelist
awk 'length($4) == 40 {print $4}' extracted_cellBC_UMI_LibIDs_LARRY_filtered.txt > filtered_LARRY_only.txt
# then # >>> use "filtered_LARRY_only.txt" as input file for 02_Create_LARRY_WL.py script to generate LARRY_WL.txt
# e.g. python 02_Create_LARRY_WL.py --input_file filtered_LARRY_only.txt --output_file L1_LARRY_WL.txt --freq_output_file L1_LARRY_freq.txt --sensitivity 5


