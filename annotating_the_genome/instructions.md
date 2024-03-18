###########
# STEP 01 #
###########

-> Review what has happened during the break!

###########
# STEP 02 #
###########

-> fastqc

cd /workspace/ensgb-livestock-genomics-2024/annotating_the_genome;
mkdir fastq/qc;
fastqc -t 2 -o fastq/qc/ fastq/texel_colon_1.fastq fastq/texel_colon_2.fastq;
fastqc -t 2 -o fastq/qc/ fastq/texel_pit_1.fastq fastq/texel_pit_2.fastq;

-> Review what fastqc produced!

cd /workspace/ensgb-livestock-genomics-2024/annotating_the_genome/fastq/qc;
for i in `ls *.zip`; do unzip $i; done;

###########
# STEP 03 #
###########



###########
# STEP 04 #
###########

###########
# STEP 05 #
###########

###########
# STEP 05 #
###########

###########
# STEP 07 #
###########

###########
# STEP 08 #
###########

###########
# STEP 09 #
###########

###########
# STEP 10 #
###########


