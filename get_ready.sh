echo "Starting Genome Index..."
STAR --runMode genomeGenerate --genomeDir genome/ --genomeFastaFiles genome/Ovis_aries_rambouillet.ARS-UI_Ramb_v2.0.dna.primary_assembly.22.fa --outFileNamePrefix chr22_oar

echo "Starting STAR Alignment..."
STAR --runThreadN 8 --genomeDir genome/ --readFilesIn fastq/texel_pit_1.fastq fastq/texel_pit_2.fastq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/texel_pit
STAR --runThreadN 8 --genomeDir genome/ --readFilesIn fastq/texel_colon_1.fastq fastq/texel_colon_2.fastq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/texel_colon
