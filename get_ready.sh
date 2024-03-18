echo "Starting Genome Index..."
STAR --runMode genomeGenerate --genomeDir genome/ --genomeFastaFiles genome/Ovis_aries_rambouillet.ARS-UI_Ramb_v2.0.dna.primary_assembly.22.fa --outFileNamePrefix genome/chr22_oar_

echo "Starting STAR Alignment..."
STAR --runThreadN 8 --genomeDir genome/ --readFilesIn fastq/texel_pit_1.fastq fastq/texel_pit_2.fastq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/texel_pit_
STAR --runThreadN 8 --genomeDir genome/ --readFilesIn fastq/texel_colon_1.fastq fastq/texel_colon_2.fastq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/texel_colon_
