# Ensembl genome annotation (practicum)

Jose Maria Gonzalez Perez-Silva  
Bioinformatician  
EMBL-EBI

Ensembl RNASeq Practical
========================

The aim of this practical session is to use STAR to align 2 lanes of sheep (texel breed) paired end Illumina RNASeq reads to chromosome 22 of the reference sheep (ramboulliet) ARS-UI_Ramb_v2.0 assembly. We have restricted the analysis to chromosome 22 in order to speed up the alignment process. 

Once the reads have been aligned you can assemble transcript models using Scallop and finally display the transcript models in the Ensembl browser alongside some of our own annotation.

## Overview:

We are going to use "terminal" application in gitpod to execute various commands.

Command-line actions usually start with 4 spaces ("    ") and should be rendered as code blocks in GitHub markdown preview, which you can copy using the "copy" icon top right of each block. Please note that many commands are wrapped across more than one line.

There are questions throughout this document, indicated with a "->", try to answer these, they should help you understand the process.


Step 0: Familiarise yourself with the application
-------------------------------------------------

Go to the terminal (should be in the bottom of your screen)
You should be in cd /workspace/ensgb-livestock-genomics-2024/annotating_the_genome, if not:

    cd /workspace/ensgb-livestock-genomics-2024/annotating_the_genome

To check which directory you are in, please execute:

    pwd

> :warning: At the start of each step, you should be in cd /workspace/ensgb-livestock-genomics-2024/annotating_the_genome so run the `cd /workspace/ensgb-livestock-genomics-2024/annotating_the_genome` command if you get any "file / path not found" errors.
> :warning: You should always see `(annotation)` at the beguining of your line (before the green `gitpod` part). If you can't see it, run `conda activate annotation` and this will prevent many errors.

To check what is in the directories, please execute:

    ls -lh


Step 1: Check the quality of the data
-------------------------------------

This time, the data has been previously checked, of course, and due to the vicissitudes of the practicum, the results will not be great, but it's a good opportunity to practise.

FastQC command to run:

    mkdir fastq/qc;
    fastqc -t 2 -o fastq/qc/ fastq/texel_colon_1.fastq fastq/texel_colon_2.fastq;
    fastqc -t 2 -o fastq/qc/ fastq/texel_pit_1.fastq fastq/texel_pit_2.fastq;

In this way, we are creating a folder to store the results, and then running the software FastQC to check for the many variables that the program validates.

We'll have 2 results per read (so 8 in total): a `.zip` and a `.html`. To extract the zip files and be able to read them, we can execute the following bash one-liner:

    cd /workspace/ensgb-livestock-genomics-2024/annotating_the_genome/fastq/qc;
    for i in `ls *.zip`; do unzip $i; done;

This will take us to the results folder and then execute a loop that will extract one by one all the zip files availables.

The FastQC documentation includes a lot of explanations and guides to understand the results.

> :warning: remember to return to the original folder!!! You can run `cd ../../` or the comand in the Step 0.

Step 2: Index the genome file
-----------------------------

Now that we know the quality of our reads, we need to index the genome file so that STAR can use it. We do this using the STAR genomeGenerate command. 

STAR command to index the genome fasta file:

    STAR --runMode genomeGenerate --genomeDir genome/ --genomeFastaFiles genome/Ovis_aries_rambouillet.ARS-UI_Ramb_v2.0.dna.primary_assembly.22.fa --outFileNamePrefix genome/chr22_oar;

Output will be printed in your terminal, STAR commands complete successfully with the output like "[main] Real time: 10 sec; CPU 15 sec".

* -> 1. What output did this command create? hint: use `ls ./Annotating_the_genome/Genome/` or use the file browser on the left to see what new files have appeared.

* -> 2. Why do we index the genome before aligning to it? 


Step 3: Align the reads to the genome
-------------------------------------

When the genome has been indexed, we can start the alignment. We are using 2 lanes of paired end reads, that gives us 4 files, 2 for each lane containing the 1st (or forward) and 2nd (or reverse) reads respectively.

To test how different tissues might impact annotation, we are using data from pituitary gland (pit henceford) and colon, both from young texel sheep.

Command to align pit data:

    STAR --runThreadN 8 --genomeDir genome/ --readFilesIn fastq/texel_pit_1.fastq fastq/texel_pit_2.fastq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/texel_pit_;

Command to align colon data:

    STAR --runThreadN 8 --genomeDir genome/ --readFilesIn fastq/texel_colon_1.fastq fastq/texel_colon_2.fastq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/texel_colon_;

* -> 3. What files did you create this time? hint: use "ls ./Annotating_the_genome/Bam".



Step 4: Get basic statistics about the alignments
-------------------------------------------------

Flagstat command for pit data:

    samtools flagstat bam/texel_pit_Aligned.sortedByCoord.out.bam;

Flagstat command for colon data:

    samtools flagstat bam/texel_colon_Aligned.sortedByCoord.out.bam;

* -> 4. What percentage of your reads were successfully mapped?
  

Step 5: Index the bam files
----------------------------

In order to view the bam files in the Ensembl browser (Step 7) you will need to index them.

Command to index the pit BAM file:

    samtools index bam/texel_pit_Aligned.sortedByCoord.out.bam;

Command to index the colon BAM file:

    samtools index bam/texel_colon_Aligned.sortedByCoord.out.bam;

Step 6: Assemble the transcripts
--------------------------------

When the reads have been aligned, we can assemble the transcripts. 

Note: scallop writes a lot to STDOUT, so it's best to redirect this output to a file (we add this to the end of the command `> gtf/texel_XXXX_scallop.out`). We also need a folder to store the results, so:

    mkdir -p ./Annotating_the_genome/GTF

Command to assemble transcripts from pit data:

    scallop -i bam/texel_pit_Aligned.sortedByCoord.out.bam -o gtf/texel_pit_scallop.gtf > gtf/texel_pit_scallop.out

Command to assemble transcripts from colon data:

    scallop -i bam/texel_colon_Aligned.sortedByCoord.out.bam -o gtf/texel_colon_scallop.gtf > gtf/texel_colon_scallop.out

* -> 5. What files did you create this time? hint: use "ls ./Annotating_the_genome/GTF".
* -> 6. What is stored in the GTF file?


Step 7: View your results in the Ensembl browser
------------------------------------------------

BAM files are often quite large and are unsuitable for uploading to a website, so in order to view the alignments in the Ensembl browser you need to host the sorted and indexed BAM files on either a webserver or an ftp site. Both the sorted bam file and the index file (ending .bai) are needed to view the alignments on the website. The website requires the bam file URL to be entered, it then looks for a .bai file with the same name in the same directory.

For convenience we have already set up an ftp site that contains all the files you just created.

To view the ftp, enter the following URL into your web browser:

https://ftp.ebi.ac.uk/pub/databases/ensembl/ereboperezsilva/livestock_genomics_2024/
 
We have chosen NT5C2 as a good example of a chr22 gene with differential expression highlighted by the different tissues. 

http://www.ensembl.org/Ovis_aries_rambouillet/Location/View?db=core;g=ENSOARG00020022666;r=22:23420325-23451572

To load the alignments from the Bam files:

* Go to the browser: [www.ensembl.org](http://www.ensembl.org/Ovis_aries_rambouillet/Location/View?db=core;g=ENSOARG00020022666;r=22:23420325-23451572)
* If you are on the Gene view, click on the Location tab at the top left hand side of the page to take you to the view of the region on the chromosome.
* Now we can load our BAM files.
* Click on the “Configure this page” button on the left panel, this opens a configuration panel.
* To load the BAM files click on “Personal data” tab at the top.
* Name the track (if you have an Ensembl account you can store the track in your account to use another time).
* Add the data using the URL of the BAM files (https://ftp.ebi.ac.uk/pub/databases/ensembl/ereboperezsilva/livestock_genomics_2024/bam/)
* Click on the tick in the top right hand corner of the "Configure this page" window to return to the region view browser.
* Once you have attached the remote files you should be able to see them in the region view browser, if they do not show up you may need to turn them on by going to "Configure this page" -> "Your data" and set it to “Normal”.

Done!

* -> 7. What does your BAM track look like compared to the annotated genes? Where do the reads align? hint: look for stacks of reads.
 
To load the transcript models from the GTF files:

* Click on the “Configure this page” button on the left panel, this opens a configuration panel.
* To load the GTF files click on “Personal data” tab at the top.
* Name the track (if you have an Ensembl account you can store the track in your account to use another time).
* Add the data using the URL of the GTF files (https://ftp.ebi.ac.uk/pub/databases/ensembl/ereboperezsilva/livestock_genomics_2024/gft/)
* Click on the tick in the top right hand corner of the "Configure this page" window to return to the region view browser.
* Once you have attached the remote files you should be able to see them in the region view browser, if they do not show up you may need to turn them on by going to "Configure this page" -> "Your data" and set it to “Normal”.

Done.

* -> 8. What does your GTF track look like compared to the Ensembl annotated genes?

Step 9: View a different example
--------------------------------

For the result of slightly more complex (and more time consuming) data, repeat the previous step, but this time use the data in the FTP folder:

https://ftp.ebi.ac.uk/pub/databases/ensembl/ereboperezsilva/bga23/

Load that in a similar fashion, but this time, you MUST use a different genome:
* Go to the browser: www.ensembl.org
* Enter ENSDARG00000055381 into the search box to take you to the gene view page, it is a gene called "bambia", and the genome should be the zebra fish. We want to go on the Location view, so please select "Region in Detail", and...
* Load the data as before. You will observe many diferences.

* -> How do both result compare to each other?
