# AIMS-tag
Automated workflow for de novo transcriptome assembly and Tag-Seq gene expression analysis
AIMS-tag is a set of scripts to be run to assemble a de novo transcriptome from paired end read RNA-seq data (scripts 1-10).  The assembled transcriptome is then used as a reference to blast Tag-seq single end reads from multiple samples to provide the data to run gene expression analysis using R.

01_prepare.job

This script takes uses fastq format RNA-seq data to assemble the de novo transcriptome.  If you have multiple files for the one sample it is best to concatenate them into one file for left and right.  It firstly trims the reads removing the tags used for sequencing using the fastx clipper, then runs a quality filter using fastq quality filter.  This step outputs files that correspond to the original file name with a ‘.trim’ added to the file name and stored in 01_results folder.  A count of the lines in each sample file is also undertaken and stored in 01_counts.trim file in the ‘counts’ folder.

All of the left and right reads from the RNA-seq data that has undergone trimming and filtering are then concatenated into R1.trim.cat and R2.trim.cat stored in the 01_results folder.

These files are then repaired and will be named R1_R1.trim.cat, R2_R2.trim.cat and Unp_R1.trim.cat_R2.trim.cat stored in the 01_results folder.  These files then will be counted and the results stored in the counts folder under 02_counts.repair.dedup

These files will then undergo de duplication to remove any duplicate reads using the deupTranscriptome.pl script and the lines will be again counted and added to the 01_counts.repair.dedup file.  The deduplicated files are named R1_R1.trim.cat.dedup, R2_R2.trim.cat.dedup and Unp_R1.trim.cat_R2.trim.cat.dedup.

To prepare the files for assembly a suffix needs to be added to the header of each read in the file with /1 for left and unpaired, /2 for right reads.  This step concatenates the R1 and Unp files into R1p_suf1_R1 and the R2 files into R2_suf2.R2.

02_assembly.job

This script invokes the Trinity package to undertake assembly of the output files from 01_prepare.job (R1p_suf1_R1 and R2_suf2.R2) using the normalize reads option.  This output file is stored as Trinity.fasta in the 02_results_trinity folder.

The stats on the assembly are undertaken using the seq_stats.pl script and output into the counts folder in the file 02_seq_stats_assm

Sequences smaller than 400bp are removed using the removesmalls.pl script and output to the file 02_Trinity-l400.fasta in 02_results_trinity and the stats on this file are then output into the 02_seq_stats_assm_removesmal file in the counts folder.

The GC content on the 02_Trinity-l400.fasta file are undertaken and output into the 02_GC_bbmap_Trinisty-l400 file in the counts folder using bbmap stats.sh.

03a_uniprot.job

this script fistly makes a unprot databade to blast the tanscriptome against.

It then takes the 02_Trinity-l400.fasta file from the previous step and splits it into chunks (the size of the chunks can be changed in the configfile.txt by changing the number after CHUNKSIZE.  The chunk size is currently set to 200bp.  These chunks are stored in the 03_results folder.
This script then submits a job for the 3b_uniprot_blast.job which blasts each chunk against the data base. The output files are then stored in 03_results under All.uni.br

04_isogroup.job
this script then creates the files where the sequences are assigned to their isogroups.  
The 04_seq2iso.tab is used for analysis in R and 04_ALL_iso_corrected.fasta is used in the subsequent steps as it is the complete transcriptome.

05_getgenenames.job
this script uses the getGeneNameFromUniProtKB.pl script to find the gene names associated with the isogroups from the UniProt database. It uses the All.uni.br and the 04_ALL_iso_corrected.fasta as input files.

06_contiguity.job
this script calculates the contiguity of the transcriptome using  the 04_ALL_iso_corrected_hits.tab input file.

07_GO.job
this script extracts the Gene Ontology terms associated with the reads.  It uses the All.uni.br and the 04_ALL_iso_corrected.fasta as input files.

08a_keg.job
This is the first step to the KEG annotations to the transcriptome using the 04_ALL_iso_corrected.fasta input file. The output from this step need to be copied to your computer and then:
####use web browser to submit transcriptome_4kegg.fasta file to KEGG's KAAS server 
####( http://www.genome.jp/kegg/kaas/ )

####select SBH algorithm 
####upload nucleotide query (transcriptome_4kegg.fasta)
####select representative genes for model inverts
#### copy this into the	representetive genes for inverts: dme, cel, nve, hmg, aqu, cin
####(Drosophila, C.elegans, N.vec, Hydra, A.queenslandica, Ascidian)

##### Once it is done, download the 'text' output from KAAS, name it query.ko (default)
#####put this file back into your 08_results folder"

08b_keg.job
is the second processes for KEG annotation.  Once you have downloaded the query.ko file and put it back into the results folder then run the script.

09_kog_cegma.job
this script creates a prots database downloaded from http://korflab.ucdavis.edu/Datasets/genome_completeness/core/248.prots.fa.gz
then blasts the transcriptome against this dataset and stores the results in 09_results.
The percentage of KOG terms represented is then calculated.  It uses the 04_ALL_iso_corrected.fasta as an input file.

10a_kog .job
this script uses the All.uni.br and the 04_ALL_iso_corrected.fasta as input files.  The output file *.PRO.fas needs to be downloaded to your computer and then upload the file onto the website http://weizhong-lab.ucsd.edu/webMGA/server/kog/ with default e value = 0.001.  once downloaded go onto script 10b.

10b_kog.job
In this script you need to change the ‘JOBID=’ number to the number for your job from the website in step 10a.  then run the script which will download your results, unzip them and store them in the 10_results folder.


20_read.counts.job
This file reads in the Tag-seq fastq files.  It is best to concatenate all the reads for the one sample if there are multiple files for each sample.  This file then uses fastx_clipper and quality filter to remove the tags and filter out poor quality reads.  The output files are stored in 20_results.  counts of these files are then undertaken and stored in 20_tagseq.counts.trim in the counts folder.

21_mapping.job
this script creates the transcriptome for blasting using bowtie2 and uses the input file 04_ALL_iso_corrected.fasta and renames it to transcriptome.fasta.  this script then reads in the quality filtered and trimmed fastq files for each sample from the 20_results folder and blasts them against the transcriptome, and outputs files ending in .sam for each sample stored in the 21_results folder.  This step also has two output files called mapping.e### and mappint .o### these are important files that contain the percentage of allignement for the samples.

22_sc.job
this script then uses the .sam files and the samcount.pl script to calculate the representation of isogroups in the transcriptome from each sample using the 04_seq2iso.tab file.  The output from this step is stored in 22_results as .counts files.  The expressions_compiler.pl script is then used on all the .counts files in 22_results and combined into an allcounts.txt file.  Some changes to the allcounts.txt files are then made and the outputs of allcontsTankExpt.txt and allcountsSSCExpt.txt are made.


 


