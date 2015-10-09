# FluSeqID
This collection of short scripts forms a pipeline for the detection and extraction of accurate whole genome consensus of Influenza virus from clinical samples, tissue culture or passaged material.  The pipeline has been developed for use with Illumina paired-end data.  

A key prerequisite for this is a properly formatted blast database.  Essentilly, a separate batabase for each of the eight genome segments is required.  There is a script to do this automatically following the download of whole genomes from http://www.fludb.org

In the main script the following steps are run automatically:

1.	Map raw sequence data to host genome (BWA)
2.	Extract reads that do not map to the host (Samtools)
3.	Assemble non host reads (Velvet)
4.	Identify closest match for each genome segment (BLAST)
5.	Map original data to top reference segments (BWA)
6.	Call new consensus (vcf2consensus.pl)
7.	Perform further iterations of steps 5 and 6 to improve new consensus (IterMap)
8.	Output final genome consensus
