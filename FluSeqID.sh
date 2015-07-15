#!/bin/bash
set -e

# Pipeline to generate the consensus sequence for an unknown influenza virus genome.  Based on Sherlock.sh and adapted to work on genome segments independently.

# Pipeline by ellisrichardj making use of a variety of standard bioinformatic tools

# Prerequisites:
	# bwa - 0.7.5a or above
	# samtools 
	# velvet (must be compiled for multithreading and use of k-mers of at least 101)
	# picard
	# GATK 
	# BLAST+ 
	# IterMap.sh (available at https://github.com/ellisrichardj/csu_scripts)
	# vcf2consensus.pl (available at https://github.com/ellisrichardj/csu_scripts)

# Ensure that each of these are in PATH or symbolic links are in a folder that is in your PATH

# Version 0.2.0	01/06/05	Initial version (FluSearch) adapted from Sherlock version 0.2.0
# Version 0.2.1 22/06/15	Minor edits and new name
# Version 0.2.2 07/08/15	Makes use of pre-formated BLAST db - see FormatFlu_db.sh
# Version 0.2.3 13/07/15	Clarify usage statement; retain original filename for host reference when generating a
#							new index; remove intermediate files from IterMap process
# Version 0.3.0 15/07/15	Major change to embed iterative mapping into this single script

# set defaults for the options
KVALUE=101
CUTOFF=auto
Blast_e_value=0.0001

Start=$(date +%s)

# parse the options
while getopts 'H:e:r:c:k:' opt ; do
  case $opt in
    H) PathToHOST=$OPTARG ;;
    e) Blast_e_value=$OPTARG ;;
    r) PathToSearchData=$OPTARG ;;
    c) CUTOFF=$OPTARG ;;
    k) KVALUE=$OPTARG ;;
  esac
done
# skip over the processed options
shift $((OPTIND-1)) 
# check for mandatory positional parameters
if [ $# -lt 2 ]; then
  echo "
Usage: $0 [options] 
	-H path to HOST genome [required]
	-r path to Reference Database Directory [required]
	Read1.fastq.gz Read2.fastq.gz
"
  echo "Reference database should be a directory containing a separate fasta file for each virus segment
"
  echo "Options: De Novo assembly: -k Velvet kmervalue (default = 101)
		 	-c Velvet cov_cutoff (default = auto)
"
  echo "Options: Blast Output: -e Blast e value (default = 0.0001)
"
  exit 1
fi

LEFT="$(readlink -f "$1")"
RIGHT="$(readlink -f "$2")"

threads=$(grep -c ^processor /proc/cpuinfo)

# Get sample/reference/host names from file names
	sfile1=$(basename "$LEFT")
	sfile2=$(basename "$RIGHT")
	samplename=${sfile1%%_*}

	ref=$(basename "$PathToSearchData")

	host=$(basename "$PathToHOST")
	hostname=${host%%.*}

# Create output OutputDirectory
OutputDir="$samplename"_FluSeqID_"$ref"
mkdir "$OutputDir"
echo "Getting host genome"

# Generate bwa index of host genome if it doesn't already exist
if [ -s "$PathToHOST".sa ] && [ -s "$PathToHOST".amb ] && [ -s "$PathToHOST".ann ] && [ -s "$PathToHOST".bwt ] && [ -s "$PathToHOST".pac ]
then 	echo "Using existing host index"
else 	echo "Generating local index of host genome" 
	ln -s "$(readlink -f "$PathToHOST")" "$OutputDir"/"$hostname".fa
	bwa index "$OutputDir"/"$hostname".fa
	PathToHOST="$OutputDir"/"$hostname".fa
fi

# Map raw data to host genome
echo "Mapping raw reads to host reference genome"
bwa mem -t "$threads" "$PathToHOST" "$LEFT" "$RIGHT" | samtools view -Su - | samtools sort - "$OutputDir"/"$samplename"_"$hostname"_map_sorted

# Extract sequence reads that do not map to the host genome
samtools view -b -f 4 "$OutputDir"/"$samplename"_"$hostname"_map_sorted.bam > "$OutputDir"/"$samplename"_nonHost.bam

# Denovo assembly of non-host reads
velveth "$OutputDir"/"$samplename"_nonHost_"$KVALUE" "$KVALUE" -shortPaired -bam "$OutputDir"/"$samplename"_nonHost.bam
velvetg "$OutputDir"/"$samplename"_nonHost_"$KVALUE" -exp_cov auto -cov_cutoff "$CUTOFF" -ins_length 300 -clean yes -unused_reads yes

# Now run Blast of output contigs against Influenza database of choice.

cd "$OutputDir"
echo "BLASTing contigs against selected database(s)"
for file in "$PathToSearchData"/*.fasta
do
	seg=$(basename "$file")
	segname=${seg%%.*}

	blastall -b 5 -p blastn -d "$file" -i "$samplename"_nonHost_"$KVALUE"/contigs.fa -o "$samplename"_"$segname"_crunch.txt -m 8 -e "$Blast_e_value" &
done
wait

# Extract the single sequence from blast output corresponding to the highest scoring match for each genome segment

for crunch in *.txt
do	
	topseg=${crunch%%.*}
	part=${topseg#*_}
	segname=${part%_*}
	sort -k12,12 -rn "$crunch" | head -1 - | awk '{print $2}' - > "$topseg"_match.txt
#	if [ -s "$topseg"_match.txt ] ; then
		AccNo=$(sed -e 's/.*\([A-Z][A-Z][0-9][0-9][0-9][0-9][0-9][0-9]\).*/\1/g'  "$topseg"_match.txt)
		blastdbcmd -db "$PathToSearchData"/"$segname".fasta -entry_batch "$topseg"_match.txt | sed '1s/.*/>'"$AccNo"-"$segname"'/g' -> "$segname"_match.fas	
#	else shuf 
done

cat *_match.fas > top_matches.fa

# Map original data to selected reference sequences and generate new consensus (iteratively)
#echo "Now mapping to generate new consensus"
#IterMap.sh -i4 top_matches.fas "$(readlink -f "$LEFT")" "$(readlink -f "$RIGHT")"

rfile="$(readlink -f top_matches.fa)"

	ref=$(basename "$Ref")
	refname=${ref%%_*}
	reffile=${ref%%.*}

mkdir "$samplename"_IterMap"$iter"
cd "$samplename"_IterMap"$iter"

iter=4
count=1

if [ $minexpcov -lt 5 ]; then depth=1; else depth=10; fi

threads=$(grep -c ^processor /proc/cpuinfo)

# Log commands that are run
echo -e "$now\n\tItermap v0.3.0 running with $threads cores\n\tThe following commands were run:\n"  > "$samplename"_IterMap"$iter".log
# Define function to log commands and then run them
LogRun(){
echo -e "\n$@" >> "$samplename"_IterMap"$iter".log
eval "$@"
}

	while (($count <= $iter))
	do
	# Set reduced mapping stringency for first iteration
	if [ $count == 1 ] && [ $iter != 1 ]; then
		mem=16
		mmpen=2
		gappen=4
	else # bwa mem defaults
		mem=18
		mmpen=4
		gappen=6
	fi

	# mapping to original reference or most recently generated consensus
	LogRun bwa index "$rfile"
	LogRun samtools faidx "$rfile"
	LogRun picard-tools CreateSequenceDictionary R="$rfile" O=${rfile%%.*}.dict
	LogRun bwa mem -T10 -t "$threads" -k "$mem" -B "$mmpen" -O "$gappen" -R '"@RG\tID:"$samplename"\tSM:"$samplename"\tLB:"$samplename""' "$rfile" "$LEFT" "$RIGHT" |
	 samtools view -Su - |
	 samtools sort - "$samplename"-"$reffile"-iter"$count"_map_sorted
	samtools index "$samplename"-"$reffile"-iter"$count"_map_sorted.bam
	LogRun GenomeAnalysisTK.jar -T RealignerTargetCreator -nt "$threads" -R "$rfile" -I "$samplename"-"$reffile"-iter"$count"_map_sorted.bam -o indel"$count".list
	LogRun GenomeAnalysisTK.jar -T IndelRealigner -R "$rfile" -I "$samplename"-"$reffile"-iter"$count"_map_sorted.bam -targetIntervals indel"$count".list -maxReads 50000 -o "$samplename"-"$reffile"-iter"$count"_realign.bam
	rm "$samplename"-"$reffile"-iter"$count"_map_sorted.bam
	rm "$samplename"-"$reffile"-iter"$count"_map_sorted.bam.bai

if [ $count == $iter ]; then
	# generate and correctly label consensus using cleaned bam on final iteration
# rmdup does not work in v1.x of samtools so replace this with Picard equivalent	LogRun samtools rmdup "$samplename"-"$reffile"-iter"$count"_realign.bam "$samplename"-"$reffile"-iter"$count"_clean.bam
	LogRun picard-tools MarkDuplicates INPUT="$samplename"-"$reffile"-iter"$count"_realign.bam OUTPUT="$samplename"-"$reffile"-iter"$count"_clean.bam REMOVE_DUPLICATES=true METRICS_FILE=dup_metrics"$count".txt
	LogRun samtools view -bF4 -o "$samplename"-"$reffile"-iter"$count"_clean_mapOnly.bam "$samplename"-"$reffile"-iter"$count"_clean.bam
	samtools index "$samplename"-"$reffile"-iter"$count"_clean_mapOnly.bam

	LogRun samtools mpileup -L 10000 -Q1 -AEupf "$rfile" "$samplename"-"$reffile"-iter"$count"_clean_mapOnly.bam |
	 bcftools call -c - > "$samplename"-"$reffile"-iter"$count".vcf
	LogRun vcf2consensus.pl consensus -d "$minexpcov" -Q "$minQ" -f "$rfile" "$samplename"-"$reffile"-iter"$count".vcf |
	 sed '/^>/ s/-iter[0-9]//;/^>/ s/$/'-iter"$count"'/' - > "$samplename"-"$reffile"-iter"$count"_consensus.fa

	# mapping statistics
	LogRun samtools flagstat "$samplename"-"$reffile"-iter"$count"_realign.bam > "$samplename"-"$reffile"-iter"$count"_MappingStats.txt
	rfile="$samplename"-"$reffile"-iter"$count"_consensus.fa

else

	LogRun samtools mpileup -L 10000 -Q1 -AEupf "$rfile" "$samplename"-"$reffile"-iter"$count"_realign.bam |
	 bcftools call -c - > "$samplename"-"$reffile"-iter"$count".vcf
	LogRun vcf2consensus.pl consensus -d "$minexpcov" -Q "$minQ" -f "$rfile" "$samplename"-"$reffile"-iter"$count".vcf |
	 sed '/^>/ s/-iter[0-9]//;/^>/ s/$/'-iter"$count"'/' - > "$samplename"-"$reffile"-iter"$count"_consensus.fa

	# mapping statistics
	LogRun samtools flagstat "$samplename"-"$reffile"-iter"$count"_realign.bam > "$samplename"-"$reffile"-iter"$count"_MappingStats.txt
	rfile="$samplename"-"$reffile"-iter"$count"_consensus.fa

fi

	((count=count+1))
	echo "New Consensus: "$rfile""
done

complete=$(date '+%x %R')
echo -e "Completed processing $complete"  >> "$samplename"_IterMap"$iter".log

# Clean up intermediate files
cd ..
mkdir "$samplename"_FinalResults
mv "$samplename"_IterMap"$iter"/"$samplename"-"$reffile"-iter4_MappingStats.txt  "$samplename"_FinalResults/"$samplename"-"$reffile"-iter4_MappingStats.txt
mv "$samplename"_IterMap"$iter"/"$samplename"-"$reffile"-iter4.vcf  "$samplename"_FinalResults/"$samplename"-"$reffile"-iter4.vcf
mv "$samplename"_IterMap"$iter"/"$samplename"-"$reffile"-iter4_consensus.fa  "$samplename"_FinalResults/"$samplename"-"$reffile"-iter4_consensus.fa
mv "$samplename"_IterMap"$iter"/"$samplename"-"$reffile"-iter3_consensus.fa  "$samplename"_FinalResults/"$samplename"-"$reffile"-iter3_consensus.fa
mv "$samplename"_IterMap"$iter"/"$samplename"-"$reffile"-iter4_clean_mapOnly.bam  "$samplename"_FinalResults/"$samplename"-"$reffile"-iter4_clean_mapOnly.bam
mv "$samplename"_IterMap"$iter"/"$samplename"-"$reffile"-iter4_clean_mapOnly.bam.bai  "$samplename"_FinalResults/"$samplename"-"$reffile"-iter4_clean_mapOnly.bam.bai

rm -rf "$samplename"_IterMap"$iter"
rm "$samplename"_"$hostname"_map_sorted.bam

End=$(date +%s)
TimeTaken=$((End-Start))
echo "Results are in: "$OutputDir"/"$samplename"_FinalResults"
echo  | awk -v D=$TimeTaken '{printf "Performed FluSeqID Analysis in: %02d'h':%02d'm':%02d's'\n",D/(60*60),D%(60*60)/60,D%60}'
