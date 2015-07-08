#!/bin/bash
set -e

# Simple script by ellisrichardj to separate downloaded Influenza genomes one fasta file for each segment
# This will correctly format the database for use in the FluSeqID pipeline

# The multi-segment fasta file can be downloaded from the Influenza Research Database (http://www.fludb.org)
# via the Nucleotide Sequence Search page and selecting 'Complete Genomes Only' radio button
# Squires et al.  Influenza Other Respir Viruses. 2012 Nov;6(6):404-16. http://dx.doi.org/10.1111/j.1750-2659.2011.00331.x

# Prerequisites:
#	BLAST+

# Version 0.0.1	08/07/15	Initial version

OutpurDir="$PWD"

# parse the options
while getopts 'o:' opt ; do
  case $opt in
    o) OutputDir=$OPTARG ;;
  esac
done
# skip over the processed options
shift $((OPTIND-1))

# check for mandatory positional parameters
if [ $# -lt 1 ]; then
  echo "
Usage: $0 <Path to multi-segment Influenza fasta file>
	Options -o Path to Output Directory (default: use current directory)
"
  exit 1
fi

Input="$(readlink -f "$1")"
cd "$OutputDir"

mkdir InfluenzaDB

echo "Splitting database according to genome segment"
awk 'BEGIN {RS=">"} /Segment:1/ {print ">"$0}' $Input > InfluenzaDB/Influenza_seg1.fasta &
awk 'BEGIN {RS=">"} /Segment:2/ {print ">"$0}' $Input > InfluenzaDB/Influenza_seg2.fasta &
awk 'BEGIN {RS=">"} /Segment:3/ {print ">"$0}' $Input > InfluenzaDB/Influenza_seg3.fasta &
awk 'BEGIN {RS=">"} /Segment:4/ {print ">"$0}' $Input > InfluenzaDB/Influenza_seg4.fasta &
awk 'BEGIN {RS=">"} /Segment:5/ {print ">"$0}' $Input > InfluenzaDB/Influenza_seg5.fasta &
awk 'BEGIN {RS=">"} /Segment:6/ {print ">"$0}' $Input > InfluenzaDB/Influenza_seg6.fasta &
awk 'BEGIN {RS=">"} /Segment:7/ {print ">"$0}' $Input > InfluenzaDB/Influenza_seg7.fasta &
awk 'BEGIN {RS=">"} /Segment:8/ {print ">"$0}' $Input > InfluenzaDB/Influenza_seg8.fasta &
wait

for file in InfluenzaDB/*.fasta
	do
		seg=$(basename "$file")
		segname=${seg%%.*}
		makeblastdb -in InfluenzaDB/"$segname".fasta -parse_seqids -dbtype nucl &
	done
wait
exit
