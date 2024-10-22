#!/bin/bash

# NOTE: You can redefine this list based on the segments that you need for your analysis
declare -a segments=("PB2" "PB1" "PA" "HA" "NP" "NA" "MP" "NS")

# Required arguments:
main_fasta="$1"  # Provide a path to a file with all segments
ref_seg="$2"  # Name of the segment to use as the reference (typically - HA)
outdir="$3"  # Path to the directory to store the results

rm -r $outdir  # Clear out the directory
mkdir $outdir  # Re-create the directory

name=${main_fasta##*/}

# Split out the segments and align them
for seg in "${segments[@]}"
do
	cat $main_fasta | smof grep "|${seg}|" > ${outdir}/${seg}-${name}
	echo "Aligning ${seg}"
	mafft --thread 6 ${outdir}/${seg}-${name} | sed "s/|${seg}|/|/g"> ${outdir}/${seg}-${name}.aln
	rm ${outdir}/${seg}-${name}
done

# Build fasttree trees in parallel
echo "Building trees in parallel"
for seg in "${segments[@]}"
do
	fasttree -nt -gtr -gamma ${outdir}/${seg}-${name}.aln > ${outdir}/${seg}-${name}.tre &
done
wait  # Wait to finish.

# Root the trees with a custom rooting script (in parallel)
for seg in "${segments[@]}"
do
	python treetime-root.py ${outdir}/${seg}-${name}.tre ${outdir}/${seg}-${name}.aln &
done
wait

# Create a descriptor file
descriptor=${outdir}/descriptor.csv
for seg in "${segments[@]}"
do
	if [ $seg == $ref_seg ]; then
		echo -n "*" >> $descriptor
	fi
	echo "${seg},${seg}-${name}.aln,${seg}-${name}.aln.rooted.tre" >> $descriptor
done

