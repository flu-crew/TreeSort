#!/bin/bash

# Usage: ./prepare_dataset.sh [--segments "..." --fast] fasta_path reference_segment outdir
# Using --fast will make all trees to be inferred with FastTree.
# By default (without --fast) the reference tree is inferred with IQ-Tree, which is recommended for better accuracy.
# Example usage: ./prepare_dataset.sh --segments "HA,NA" segments.fasta HA myoutdir
# Example with default segments:  ./prepare_dataset.sh segments.fasta HA myoutdir

# These are the default segment names
declare -a segments=("PB2" "PB1" "PA" "HA" "NP" "NA" "MP" "NS")
FAST=0

POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
	case $1 in
		--segments)
			SEGMENTS_STR="$2"
			segments=(${SEGMENTS_STR//,/ })
			shift  # past argument
			shift  # past value
			;;
		--fast)
			FAST=1
			shift
			;;
		-*|--*)
			echo "Unrecognized option $1"
			exit 1
			;;
		*)
			POSITIONAL_ARGS+=("$1")  # save positional arg
			shift  # past argument
			;;
	esac
done

set -- "${POSITIONAL_ARGS[@]}"

# Required arguments:
main_fasta="$1"  # Provide a path to a fasta file with all segments
ref_seg="$2"  # Name of the segment to use as the reference (typically - HA)
outdir="$3"  # Path to the directory to store the results

rm -r $outdir  # Clear out the directory
mkdir $outdir  # Re-create the directory

name=${main_fasta##*/}

# Split out the segments and align them
for seg in "${segments[@]}"
do
	cat $main_fasta | smof grep "|${seg}|" > ${outdir}/${seg}-${name}
	echo "Aligning ${seg}..."
	mafft --thread 6 ${outdir}/${seg}-${name} | sed "s/|${seg}|/|/g"> ${outdir}/${seg}-${name}.aln
	rm ${outdir}/${seg}-${name}
done

if [ $FAST -eq 0 ]; then
	# Build fasttree trees in parallel for non-reference segments
	echo "Building non-reference trees in parallel with FastTree..."
	for seg in "${segments[@]}"
	do
		if [ $seg != $ref_seg ]; then
			fasttree -nt -gtr -gamma ${outdir}/${seg}-${name}.aln > ${outdir}/${seg}-${name}.tre &
		fi
	done
	wait  # Wait to finish.

	# Build an IQ-Tree tree for the reference segment. We use the GTR+F+R5 model by default which can be changed
	echo "Building the reference tree with IQ-Tree..."
	iqtree2 -s ${outdir}/${ref_seg}-${name}.aln -T 6 --prefix "${outdir}/${ref_seg}-${name}" -m GTR+F+R5
	mv ${outdir}/${ref_seg}-${name}.treefile ${outdir}/${ref_seg}-${name}.tre
else
	# Build all trees with FastTree in parallel.
	echo "Building trees in parallel with FastTree..."
	for seg in "${segments[@]}"
	do
		fasttree -nt -gtr -gamma ${outdir}/${seg}-${name}.aln > ${outdir}/${seg}-${name}.tre &
	done
	wait  # Wait to finish.
fi

# Root the trees with a custom rooting script (in parallel)
echo "Rooting trees with TreeTime..."
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
echo "The descriptor file was written to ${descriptor}"
