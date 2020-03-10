#!/bin/bash
for file in Downsampling/BAM_input/*; do
	out="Downsampling/BAM_downsampled/"
	in="Downsampling/BAM_input/"
	echo "processing.." ${in}"$(basename "$file")"
	bbmap/reformat.sh in=${in}"$(basename "$file")" out=${out}"$(basename "$file")".downsampled.bam samplereads=350000
    echo "processed." ${out}"$(basename "$file")".downsampled.bam
    
done

Rscript bpr_downsamp.R