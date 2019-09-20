#!/bin/bash

#simple wrapper for running aligner program and piping output to samtools view/sort

if [ -z "$3" ]; then
    echo 'Usage: sam2bam.sh "aligner_command" bam_threads bam_output'
    echo '**The double quotes are required around aligner command**'
    exit
fi

#construct the command
cmd="$1 | samtools view -@ $2 -bS - | samtools sort -@ $2 -o $3 -"

#run the command
eval $cmd
