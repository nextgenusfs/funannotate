#!/bin/bash

#simple wrapper for running aligner program and piping output to samtools view/sort

if [ -z "$4" ]; then
    echo 'Usage: sam2bam.sh "aligner_command" assembly bam_threads bam_output'
    echo '**The double quotes are required around aligner command**'
    exit
fi

#construct the command
cmd="$1 | samtools view -T $2 -@ $3 -bS - | samtools sort --reference $2 -@ $3 -o $4 -"

echo $cmd

#run the command
eval $cmd
