#!/usr/bin/env perl -w

#This script is written by Jason Stajich, lifted from https://github.com/hyphaltip/genome-scripts/blob/master/gene_prediction/maker2evm.pl
use strict;
use warnings;

open(my $tr => ">transcript_alignments.gff3")|| die $!;
open(my $gene => ">gene_predictions.gff3")|| die $!;
open(my $pep => ">protein_alignments.gff3")|| die $!;

while(<>) {
    if (/^\#/) {
	next;
    } elsif (/^>/ ) {
	last;
    }
    chomp;
    my @row = split(/\t/,$_);
    if( $row[1] =~ /genemark|snap_masked|augustus_masked/ ) {
	$row[1] =~ s/_masked//;
	if( $row[2] eq 'match' ) {	    
	    $row[2] = 'gene';
	    print $gene join("\t",@row),"\n";
	    $row[2] = 'mRNA';
	} elsif( $row[2] eq 'match_part' ) {
	    $row[2] = 'CDS';
	} else {
	    warn("unknown field type for $row[2]\n");
	}
	print $gene join("\t", @row), "\n";
    } elsif( $row[1] eq 'est2genome' ) {
	next if $row[2] eq 'expressed_sequence_match';
	if( $row[2] eq 'match_part' ) {
	    $row[2] = 'EST_match';
	}
	print $tr join("\t", @row), "\n";
    } elsif( $row[1] eq 'protein2genome' ) {
	next if $row[2] eq 'protein_match';
	if( $row[2] eq 'match_part' ) {
	    $row[2] = 'nucleotide_to_protein_match';
	}
	print $pep join("\t", @row), "\n";
    } elsif( $row[1] =~ /maker|repeatmasker|blast[xn]|repeatrunner|\./) {
	next; # skipping these
    } else {
	warn("unknown type for $row[1] $row[2]\n");
	next;
    }
#    print join("\t", @row),"\n";
}