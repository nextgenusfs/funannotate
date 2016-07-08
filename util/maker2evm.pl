#!/usr/bin/env perl -w
use strict;
use warnings;

#written by Jason Stajich
#https://github.com/hyphaltip/genome-scripts/blob/master/gene_prediction/maker2evm.pl
#modified by Jon Palmer
#while Jason originally wrote this to pull out predictions from the de novo predictors, I'm changing to pulling out the models predicted by maker and then running through EVM.  This might not change any of the results, but for funannotate that is kind of the point.  Could maybe make this script more flexible to do a few different things.

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
    if( $row[1] =~ /maker/ ) {
	if( $row[2] eq 'gene' ) {
	} elsif( $row[2] eq 'mRNA' ) {
	} elsif( $row[2] eq 'exon' ) {
	} elsif( $row[2] eq 'CDS' ) {
	} elsif( $row[2] eq 'tRNA' ) {
	    print $gene join("\t",@row),"\n";
	} else {
	    warn("unknown field type for $row[2]\n");
	}
	print $gene join("\t", @row), "\n";
    } elsif( $row[1] eq 'est2genome' ) {
	next if $row[2] eq 'expressed_sequence_match';
	if( $row[2] eq 'match_part' ) {
	    $row[2] = 'EST_match';
	}
	} elsif( $row[1] eq 'cdna2genome' ) {
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
    } elsif( $row[1] =~ /augustus_masked|genemark|snap_masked|repeatmasker|blast[xn]|tblastx|repeatrunner|\./) {
	next; # skipping these
    } else {
	warn("unknown type for $row[1] $row[2]\n");
	next;
    }
#    print join("\t", @row),"\n";
}