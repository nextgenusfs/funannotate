#!/usr/bin/env perl 

#This script is lifted from the Maker 2.31.8 package.  Which I believe to have been written by Carson Holt. 

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use warnings;
use Getopt::Long;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "

Synopsis:

map_gff_ids genome.all.id.map genome.all.gff

Description:

This script takes a id map file and changes the name of the
appropriate attributes in a gff file.  The map file is tab delimited
file with two columns: old_name and new_name.The script requires that
old_name (or the part of old_name before the first colon) be equal to
an attribute value before it will map it.

";


my ($help);
my $opt_success = GetOptions('help'    => \$help,
			      );

my ($map_file, $gff_file) = @ARGV;
die $usage unless $map_file && $gff_file;

# Read map file and create a map hash
open (my $MAP, '<', $map_file) or die "Can't open $map_file for reading\n$!\n";
my %map;
map {my ($old, $new) = split;$map{$old} = $new} (<$MAP>);
close $MAP;

# Open $gff_file for reading, unlink it to avoid clobbering it and then open
# the same file for writing.  This allows in-place editing.
open (my $IN, '<', $gff_file) or die "Can't open $gff_file for reading\n$!\n";
unlink($gff_file);
open(my $OUT, '>', $gff_file) or die "Can't open $gff_file for writing\n$!\n";

# Set the order in which attributes will appear in the attribute text.
# Attributes not listed here will be sorted alphabetically.
my %order = (ID     => 1,
	     Parent => 2,
	     Name   => 3,
	     Alias  => 4,
	    );

my $fasta_start;
LINE:
while (my $line = <$IN>) {
	$fasta_start++ if $line =~ /^\#\#FASTA/;
	my ($seq, $source, $type, $start, $end, $score,
	    $strand, $phase, $attrb_text) = split /\t/, $line;
	if (! $fasta_start    && # If we're not in the fasta section...
	    $line !~ /^\s*\#/ && # and we're not a comment line...
	    $start =~ /^\d+$/ && # and we have a valid start column...
	    $end =~ /^\d+$/   && # and a valid end column...
	    $attrb_text          # and attribute text
	   ) {
		chomp $attrb_text;
		my %attrb = parse_attributes($attrb_text);

		# Create an alias with the old_name for gene and mRNA features.
		my $alias;
		if ($type eq 'mRNA' || $type eq 'gene') {
			$alias = $attrb{ID}[0] if $type =~ /^(gene|mRNA)$/;
			#avoids redundancy in existing Aliases
			if($attrb{Alias}){
			    @{$attrb{Alias}} = grep {$_ ne $attrb{ID}[0]} @{$attrb{Alias}};  
			}
			#don't make aliases to Apollo's temporary ID's
			$alias = undef if($alias =~ /^[^\:]+\:temp\d+\:|^CG\:/);
		}

		# Collect all IDs that should be mapped for this feature
		# (the values of Parent and ID tags)
		my @ids;
		push @ids, @{$attrb{Parent}} if $attrb{Parent};
		push @ids, @{$attrb{ID}}     if $attrb{ID};
		push @ids, @{$attrb{Alias}}  if $attrb{Alias};

		# Keep only the IDs that are in the mapping file.
		my @map_ids;
		map {push @map_ids, $_ if $map{$_}}  @ids;

		# If there is nothing to map, print the line and go to
		# the next feature.
		if (! @map_ids) {
			print $OUT $line;
			next LINE;
		};

		# Map old_name to new_name for each tag value...
		for my $tag (keys %attrb) {
		    #avoid replacing values in other tags
		    next unless($tag eq 'ID' || $tag eq 'Name' || $tag eq 'Parent' || $tag eq 'Derives_from');
			for my $value (@{$attrb{$tag}}) {
				for my $id (@map_ids) {
					# Only if the value (or the portion preceding
					# the first colon) is equal to the map key.
					next unless ($value eq $id || $value =~ /^$id:/);
					$value =~ s/$id/$map{$id}/ unless($tag eq 'Name' && $id !~ /\-gene\-\d+\.\d+|^CG\:|^....\:|^[^\:]+\:temp\d+\:/);
				}
			}
		}
		# Now that we're done mapping, add the alias to the attributes
		# if we have one.
		$alias = undef if($alias && ($alias eq $attrb{ID}[0] || $alias eq $attrb{Name}[0]));
		push(@{$attrb{Alias}}, $alias) if $alias;

		# Make sure we have all attribute keys added to sort order to
		# avoid undef in <=>
		map {$order{$_} ||= 99} keys %attrb;

		# Create new attribute text.
		my $new_attrb_text;
		for my $tag (sort {$order{$a} <=> $order{$b} || $a cmp $b} keys %attrb) {
			my $value_text = join ',', @{$attrb{$tag}};
			$new_attrb_text .= "$tag=$value_text;"
		}
		$new_attrb_text .= "\n";

		# Create a new $line.
		$line = join "\t", ($seq,
				 $source,
				 $type,
				 $start,
				 $end,
				 $score,
				 $strand,
				 $phase,
				 $new_attrb_text,
				 );
	}
	print $OUT $line;
}
#-----------------------------------------------------------------------------
#-------------------------- Subroutines --------------------------------------
#-----------------------------------------------------------------------------
sub parse_attributes {
	my $attrb_text = shift;

	my %attrb;
	my @attrb_list = split /\s*;\s*/, $attrb_text;
	for my $attrb_pair (@attrb_list) {
		my ($tag, $value_text) = split /\s*=\s*/, $attrb_pair;
		my @values = split /\s*,\s*/, $value_text;
		$attrb{$tag} = \@values;
	}
	return %attrb;
}
