#!/usr/bin/perl
use strict;
use warnings "all";

if (!defined($ARGV[0]) || $ARGV[0] eq "-h" || $ARGV[0] =~ /^-?-help/) {
	print STDERR "proteinortho5_singletons.pl FASTA1 FASTA2 FASTAN <PROTEINORTHO_OUTFILE >SINGLETON_GENES\n";
	print STDERR "Reads Proteinortho outfile and its source fasta files to determin entries which occure once only\n\n";
	exit;
}

my %present; # Genes present in the matrix

### Parse matrix, store present genes and species order
my %order;
while (<STDIN>) {
	if ($_ =~ /#\sSpecies/) {
		my @species = split(/\s+/);
		shift @species;shift @species;shift @species;shift @species;
		for (my $i = 0; $i < scalar(@species); $i++) {
			$order{$species[$i]} = $i;
		}
		next;
	}
	if ($_ =~ /^#/ || length($_) < 4) {next;}
	chomp;
	my @row = split(/\s|,/,$_);
	for (my $i = 3; $i < scalar(@row); $i++) {
		$present{$row[$i]} = 1;
	}
}

### For each fasta file, parse gene-IDs
foreach my $file (@ARGV) {
	# print "\# $file\n";
	my $pos;
	unless (defined($order{$file})) {
		my $short_version = $file;
		$short_version =~ s/^.*\///;
		unless (defined($order{$short_version})) {
			die("Species $file is not in the matrix\n");
		}
		$pos = $order{$short_version};		
	}
	else {$pos = $order{$file};}
	open(FILE,"<$file") || die("Error, could not open file $file: $!");
	while(<FILE>) {
		if ($_ =~ /^>([^\s]+)/) {
			my $id = $1;
			if (exists($present{$id})) {next;}
			
			# add to matrix
			print "1\t1\t0";
			for (my $i = 0; $i < scalar(keys %order); $i++) {
				print "\t";
				if ($i == $pos) {
					print $id;
				}
				else {
					print "*";
				}
			}
			print "\n";
		}
	}
	close(FILE);
}

