#!/usr/bin/perl

use strict;
use warnings "all";

if (!defined($ARGV[1])) {
	print STDERR "proteinortho5_clean_edges.pl RM_LIST EDGELIST(S) >CLUSTERED_EDGELIST \nReads Proteinortho outfile and its initial edge list and calculates\nthe remaining edge list after the step of clustering.\n";
	exit;
}


# Proteinortho Output
print STDERR "Cleaning edge list...\n";
my %map;
open(IN,"<$ARGV[0]") || die("Could not open file $ARGV[0]: $!");
while (<IN>) {
	chomp;
	my ($a, $b) = sort split(' ',$_,2);
	unless ($b) {die("Line does not match filter list pattern\n");}
	$map{$a.' '.$b} = 1;
}
close(IN);

# Edgelist
shift(@ARGV);
my $rm = 0;
my $all = 0;
foreach my $file (@ARGV) {
	open(IN,"<$file") || die("Could not open file $file: $!");
	while (<IN>) {
		if ($_ =~ /^#/) {print $_; next;}
		my ($a,$b,undef) = split("\t",$_,3);
		unless ($b) {die("Line does not match Proteinortho graph format\n");}
		$all++;
		($a,$b) = sort($a, $b);

		if (exists($map{$a.' '.$b})) {$rm++; next;}
		print $_;
	}
	close(IN);
}

print STDERR "Removed $rm / $all edges\n";
print STDERR "Done.\n";
