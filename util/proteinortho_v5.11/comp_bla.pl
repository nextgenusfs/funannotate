#!/usr/bin/perl

use warnings;
use strict;

unless ($ARGV[1]) {
	print STDERR "Usage: comp_bla.pl FILE_A FILE_B\n\nCompares two Proteinortho-graph files and reports additional and different entrys.\n D = different\n O = only here\n";
	exit;
}

my %a;
open(FILE,"<$ARGV[0]") || die("Error, could not open file $ARGV[0]: $!");
while(<FILE>) {
	chomp;
	if ($_ =~ /^#/) {next;}
	my @row = split(/\s+/);
	my $key = $row[0]." ".$row[1];
	my $val = $_;
	if ($a{$key}) {die("Key $key is defined twice in $ARGV[0]!!!\n")}
	$a{$key} = $val;
}
close(FILE);

my %b;
my %a_b;
open(FILE,"<$ARGV[1]") || die("Error, could not open file $ARGV[1]: $!");
while(<FILE>) {
	chomp;
	if ($_ =~ /^#/) {next;}
	my @row = split(/\s+/);
	my $key = $row[0]." ".$row[1];
	my $val = $_;
	if ($b{$key}) {die("Key $key is defined twice in $ARGV[1]!!!\n")}
	$b{$key} = $val;

	if (exists $a{$key}) {
		if ($a{$key} eq $b{$key}) {
			$a_b{$key} = 1;
		}
		else {
			print "D: $ARGV[0]\t$a{$key}\nD: $ARGV[1]\t$b{$key}\n";
		}
	}
}
close(FILE);

foreach my $key (keys %a) {
	unless (exists $a_b{$key}) {
		print "O: $ARGV[0]\t$a{$key}\n";
	}
}

foreach my $key (keys %b) {
	unless (exists $a_b{$key}) {
		print "O: $ARGV[1]\t$b{$key}\n";
	}
}
