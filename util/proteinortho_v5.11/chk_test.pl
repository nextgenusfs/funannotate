#!/usr/bin/perl

use warnings;
use strict;

open(FILE,"<$ARGV[0]") || die("Error, could not open file $ARGV[0]: $!");
my @in = <FILE>;
close(FILE);

if (scalar(@in) < 30) {
	print STDERR "Something went wrong. Test failed!\n";
	exit 1;
}
