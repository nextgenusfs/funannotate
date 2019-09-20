#!/usr/bin/env perl

#script posted online here: https://groups.google.com/forum/#!topic/blast2go/ZX6565FqSZw via a user named Christoph

use strict;
use warnings;
use Getopt::Long;

my $total_xml = $ARGV[0];
my $destination = $ARGV[1];
my $help;
my $switch = 0;
my (@array, @title, @small_title);

my $USAGE = "\nThis script splits the xml output of Interproscan 5 into individually numbered xml files for each protein
		USAGE: ./prepare_ind_xml.pl *.xml <save to>\n\n";

GetOptions (	"help!" => \$help) or die "Incorrect usage!\n$USAGE";


if ((!$ARGV[0]) || (!$ARGV[1]) || ($help)){
	print "$USAGE\n";
	exit;
}

open (XML_TOTAL,"<$ARGV[0]") or die $!;
while (my $line = <XML_TOTAL>){
	if (($line =~ /<protein>/) && ($switch == 0)){
		$switch = 1;
		push (@array, "<EBIInterProScanResults><interpro_matches>");
		push (@array, $line);
	}elsif (($switch == 1) && !($line =~ /<\/protein>/)){
		push (@array, $line);
		if ($line =~ /<xref id="(.+?)"/){
			@title = split (/\"/, $line);
			@small_title = split (/ /, $title[1]);
			print $small_title[0] . "\n";
			open (XML,">$destination/$small_title[0].xml") or die $!;
		}
	}elsif (($line =~ /<\/protein>/) && ($switch == 1)){
		push (@array, $line);
		push (@array, "</interpro_matches></EBIInterProScanResults>");
		$switch = 0;
		for (@array){
			chomp;
			print XML $_ . "\n";
		}
		close XML;
		undef @array;
	}	

}

close XML_TOTAL;
