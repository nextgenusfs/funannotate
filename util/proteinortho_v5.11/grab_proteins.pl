#!/usr/bin/perl

## Grabs proteins and tells missing species from groups
use strict;
use warnings;

my $proteins_path = "./";
my $genomes_path = "./";

unless (scalar(@ARGV)) {
	print STDERR "Fetches proteins from orthologous groups\nUsage: grab_proteins.pl [OPTIONS] Proteinortho_Output_Table\nOptions: -path=[PATH2PROTEINS] path to protein files\n (set if full path is not present in header of the Proteinortho output table)\n";
	exit;
}

foreach (@ARGV) {
	if 	($_ =~ /-path=(.+)/) 	{ $proteins_path = $1; }
}

my @paths = ();
my $line = -2;
my %genomes;
my $file = pop(@ARGV);
open(IN,"<$file") || die ("Could not open '$file': $!");
while(<IN>){
	$line++;
	if ($line >= 0) {print STDERR "Writing Set $line... ";}
	chomp;
	my @data = split(/\t+/,$_);
	if ($line == -1) {
		@paths = @data;		# no we know the file for each entry
		my $error = 0;
		print STDERR "Loading data...\n";
		for (my $i = 3; $i < scalar(@paths); $i++) {
			my $genome = $paths[$i];
			$genome =~ s/\....?$/\.fna/;
			if (!-e "$proteins_path/$paths[$i]") 	{$error++; print STDERR "Could not find PROTEINs for '$proteins_path/$paths[$i]'\n";}
			my %hash = fasta2hash("$proteins_path/$paths[$i]",1);
			$genomes{$genome} = \%hash;
		}
		if ($error > 0) {die("Could not proceed. $error errors occured!\n");}
	}
	if ($data[0] =~ /^#/) {next;}
	
	# Fetch all files of
	my $genes = 0;
	open(OUT,">$line.fasta");
	for (my $i = 3; $i < scalar(@data); $i++) {
		$data[$i] =~ s/,.*//;
		my $name = $paths[$i];
		$name =~ s/\....?$//;
		$name =~ s/\s//g;
		if (length($name) < 1) {die("Error at line $line, entry $i\n");}
#		print " $name";
		my $genome = $name;		
		$genome .= ".fna";
		if ($data[$i] eq "_GENOME_HIT_") {

		}
		elsif ($data[$i] eq "*") {
			# SET GENOME

#			print " -> *\n";
		}
		else {
		# open proteome
		if (!defined($genomes{$genome}{$data[$i]})) {die("Could not find '$data[$i]' in '$proteins_path/$paths[$i]'\n ");}
		print OUT ">$data[$i] $name\n$genomes{$genome}{$data[$i]}\n";
		$genes++;
		}
	}
	close(OUT);

	print STDERR "$line.fasta ($genes genes)\n";
}
close(IN);



# Reads a fasta file, returns two array pointers for header an sequence
# Access e.g. via $results[x]->[y]
#my %hash = fasta2hash($file);
sub fasta2hash {
	my $file = shift;
	my $cut = shift;
	if (!defined($cut)) {$cut = 0;}
	my %hash;
	open(FILE,"<$file") || die("fasta2hash(): Could not open file: $file\n");
	my $last_hash = "";
	while(<FILE>) {
		$_ =~ s/\r/\n/g;
		chomp;
		if ($_ =~ /^>/) {
			$_ =~ s/^>//;
			if ($cut) {
				$_ =~ s/ .*$//s;
			}
			$hash{$_} = "";
			$last_hash = $_;
			next;
		}
		$hash{$last_hash} .= $_;
	}
	close(FILE);
	return %hash;
}


