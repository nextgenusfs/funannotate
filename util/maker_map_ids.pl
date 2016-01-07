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

maker_map_ids --prefix PYU1_ --justify 8 genome.all.gff > genome.all.id.map

Description:

This script wil take a GFF3 file and create a mapping file of gene and
transcript IDs to more numerically incremented human friendly unique
IDs.

Options:

  --prefix      The prefix to use for all IDs (default = 'MAKER_')
  --suffix      A suffix to use for all transcript IDs (default = '-R').
                Specifying --suffix will also turn on --iterate.
  --initial     The initial value to start with for ID count. If --initial
                is supplied more than once than the first value will be used
                for genes and the second for transcripts (default = 1)
  --abrv_gene   Optional abreviation added to IDs for genes (i.e. = 'G')
  --abrv_tran   Optioanl abreviation added to IDs for transcripts (i.e. = 'T')
  --iterate     Transcript iteration to add to IDs. Value can be '0', '1',
                or 'A', i.e. mRNA-0 or mRNA-1 or mRNA-A (default = 'A')
  --justify     The unique integer portion of the ID will be right justified
                with '0's to this length (default = 8)
  --sort_order  A tab delimited file containing two columns: contig_id
                and sort_order.  Genes and transcripts will be named
                in consecutive order along the contigs, and the
                contigs will be sorted in the order specified by the
                file.  If sort_order is not given and there are
                ##sequence-region directives at the top of the gff
                file then naming will be ordered by decreasing contig
                length.

Formating:
  By default gene and transcript abreviations will appear at the end of the
  suffix and the iterate value will appear at the end of the prefix (only on
  transcripts).  You can specify an alternate location for the abbreviation by
  placing a '?' character as part of the prefix or the sufix. An alternate
  location for the iterator can be given by adding a % character to the
  suffix or prefix.

";


my $prefix = 'MAKER_';
my $suffix = '-R%';
my $abrv_gene = '';
my $abrv_tran = '';
my $iterate = 'A';
my $justify = 8;
my $sort_order;
my $apollo;
my @initial;

my $opt_success = GetOptions('prefix=s'     => \$prefix,
			     'suffix=s'     => \$suffix,
			     'abrv_gene=s'  => \$abrv_gene,
			     'abrv_tran=s'  => \$abrv_tran,
			     'iterate=s'    => \$iterate,
			     'initial=s'    => \@initial,
			     'justify=s'    => \$justify,
			     'sort_order=s' => \$sort_order,
			     'apollo'       => \$apollo,
			     'help'         => sub{print $usage; exit(0);},
			      );

my $gff_file = shift;

if (! $opt_success || ! $gff_file){
    print $usage;
    exit(0);
}

if($suffix ne '' && $iterate eq ''){
    $iterate = 'A';
}

if($iterate ne '' && $iterate !~ /^0$|^1$|^A$/i){
    die "ERROR: Invalid value for iterate: $iterate\n";
}

$prefix .= '?' if($prefix !~ /\?/ && $suffix !~ /\?/);
$suffix .= '%' if($suffix !~ /\%/ && $prefix !~ /\%/);
$iterate = uc($iterate);
$justify ||= 8;
@initial = (1) if(! @initial);

my $sort_map = parse_sort_order($sort_order, $gff_file);

#initialize counter
my %counter;
$counter{gene} = $initial[0];
$counter{mRNA} = (@initial > 1) ? $initial[1] : $initial[0];

#teststring to see if current ID already matches this format
my $test_g_id = $prefix.'(\d+)';
$test_g_id =~ s/\?/$abrv_gene/g; #add abrevition
$test_g_id =~ s/\%//g; #add iterator

my $test_t_id = $prefix.'(\d+)'.$suffix;
$test_t_id =~ s/\?/$abrv_tran/g; #add abrevition
$test_t_id =~ s/\%/\(\[A\-Z\]+\)/g if($iterate eq 'A'); #add iterator
$test_t_id =~ s/\%/\(\[0\-9\]+\)/g if($iterate eq '0'); #add iterator
$test_t_id =~ s/\%/\(\[1\-9\]\[0\-9\]*\)/g if($iterate eq '1'); #add iterator
$test_t_id =~ s/\%//g if($iterate eq ''); #add iterator

# Build a hash of ids that need to be mapped;
my %ids;
my %parents;
my %exists;
my %gene_index;
my %mRNA_index;

open (my $IN, '<', $gff_file) or die "ERROR: Can't open $gff_file for reading\n$!\n";
while (<$IN>) {
    last if /^\#\#FASTA|^>/;
    next if /^\s*\#/;
    
    my ($seq, $source, $type, $start, $end, $score, $strand,
	$phase, $attrb_text) = split /\t/, $_;
    
    # Here we explicity limit the mapping to genes and mRNAs
    if ($type && ($type =~ /^(gene|pseudogene)$/ || $type =~ /(RNA|transcript)$/)) {
	$type = 'gene' if($type eq 'pseudogene');
	$type = 'mRNA' if($type =~ /RNA|transcript/);

	my ($id) = $attrb_text =~ /ID=([^\;\n]*)/;
	my ($name) = $attrb_text =~ /Name=([^\;\n]*)/;
	my ($parent) = $attrb_text =~ /Parent=([^\;\n]*)/;	
	
	if($type eq 'gene'){
	    $ids{$seq}{$id}{ID}     = $id;
	    $ids{$seq}{$id}{name}   = $name;
	    $ids{$seq}{$id}{type}   = $type;
	    $ids{$seq}{$id}{start}  = $start;
	    $ids{$seq}{$id}{end}    = $end;
	    $ids{$seq}{$id}{strand} = ($strand eq '+') ? 1 : -1;
	    
	    if($id =~ /^$test_g_id$/ || ($name && $name =~ /^$test_g_id$/)){
		$gene_index{$id}{count} = $1;
		$counter{gene} = $1 + 1 if($1 >= $counter{gene}); #make one above (next usable value)
	    }
	}
	elsif($type eq 'mRNA'){
	    my %mRNA = (ID     => $id,
			name   => $name,
			type   => $type,
			start  => $start,
			end    => $end,
			parent => $parent,
			strand => ($strand eq '+') ? 1 : -1
			);

	    $ids{$seq}{$parent}{mRNAs}{$id} = \%mRNA;
	    
	    if($id =~ /^$test_t_id$/ || ($name && $name =~ /^$test_t_id$/)){
		if($iterate ne ''){
		    if(!$gene_index{$parent}{iterate} ||
		       ($iterate eq 'A' && $2 ge $gene_index{$parent}{iterate}) ||
		       ($iterate ne 'A' && $2 >= $gene_index{$parent}{iterate})
		       ){
			$gene_index{$parent}{iterate} = $2;
			$gene_index{$parent}{iterate}++; #make one above (next usable value)
		    }

		    $mRNA_index{$id}{iterate} = $2;
		}

		$counter{mRNA} = $1 + 1 if($1 >= $counter{mRNA}); #make one above (next usable value)
		$mRNA_index{$id}{count} = $1;
	    }
	}
	
	$exists{$id}++;
	$exists{$name}++ if($name);
    }
}

# Create the new ID map.  Sort contigs by $sort_map in outer loop and 
# features by start in second loop.
foreach my $contig_id (sort {sort_contigs($a, $b, $sort_map)} keys %ids) {
    my $contig = $ids{$contig_id};
    #first build IDs for genes
    foreach my $gene (sort {$a->{start} <=> $b->{start} || $b->{strand} <=> $a->{strand} || $a->{end} <=> $b->{end}} values %{$contig}) {
	my $gid = $gene->{ID};
	my $gcount;
	my $new_id;

	#create new ID
	while(! $new_id){
	    if(! $gene_index{$gid}){
		$gene_index{$gid}{count} = sprintf '%0'.$justify.'s', $counter{gene}++;
	    }
	    if($iterate ne '' && ! defined $gene_index{$gid}{iterate}){
		$gene_index{$gid}{iterate} = $iterate;
	    }

	    my $gcount = $gene_index{$gid}{count};

	    #create new ID
	    $new_id = $prefix.$gcount;
	    $new_id =~ s/\?/$abrv_gene/g; #add abrevition
	    $new_id =~ s/\%//g; #add iterator empty for genes
	}

	#print gene_id
	print "$gid\t$new_id\n";

	#now create IDs for mRNAs of the gene
	foreach my $mRNA (sort {$a->{start} <=> $b->{start} || $a->{end} <=> $b->{end}} values %{$gene->{mRNAs}}){
	    my $tid = $mRNA->{ID};
	    my $parent = $mRNA->{parent};
	    my $tcount;
	    my $iter;
	    $new_id = '';

	    while(! $new_id){
		if($mRNA_index{$tid} && $iterate ne '' && $mRNA_index{$tid}{count} == $gene_index{$parent}{count}){
		    $tcount = $mRNA_index{$tid}{count};
		    $iter = $mRNA_index{$tid}{iterate};
		}
		elsif($mRNA_index{$tid} && $iterate eq ''){
		    $tcount = $mRNA_index{$tid}{count};
		    $iter = '';
		}
		elsif($iterate ne ''){
		    $tcount = $gene_index{$parent}{count};
		    $iter = $gene_index{$parent}{iterate}++;
		}
		else{
		    $tcount = sprintf '%0'.$justify.'s', $counter{mRNA}++;
		    $iter = '';
		}

		#create new ID
		$new_id = $prefix.$tcount.$suffix;
		$new_id =~ s/\?/$abrv_tran/g; #add abrevition
		$new_id =~ s/\%/$iter/g; #add iterator

	    }

	    #print gene_id
	    print "$tid\t$new_id\n";
	}
    }
}
#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------
sub parse_sort_order {
    my ($sort_order_file, $gff_file) = @_;
    
    my %sort_order;
    
    if ($sort_order_file) {
	open(my $IN, '<', $sort_order_file) or die "Can't open $sort_order_file\n$!\n";
	%sort_order = map {chomp $_; split(/\t/, $_)} <$IN>;
	close $IN;
    }
    else{
	my $flag;
	open(my $IN, '<', $gff_file) or die "Can't open $gff_file\n$!\n";
	while (my $line = <$IN>) {
	    if( $line =~ /\#\#sequence-region|\tcontig\t/){
		$flag++;
		last;
	    }
	}
	close($IN);
	    
	if ($flag){
	    open(my $IN, '<', $gff_file) or die "Can't open $gff_file\n$!\n";
	    while (my $line = <$IN>) {
		if ($line =~ /\#\#sequence-region\s+(\S+)\s+(\d+)\s+(\d+)/ || $line =~ /^(\S+)\t\S+\tcontig\t(\d+)\t(\d+)/) {
		    my ($contig_id, $start, $end) = ($1, $2, $3);
		    $sort_order{$contig_id} = ($end - $start);
		}
	    }
	    close $IN;
	}
	else {
	    open(my $IN, '<', $gff_file) or die "Can't open $gff_file\n$!\n";
	    while (<$IN>) {
		my @fields = split /\t/, $_;
		next unless @fields == 9;
		my $seq = $fields[0];
		$sort_order{$seq}++;
	    }
	    close($IN);
	}
    }

    return \%sort_order;
}
#-----------------------------------------------------------------------------
sub sort_contigs {
	my ($a, $b, $sort_map) = @_;
	return ($sort_map->{$a} <=> $sort_map->{$b}) if defined $sort_order;
	return ($sort_map->{$b} <=> $sort_map->{$a});
}
