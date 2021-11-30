#!/usr/bin/env perl 

#script from the Maker 2.31.8 distribution

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use warnings;
use strict;

#usage statement
my $usage = "
USAGE:
      genemark_gtf2gff3 <filename>

      This converts genemark's GTF output into GFF3 format.
      The script prints to STDOUT. Use the '>' character to
      redirect output into a file.

";

my $file = shift;

#error checking
if (! $file ){
    print $usage;
    exit;
}

if(! -e $file){
    warn "ERROR: The file $file does not exist\n";
    print $usage;
}

#parse file
open(IN, "< $file");
my %genes;
while(my $line = <IN>){
    chomp $line;
    my @F = split(/\t/, $line);
    next if(@F < 8);
    next if($F[2] ne 'CDS');

    #genemark by default only fills in the ids and not the names
    my ($g) = $F[8] =~ /gene_id \"([^\"]+)\"/;
       ($g) = $F[8] =~ /gene_name \"([^\"]+)\"/ if(! defined $g);
    my ($t) = $F[8] =~ /transcript_id \"([^\"]+)\"/;
       ($t) = $F[8] =~ /transcript_name \"([^\"]+)\"/ if(! defined $t);

    die "ERROR: Cannot understand format\n".
	"expecting -> gene_id \"xxxx\"\; transcript_id \"xxxx\"\;\n"
	if(! defined $g || ! defined $t);

    #get cintig name
    my $s = $F[0];
    
    #set needed column information
    $genes{$s}{$g}{seqid}  = $F[0] if(! $genes{$s}{$g}{seqid});
    $genes{$s}{$g}{source} = $F[1] if(! $genes{$s}{$g}{source});
    $genes{$s}{$g}{strand} = $F[6] if(! $genes{$s}{$g}{strand});

    $genes{$s}{$g}{mRNA}{$t}{seqid}  = $F[0] if(! $genes{$s}{$g}{mRNA}{$t}{seqid});
    $genes{$s}{$g}{mRNA}{$t}{source} = $F[1] if(! $genes{$s}{$g}{mRNA}{$t}{source});
    $genes{$s}{$g}{mRNA}{$t}{strand} = $F[6] if(! $genes{$s}{$g}{mRNA}{$t}{strand});
    $genes{$s}{$g}{mRNA}{$t}{parent} = $g if(! $genes{$s}{$g}{mRNA}{$t}{parent});

    #set start/end of gene
    $genes{$s}{$g}{B} = $F[3] if(! defined $genes{$s}{$g}{B} || $F[3] < $genes{$s}{$g}{B});
    $genes{$s}{$g}{E} = $F[4] if(! defined $genes{$s}{$g}{E} || $F[4] > $genes{$s}{$g}{E});

    #set start/end of transcript
    $genes{$s}{$g}{mRNA}{$t}{B} = $F[3] if(! defined $genes{$s}{$g}{mRNA}{$t}{B} ||
					   $F[3] < $genes{$s}{$g}{mRNA}{$t}{B}
					   );
    $genes{$s}{$g}{mRNA}{$t}{E} = $F[4] if(! defined $genes{$s}{$g}{mRNA}{$t}{E} ||
					   $F[4] > $genes{$s}{$g}{mRNA}{$t}{E}
					   );
    
    #add CDS to transcript
    my %c = (seqid  => $F[0],
	     source => $F[1],
	     B      => $F[3],
	     E      => $F[4],
	     score  => $F[5],
	     strand => $F[6],
	     phase  => $F[7],
	     parent => $t
	     );

    push (@{$genes{$s}{$g}{mRNA}{$t}{CDS}}, \%c);
}
close(IN);


#build GFF3 structure and dump to file
print "\#\#gff-version 3\n";
gff3_contig(\%genes);

#--------------------------------------------------------------------------
#-------------------------------- SUBS ------------------------------------
#--------------------------------------------------------------------------
sub gff3_contig {
    my $hash = shift;
    foreach my $f (keys %$hash){
	gff3_gene($hash->{$f});
    }
}
#--------------------------------------------------------------------------
sub gff3_gene {
    my $hash = shift;

    foreach my $g (sort {$hash->{$a}{B} <=> $hash->{$b}{B}} keys %$hash) {
	my $gene = $hash->{$g};
	
	print join("\t",$gene->{seqid},$gene->{source},'gene',$gene->{B},
		   $gene->{E},'.',$gene->{strand},'.',sprintf('ID=%s;Name=%s',$g,$g)),"\n";
	
	gff3_mRNA($gene->{mRNA});
    }
}
#--------------------------------------------------------------------------
sub gff3_mRNA {
    my $hash = shift;

    foreach my $t (keys %$hash){
	my $mRNA = $hash->{$t};
	print join("\t",$mRNA->{seqid},$mRNA->{source},"mRNA",
		   $mRNA->{B},$mRNA->{E},'.',$mRNA->{strand},'.',
		   sprintf('ID=%s;Name=%s;Parent=%s',$t,$t,$mRNA->{parent})),"\n";
	
	gff3_CDS($mRNA->{CDS});
    }
}
#--------------------------------------------------------------------------
sub gff3_CDS {
    my $array = shift;

    #define the id
    my $i = 1;

    my @exons;
    my @CDSs;
    foreach my $c (@$array){
	#make exon line
	my $id = $c->{parent} .":exon:". $i;
	my $exon = join("\t",$c->{seqid},$c->{source},"exon",
			$c->{B},$c->{E},'.',$c->{strand},'.',
			sprintf('ID=%s;Name=%s;Parent=%s',$id,$id,$c->{parent}))."\n";
	push(@exons, $exon);
	
	#make CDS line
	$id = $c->{parent} .":CDS:". $i++;
	my $cds = join("\t",$c->{seqid},$c->{source},"CDS",$c->{B},$c->{E},$c->{score},
		       $c->{strand},$c->{phase},sprintf("ID=%s;Name=%s;Parent=%s",$id,$id,$c->{parent}))."\n";
	push(@CDSs, $cds);
    }
    
    #print all exons together then all CDSs together
    print join('', @exons);
    print join('', @CDSs);
}




