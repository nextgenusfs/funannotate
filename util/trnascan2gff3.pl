#!/usr/bin/env perl -w

use strict;
use Bio::Tools::tRNAscanSE;
use Bio::Tools::GFF;

my $out = Bio::Tools::GFF->new(-gff_version => 3);
my $parser = Bio::Tools::tRNAscanSE->new(-file => shift);

# parse the results
while( my $gene = $parser->next_prediction ) {
     $out->write_feature($gene);
     for my $trans ( $gene->get_SeqFeatures() ) {
	 $out->write_feature($trans,$trans->get_SeqFeatures());
    }

}