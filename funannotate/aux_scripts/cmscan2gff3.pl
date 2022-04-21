#!/usr/bin/env perl

#Jason Stajich (2022) 

=head1 NAME

cmscan2gff3.pl - convert raw output of cmscan to gff3

=head1 SYNOPSIS

USAGE: cmscan2gff3.pl --input=/path/to/cmscan.out

=head1 OPTIONS

B<--input,-i>
    The raw output from cmscan:

 #idx target name            accession query name           accession clan name mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc olp anyidx afrct1 afrct2 winidx wfrct1 wfrct2 description of target
 #--- ---------------------- --------- -------------------- --------- --------- --- -------- -------- -------- -------- ------ ----- ---- ---- ----- ------ --------- --- --- ------ ------ ------ ------ ------ ------ ---------------------
 1    LSU_rRNA_archaea       RF02540   NC_013790.1          -         CL00112    cm        1     2990   762872   765862      +    no    1 0.49  45.1 2763.5         0  !   ^       -      -      -      -      -      - -
 2    LSU_rRNA_archaea       RF02540   NC_013790.1          -         CL00112    cm        1     2990  2041329  2038338      -    no    1 0.48  46.1 2755.0         0  !   ^       -      -      -      -      -      - -
 3    LSU_rRNA_bacteria      RF02541   NC_013790.1          -         CL00112    cm        1     2925   762874   765862      +    no    1 0.49  45.1 1872.9         0  !   =       1  1.000  0.999      "      "      " -
 4    LSU_rRNA_bacteria      RF02541   NC_013790.1          -         CL00112    cm        1     2925  2041327  2038338      -    no    1 0.48  46.2 1865.5         0  !   =       2  1.000  0.999      "      "      " -
 5    LSU_rRNA_eukarya       RF02543   NC_013790.1          -         CL00112    cm        1     3401   763018   765851      +    no    1 0.49  41.5 1581.3         0  !   =       1  1.000  0.948      "      "      " -
 6    LSU_rRNA_eukarya       RF02543   NC_013790.1          -         CL00112    cm        1     3401  2041183  2038349      -    no    1 0.49  42.3 1572.1         0  !   =       2  1.000  0.948      "      "      " -
 7    SSU_rRNA_archaea       RF01959   NC_013790.1          -         CL00111    cm        1     1477  2043361  2041888      -    no    1 0.53   4.1 1552.0         0  !   ^       -      -      -      -      -      - -
 8    SSU_rRNA_archaea       RF01959   NC_013790.1          -         CL00111    cm        1     1477   760878   762351      +    no    1 0.54   4.1 1546.5         0  !   ^       -      -      -      -      -      - -
 9    SSU_rRNA_bacteria      RF00177   NC_013790.1          -         CL00111    cm        1     1533  2043366  2041886      -    no    1 0.53   3.7 1161.9         0  !   =       7  0.995  1.000      "      "      " -
 10   SSU_rRNA_bacteria      RF00177   NC_013790.1          -         CL00111    cm        1     1533   760873   762353      +    no    1 0.53   3.7 1156.4         0  !   =       8  0.995  1.000      "      "      " -
 11   SSU_rRNA_eukarya       RF01960   NC_013790.1          -         CL00111    cm        1     1851  2043361  2041891      -    no    1 0.53   4.6  970.4  9.9e-293  !   =       7  1.000  0.998      "      "      " -
 12   SSU_rRNA_eukarya       RF01960   NC_013790.1          -         CL00111    cm        1     1851   760878   762348      +    no    1 0.54   4.5  963.8  9.9e-291  !   =       8  1.000  0.998      "      "      " -
 13   SSU_rRNA_microsporidia RF02542   NC_013790.1          -         CL00111    cm        1     1312  2043361  2041891      -    no    1 0.53   4.6  919.9  7.7e-281  !   =       7  1.000  0.998      "      "      " -
 14   SSU_rRNA_microsporidia RF02542   NC_013790.1          -         CL00111    cm        1     1312   760878   762348      +    no    1 0.54   4.5  917.2  5.4e-280  !   =       8  1.000  0.998      "      "      " -
 15   RNaseP_arch            RF00373   NC_013790.1          -         CL00002    cm        1      303  2614544  2614262      -    no    1 0.43   0.0  184.9   1.1e-53  !   *       -      -      -      -      -      - -
 16   Archaea_SRP            RF01857   NC_013790.1          -         CL00003    cm        1      318  1064321  1064634      +    no    1 0.44   0.1  197.6   6.9e-49  !   *       -      -      -      -      -      - -
 17   FMN                    RF00050   NC_013790.1          -         -          cm        1      140   193975   193837      -    no    1 0.42   0.0  115.2   6.8e-28  !   *       -      -      -      -      -      - -
 18   tRNA                   RF00005   NC_013790.1          -         CL00001    cm        1       71   735136   735208      +    no    1 0.59   0.0   72.1   4.9e-16  !   *       -      -      -      -      -      - -
 19   tRNA                   RF00005   NC_013790.1          -         CL00001    cm        1       71  2350593  2350520      -    no    1 0.66   0.0   71.0     1e-15  !   *       -      -      -      -      -      - -
 20   tRNA                   RF00005   NC_013790.1          -         CL00001    cm        1       71  2680310  2680384      +    no    1 0.52   0.0   70.9   1.1e-15  !   *       -      -      -      -      -      - -
 21   tRNA                   RF00005   NC_013790.1          -         CL00001    cm        1       71  2351254  2351181      -    no    1 0.62   0.0   69.7   2.2e-15  !   *       -      -      -      -      -      - -
 22   tRNA                   RF00005   NC_013790.1          -         CL00001    cm        1       71   361676   361604      -    no    1 0.51   0.0   69.5   2.5e-15  !   *       -      -      -      -      -      - -
 23   tRNA                   RF00005   NC_013790.1          -         CL00001    cm        1       71  2585265  2585193      -    no    1 0.60   0.0   69.2   3.2e-15  !   *       -      -      -      -      -      - -
 24   tRNA                   RF00005   NC_013790.1          -         CL00001    cm        1       71  2585187  2585114      -    no    1 0.59   0.0   68.8   3.9e-15  !   *       -      -      -      -      -      - -
 25   tRNA                   RF00005   NC_013790.1          -         CL00001    cm        1       71  2680159  2680233      +    no    1 0.67   0.0   68.7   4.3e-15  !   *       -      -      -      -      -      - -

B<--log,-l>
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

File converter

=head1  INPUT

Input above.

=head1  OUTPUT

GFF3 to STDOUT

=head1  CONTACT

    Jason Stajich @hyphaltip
    jasonstajich.phd[at]gmail.com

=cut

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage; # could drop this if not a necessary dependency?

my %options = ();
my $results = GetOptions (\%options,
                          'input|i=s',
                          'log|l=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

## open the log if requested
my $logfh;
if (defined $options{log}) {
    open($logfh => ">".$options{log}) || die "can't create log file: $!";
}

## open the input file
open(my $ifh => "<".$options{input}) || die "can't open input file: $!";

# all output needs the gff header
print "##gff-version 3\n";

## globals
my %counts;
## parse the file
while (my $line = <$ifh>) {
    next if /^#/;
    $line =~ s/^\s+//;
    chomp ($line);		# drop trailing whitespace

    my @cols = split( /[\t]/, $line);
    chomp @cols;

    my ($idx,$targetname,$rfam,$contig,$qacc, $clan,$mdl,$mdl_start, $mdl_stop, $query_start,$query_stop, $strand, $truncated, $pass, $gc,
	$bias, $score, $evalue, $inc, $olp, $anyidx, $afrct1,$afrct2,$winidx,$wfrct1,$wfrct2,$target_description) = @cols;

    if ($contig =~ /^(.+?)\s+$/) {
        $contig = $1;
    }
    
    my $ID=sprintf("ID=%s_%d",$targetname,$counts{$rfam}++);
    my $attributes = join(";","rfam_acc=$rfam","product=$targetname","clan=$clan","mdl=$mdl","mdl_start=$mdl_start","mdl_stop=$mdl_stop",
			  "truncated=$truncated","pass=$pass","evalue=$evalue","inc=$inc","olp=$olp");
    print join("\t",$contig,'cmscan_rfam','gene',$query_start,$query_stop,$score,$strand,'.',$ID),"\n";
    print join("\t",$contig,'cmscan_rfam','ncRNA',$query_start,$query_stop,$score,$strand,'.',
	       sprintf("ID=%s_ncRNA;Parent=%s;%s",$ID,$ID,$attributes)),"\n";
    print join("\t",$contig,'cmscan_rfam','exon',$query_start,$query_stop,$score,$strand,'.',
	       sprintf("ID=%s_ncRNA_exon;Parent=%s_ncRNA",$ID,$ID)),"\n";
}

exit(0);

sub  trim { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s };

sub _log {
    my $msg = shift;
    print $logfh "$msg\n" if $logfh;
}

sub check_parameters {
    my $options = shift;
    ## make sure required arguments were passed
    my @required = qw( input );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
    ## handle some defaults
    $options{optional_argument2}   = 'foo'  unless ($options{optional_argument2});
}
