#!/usr/bin/env perl
#
#    pal2nal.pl  (v14)                                      Mikita Suyama
#
#    Usage:  pal2nal.pl  pep.aln  nuc.fasta  [nuc.fasta...]  [options]  >  output
#
#        pep.aln:    protein alignment either in CLUSTAL or FASTA format
#
#        nuc.fasta:  DNA sequences (single multi-fasta or separated files)
#
#        Options:  -output (clustal(default)|paml|fasta|codon)
#                  -nogap
#                  -html
#                  -nostderr
#                  -blockonly
#                  -nomismatch
#                  -codontable N (default=1(universal))
#
#        - IDs in pep.aln are used in the output.
#
#        - Sequence order is automatically checked (see the comment of v12).
#
#        - In-frame stop -> "*" or "_"
#
#        - Frameshift
#
#            Gene    HCDG
#            Pseudo  HC1G     2 del in pseudo
#
#            Gene    HCDG
#            Pseudo  HC2G     1 del in pseudo
#
#            Gene    ND-TY
#            Pseudo  ND1TY    1 ins in pseudo
#
#            Gene    ND-TY
#            Pseudo  ND2TY    2 ins in pseudo
#
#            Gene    EREQK
#            Pseudo  EK4QK    1 ins in pseudo
#
#
#        example:
#
#          pep.aln
#            AA1    ACDEFGARH1G-F*
#            AA2    A2DDW-A*H-G2F*
#
#          nuc.fasta
#            >NUC1
#            GCTTGTGATGAATTTGGTGCTCGTCATGGGGTTTTAA
#            >NUC2
#            GCGGGGACGACTGGGCGTAGCACGGGGGTTTTGA
#
#          output
#            NUC1   GCTTGTGATGAATTTGGTGCTCGTCATG--GGG---TTTTAA
#            NUC2   GCGGG-GACGACTGG---GCGTAGCAC---GGGGG-TTTTGA
#    
#  v14.
#     - NCBI GenBank codon table
#                1  Universal code
#                2  Vertebrate mitochondrial code
#                3  Yeast mitochondrial code
#                4  Mold, Protozoan, and Coelenterate Mitochondrial code
#                   and Mycoplasma/Spiroplasma code
#                5  Invertebrate mitochondrial
#                6  Ciliate, Dasycladacean and Hexamita nuclear code
#                9  Echinoderm and Flatworm mitochondrial code
#               10  Euplotid nuclear code
#               11  Bacterial, archaeal and plant plastid code
#               12  Alternative yeast nuclear code
#               13  Ascidian mitochondrial code
#               14  Alternative flatworm mitochondrial code
#               15  Blepharisma nuclear code
#               16  Chlorophycean mitochondrial code
#               21  Trematode mitochondrial code
#               22  Scenedesmus obliquus mitochondrial code
#               23  Thraustochytrium mitochondrial code
#
#     - Initiation Met is treated separately.
#                                                          2011/12/02
#
#  v13.
#     - file type (Win, UNIX, Mac)
#         new-line:
#            Win         CR+LF     \x0D\x0A
#            UNIX        LF        \x0A
#            Mac         CR        \x0D
#
#               undef $/;
#               s/\x0D\x0A|\x0D|\x0A/\n/g;
#
#         OS:  $^O
#                                                          2010/07/26
#
#     - run bl2seq for inconsistency check if it happens on linux.
#                                                          2008/04/23
#  v12.1
#     - If there is a mismatch, the first atg(lower-case) couldn't
#       recognized as Met.
#                                                          2009/06/22
#
#  v12.
#     - About the input files:
#       Now this script automatically check the following
#       two possibilities:
#         - If the same IDs are used in pep.aln and nuc.fasta,
#             then you don't have to care about the order of
#             the seqs.
#         - If not,
#             the same order of the seqs in pep.aln and nuc.fasta
#             is assumed.
#
#     - In case of 'inconsistency between...', run bl2seq if available
#
#     - Some '1 while' replacement was very slow if the string
#       very long (ca 10k). The routine is modified for speed-up.
#
#                                                          02/02/2007
#  v11.
#     - Mac/DOS file (^M = \r).
#
#     - bug fixed:
#         If aa = '4' (frame shift), 'codon' output
#         had 1 aa shift.
#                                                          09/08/2006
#  v10.
#     - input DNA can be upper/lower cases
#       (in older version, only upper case was allowed).
#
#     - don't report warning if the first Met is GTG
#
#     - default input aln type: clustal
#                                                          01/08/2006
#  v9.3.
#     - "codon" was added to the output options.
#     - if ($html), print "<pre>" and "</pre>".
#                                                          20/06/2006
#
#  v9.2.
#     - pn2codon is modified (to v4)
#       from:  10 residue window (incl '-')
#       to:    10 residue window (excl '-')
#                                                          19/06/2006
#  v9.1.
#     - to handle '/' followed by '*'
#
#    $aaseq[$i] =~ s/\/(-*)[A-Z]/-${1}2/g;
#                              |
#                              v
#    $aaseq[$i] =~ s/\/(-*)[A-Z\*]/-${1}2/g;
#                                                          16/06/2006
#  v9.
#     - Frameshift in tfastx/tfasty
#
#          \  ->  +1
#
#          /  ->  -1  (e.g.  PG/A -> PG-2)
#                                                          06/06/2006
#  v8.
#     - automatic detection of the input alignment format
#       (-a is no longer used.)
#                                                          14/02/2006
#     - new option:
#          -nomismatch  -> remove mismatched codons (mismatch between 
#                          pep and cDNA) from the output.
#
#          -codontable (universal(default)|vmitochondria)
#                                                          26/03/2006
#  v7.
#     - new option:
#   #       -a  input_alignment_type
#   #         clustal  (default)    CLUSTALW format
#   #         fasta
#   #         gblocks               Gblocks txt format *-gb.txt (-p=t)
#   #
#   #           in the case of "-a gblocks", only the codons in conserved
#   #           blocks are converted into codon alignment.
#                                                          27/01/2006
#
#     - if clustal alignment contains
#
#     - new option:
#          -blockonly
#                                                          03/02/2006
#
#     - show warning message only if the position is marked '#'.
#
#     - -i, -j options are deleted.
#                                                          08/02/2006
#  v6.
#   #  - new option:
#   #      -fasta        -> set this option if input pep file is
#   #                       not in *.aln format but in *.fasta format.
#   #                                                       16/08/2005
#
#  v4.
#     - new option:
#         -html         ->  for the web server (STDERR messages -> STDOUT)
#                                              (mismatch in color)
#         -nostderr     ->  for Ks,Ka calc on the web server (NO STDERR messages)
#                                                          25/04/2005
#
#  v3.
#     - renamed back to "pal2nal.pl"
#     - subroutine (pn2codon) is modified:
#         - reverse translation table (%p2c) was replaced
#         - if there is no exact match, try fragment anchors
#         - return value is changed to hash
#     - new options:
#         -output clustal|paml|fasta
#                   clustal format (default)
#                   paml
#                   fasta
#         -nogap        ->  remove all gaped codons and stop codons
#
#                                                          22/04/2005
#  v2.
#     - nuc.fasta can be either one multiple fasta file or several files.
#       The order of seqs should be the same as those in the pep.aln file.
#                                                          13/07/2004
#
#       old script name: pal2nal.pl
#                                       06/09/2000    Mikita Suyama
#
#     - modified for multiple seqs.     17/10/2000
#
#     - AA-ids and NUC-ids must be the same -> TURNED OFF  05/05/2002
#       Now, only the order is the matter.
#       AA-ids are used throughout; NUC-ids are not used at all.
#
#     - old condition: AA = NUC
#       new condition: AA <= NUC   ##  NUC(cDNA) can be longer than AA  ##
#
#     - renamed from "pal2nal.pl" to codon_aln.pl          12/07/2004
#



$| = 1;

undef $/;
$myos = $^O;

if ($#ARGV < 1) {
    & showhelp();
    exit;
} else {
    $getoutform = 0;
    $nogap = 0;
    $html = 0;
    $nostderr = 0;
    $fasta = 0;
    undef($alnfile);
    undef(@nucfiles);
    $outform = "clustal";
    $blockonly = 0;
    $nomismatch = 0;
    $getcodontable = 0;
    $codontable = 1;
    foreach $i (0..$#ARGV) {
        if ($ARGV[$i] eq "-h") {
            & showhelp();
        } elsif ($ARGV[$i] eq "-output") {
            $getoutform = 1;
        } elsif ($getoutform) {
            $outform = $ARGV[$i];
            if ($outform ne "clustal" &&
                $outform ne "paml" &&
                $outform ne "fasta" &&
                $outform ne "codon") {
                print STDERR "\nERROR:  valid output format: clustal, paml, fasta, or codon\n\n";
                exit;
            }
            $getoutform = 0;
        } elsif ($ARGV[$i] eq "-blockonly") {
            $blockonly = 1;
        } elsif ($ARGV[$i] eq "-nogap") {
            $nogap = 1;
        } elsif ($ARGV[$i] eq "-nomismatch") {
            $nomismatch = 1;
        } elsif ($ARGV[$i] eq "-codontable") {
            $getcodontable = 1;
        } elsif ($getcodontable) {
            $codontable = $ARGV[$i];
            if ($codontable != 1 && $codontable != 2 && $codontable != 3 &&
                $codontable != 4 && $codontable != 5 && $codontable != 6 &&
                $codontable != 9 && $codontable != 10 && $codontable != 11 &&
                $codontable != 12 && $codontable != 13 && $codontable != 14 &&
                $codontable != 15 && $codontable != 16 && $codontable != 21 &&
                $codontable != 22 && $codontable != 23) {
                print STDERR "\nERROR:  invalid codontable number, $codontable!!\n\n";
                exit;
            }
            $getcodontable = 0;
        } elsif ($ARGV[$i] eq "-html") {
            $html = 1;
        } elsif ($ARGV[$i] eq "-nostderr") {
            $nostderr = 1;
        } elsif (!$alnfile) {
            $alnfile = $ARGV[$i];
        } else {
            push(@nucfiles, $ARGV[$i]);
        }
    }
}

if ($html) {
    print "<pre>\n";
}

#---------------#
# check options  ("-outform codon" is not valid with -blockonly, -nogap, -nomismatch.)
#---------------#

if ($outform eq "codon" &&
    ($blockonly || $nogap || $nomismatch)) {
    if ($html) {
        print "\nERROR:  if \"codon(Output format)\" is selected, don't use \"Remove gaps, inframe stop codons\" or \"Remove mismatches\" or \"Use only selected positions\".\n\n";
    } else {
        print STDERR "\nERROR:  \"-outform codon\" is not valid with -blockonly, -nogap, -nomismatch\n\n";
    }
    exit;
}


#---------------------#
#  Get nuc sequences
#---------------------#

undef(@nucid);
undef(%id2nucseq);
$nseq = -1;
foreach $i (0..$#nucfiles) {
    open(NUCFILE, "< $nucfiles[$i]") || die "Can't open $nucfiles[$i]";
    $nucfiledata = <NUCFILE>;
    close(NUCFILE);
    $nucfiledata =~ s/\x0D\x0A|\x0D|\x0A/\n/g;
    foreach (split(/\n/, $nucfiledata)) {
        if (!/^#/ && /\S+/) {
            if (/^>(\S+)/) {
                ++$nseq;
                $tmpnucid = $1;
                push(@nucid, $tmpnucid);
            } else {
                s/[^a-zA-Z]//g;
                $id2nucseq{$tmpnucid} .= $_;
            }
        }
    }
}


#-------------------#
#  Get aa alignemt
#-------------------#

undef(@aaid);
undef(%id2aaaln);
undef(%aaidcnt);
undef(@aaseq);
undef($gblockseq);

$gettype = 1;
open(ALNFILE, "< $alnfile") || die "Can't open $alnfile";
while (<ALNFILE>) {
    chomp;
    if ($gettype && !/^#/ && /\S+/) {
        if (/^CLUSTAL/) {
            $inalntype = "clustal";
        } elsif (/^>/) {
            $inalntype = "fasta";
        } elsif (/^Gblocks/) {
            $inalntype = "gblocks";
        } else {
            $inalntype = "clustal";
        }

        $gettype = 0;
    }
}
close(ALNFILE);

open(ALNFILE, "< $alnfile") || die "Can't open $alnfile";
$getblock = 0;
$aafiledata = <ALNFILE>;
close(ALNFILE);
$aafiledata =~ s/\x0D\x0A|\x0D|\x0A/\n/g;
foreach (split(/\n/, $aafiledata)) {
    if ($inalntype eq "clustal") {
        if (!/^CLUSTAL/ && /^\S+/ && !/^#/) {
            s/\s+$//;
            @dat = split(/\s+/, $_);
            ++$aaidcnt{$dat[0]};
            push(@aaid, $dat[0]) if ($aaidcnt{$dat[0]} == 1);
            $dat[1] =~ tr/a-z/A-Z/;
            $id2aaaln{$dat[0]} .= $dat[1];

            $tmplen = length($_);

            /^\S+\s+/;
            $idspc = length($&);
            $subalnlen = length($dat[1]);

            $getblock = 1;
        } elsif ($getblock) {
            $_ .= ' ' x ($tmplen - length($_));
            $gblockseq .= substr($_, $idspc, $subalnlen);

            $getblock = 0;
        }
    } elsif ($inalntype eq "fasta") {
        if (/^>(\S+)/) {
            $tmpid = $1;
            push(@aaid, $tmpid);
        } else {
            s/\s+//g;
            tr/a-z/A-Z/;
            $id2aaaln{$tmpid} .= $_;
        }
    } elsif ($inalntype eq "gblocks") {
        $getaln = 0;
        if (/^\s+\=/) {
            $getaln = 1;
        } elsif (/^Parameters/) {
            $getaln = 0;
        } elsif ($getaln) {
            @dat = split(/\s+/, $_);
            if (/^Gblocks/) {
                $gblockseq .= $dat[1];
            } elsif (/^\S+/) {
                ++$aaidcnt{$dat[0]};
                push(@aaid, $dat[0]) if ($aaidcnt{$dat[0]} == 1);
                $dat[1] =~ tr/a-z/A-Z/;
                $id2aaaln{$dat[0]} .= $dat[1];
            }
        }
    }
}

foreach $i (0..$#aaid) {
    push(@aaseq, $id2aaaln{$aaid[$i]});
}


#------------------------------------#
#  For frameshifts in tfastx/tfasty
#------------------------------------#

foreach $i (0..$#aaid) {
    $aaseq[$i] =~ s/\\/1/g;
    $aaseq[$i] =~ s/\/(-*)[A-Z\*]/-${1}2/g;
}


#-------------------------------------#
#  Check the input seqs (pep <=> nuc)
#-------------------------------------#

if ($#aaid != $#nucid) {
    $naa = $#aaid + 1;
    $nnuc = $#nucid + 1;
    if ($html) {
        print "\nERROR: number of input seqs differ (aa: $naa;  nuc: $nnuc)!!\n\n";
    } else {
        print STDERR "\nERROR: number of input seqs differ (aa: $naa;  nuc: $nnuc)!!\n\n";
        print STDERR "   aa  '@aaid'\n";
        print STDERR "   nuc '@nucid'\n";
    }
    exit;
}


#---------------------------------#
# corresponding IDs or same order
#---------------------------------#

@commonids = & common_elem(*aaid, *nucid);

sub common_elem {
    local(*aarr, *barr) = @_;

    local(%mark, @comarr);
    grep($mark{$_}++, @aarr);
    @comarr = grep($mark{$_}, @barr);
}
        
if ($#commonids == $#aaid) {
    $idcorrespondence = "sameID";
} else {
    $idcorrespondence = "ordered";
}


#-------------------#
#  Codon sequences
#-------------------#

undef(@codonseq);
undef(%aaidpos2mismatch);
undef(@outmessage);
foreach $i (0..$#aaid) {
    if ($idcorrespondence eq "sameID") {
        $tmpnucid = $aaid[$i];
    } elsif ($idcorrespondence eq "ordered") {
        $tmpnucid = $nucid[$i];
    } else {
        print STDERR "\nERROR in ID correspondence.\n\n";
        exit;
    }

    %codonout = & pn2codon($aaseq[$i], $id2nucseq{$tmpnucid}, $codontable);
    @message = @{$codonout{'message'}};

    foreach $j (0..$#message) {
        $outl = "WARNING: $aaid[$i] $message[$j]";
        push(@outmessage, $outl);
        @dat = split(/\s+/, $message[$j]);
        $dat[1] =~ s/:$//;
        $aaidpos2mismatch{"$aaid[$i] $dat[1]"} = 1;
    }

    if ($codonout{'result'} == 1 || $codonout{'result'} == 2) {
        push(@codonseq, $codonout{'codonseq'});
    } else {
        $erraa = $aaseq[$i];
        $erraa =~ s/-//g;
        # 1 while $erraa =~ s/(.{60})(.+)/$1\n$2/;
        @frags = & my1while($erraa, 60);
        $erraa = join("\n", @frags);

        $errnuc = $id2nucseq{$tmpnucid};
        $errnuc =~ s/-//g;
        # 1 while $errnuc =~ s/(.{60})(.+)/$1\n$2/;
        @frags = & my1while($errnuc, 60);
        $errnuc = join("\n", @frags);

        if ($html) {
            print "</pre>\n";
            print "<H1>ERROR in your input files!</H1>\n";
            print "<pre>\n";
            print "#---  ERROR: inconsistency between the following pep and nuc seqs  ---#\n";

            print ">$aaid[$i]\n";
            print "$erraa\n";

            print ">$tmpnucid\n";
            print "$errnuc\n";

            print "</pre>\n";
            $tmpwhich = "tmp/tmpwhich.$$";
        } else {
            print STDERR "#---  ERROR: inconsistency between the following pep and nuc seqs  ---#\n";

            print STDERR ">$aaid[$i]\n";
            print STDERR "$erraa\n";

            print STDERR ">$tmpnucid\n";
            print STDERR "$errnuc\n";

            $tmpwhich = "tmpwhich.$$";
        }

        if ($myos eq "linux") {
            system("which bl2seq > $tmpwhich");

            $foundbl2seq = 0;
            open(TMPWHICH, "< $tmpwhich") || die "Can't open $tmpwhich";
            while (<TMPWHICH>) {
                chomp;
                $foundbl2seq = 1 if (/bl2seq$/);
            }
            close(TMPWHICH);
            unlink($tmpwhich);

            if ($foundbl2seq) {


                #------------------------------------------#
                # run bl2seq (if the command is available)
                #------------------------------------------#

                if ($html) {
                    $erraafile = "tmp/erraafile.$$";
                    $errnucfile = "tmp/errnucfile.$$";
                    $tblnout = "tmp/tbln.out.$$";
                } else {
                    $erraafile = "erraafile.$$";
                    $errnucfile = "errnucfile.$$";
                    $tblnout = "tbln.out.$$";
                }

                open(ERRAAFILE, "> $erraafile") || die "Can't open $erraafile";
                print ERRAAFILE ">$aaid[$i]\n";
                print ERRAAFILE "$erraa\n";
                close(ERRAAFILE);

                open(ERRNUCFILE, "> $errnucfile") || die "Can't open $errnucfile";
                print ERRNUCFILE ">$tmpnucid\n";
                print ERRNUCFILE "$errnuc\n";
                close(ERRNUCFILE);

                system("bl2seq -p tblastn -F F -i $erraafile -j $errnucfile -o $tblnout");

                if ($html) {
                    print "<BR>\n";
                    print "<BR>\n";
                    print "<H1>Check the following TBLASTN output.</H1><BR>\n";
                    print "<pre>\n";
                    print "      your pep -> 'Query'\n";
                    print "      your nuc -> 'Sbjct'\n";
                    print "<BR>\n";
                } else {
                    print STDERR "\n";
                    print STDERR "\n";
                    print STDERR "        ###-----   Check the following TBLASTN output:           -----###\n";
                    print STDERR "        ###-----       your pep -> 'Query'                       -----###\n";
                    print STDERR "        ###-----       your nuc -> 'Sbjct'                       -----###\n";
                    print STDERR "\n";
                }

                open(TBLNOUT, "< $tblnout") || die "Can't open $tblnout";
                $tblnoutdata = <TBLNOUT>;
                close(TBLNOUT);
                $tblnoutdata =~ s/\x0D\x0A|\x0D|\x0A/\n/g;
                foreach (split(/\n/, $tblnoutdata)) {
                    if ($html) {
                        print "$_\n";
                    } else {
                        print STDERR "$_\n";
                    }
                }

                if ($html) {
                    print "</pre>\n";
                }

                unlink($erraafile);
                unlink($errnucfile);
                unlink($tblnout);
            } else {
                if ($html) {
                    print "<BR>\n";
                    print "<BR>\n";
                    print "Run bl2seq (-p tblastn) or GeneWise to see the inconsistency.<BR>\n";
                    print "<BR>\n";
                } else {
                    print STDERR "\n";
                    print STDERR "\n";
                    print STDERR "Run bl2seq (-p tblastn) or GeneWise to see the inconsistency.\n";
                    print STDERR "\n";
                }
            }
        } else {    ####    non-linux environment
            print STDERR "\n";
            print STDERR "\n";
            print STDERR "Run bl2seq (-p tblastn) or GeneWise to see the inconsistency.\n";
            print STDERR "\n";
        }
        exit;
    }
}


#-------------------#
#  Warning in '#'?
#-------------------#

if ($gblockseq =~ /#/ && $blockonly) {
    undef(@newoutmessage);
    foreach $i (0..$#outmessage) {
        $mpos = (split(/\s+/, $outmessage[$i]))[3];
        $mpos =~ s/://;

        if (substr($gblockseq, $mpos - 1, 1) eq "#") {
            push(@newoutmessage, $outmessage[$i]);
        }
    }
    @outmessage = @newoutmessage;
}


#--------------------#
#  Warning messages
#--------------------#

if ($codontable != 1) {
    # push(@outmessage, "Codon table $codontable is used");
    splice(@outmessage, 0, 0, "Codontable $codontable is used");
}

if (!$nostderr) {
    foreach $j (0..$#outmessage) {
        if ($html) {
            if ($j == 0 && !$nomismatch) {
                print "#------------------------------------------------------------------------#\n";
            }
            if ($nomismatch) {
                # print "#  $outmessage[$j]  (excluded from the output)\n";
            } else {
                print "#  $outmessage[$j]\n";
            }
            if ($j == $#outmessage && !$nomismatch) {
                print "#------------------------------------------------------------------------#\n\n";
            }
        } else {
            if ($j == 0 && !$nomismatch) {
                print STDERR "#------------------------------------------------------------------------#\n";
                print STDERR "#  Input files:  $alnfile @nucfiles\n";
            }
            if ($nomismatch) {
                # print STDERR "#  $outmessage[$j]  (exlucded from the output)\n";
            } else {
                print STDERR "#  $outmessage[$j]\n";
            }
            if ($j == $#outmessage && !$nomismatch) {
                print STDERR "#------------------------------------------------------------------------#\n\n";
            }
        }
    }
}

undef(%errorpos);
foreach $j (0..$#outmessage) {
    @dat = split(/\s+/, $outmessage[$j]);
    $tmperrpos = $dat[3];
    $tmperrpos =~ s/:$//;
    $tmperrpos -= 1;
    $errorpos{$tmperrpos} = 1;
}


#----------------------------------#
#  Make an AA-based NUC alignment
#----------------------------------#

undef(@tmppos);
undef(@codonaln);
undef(@coloraln);
undef($maskseq);
foreach $i (0..length($aaseq[0]) - 1) {
    $tmpmax = 0;
    $apos = $i + 1;

    #-----------#
    # gblocks ?
    #-----------#

    $putcodon = 1;
    if ($gblockseq =~ /#/ && substr($gblockseq, $i, 1) ne "#" && $blockonly) {
        $putcodon = 0;
    }
    if ($nomismatch && $errorpos{$i}) {
        $putcodon = 0;
    }

    foreach $k (0..$#aaid) {
        $tmpaa = substr($aaseq[$k], $i, 1);
        if ($tmpaa !~ /\d/) {
            $tmplen = 3;
        } else {
            $tmplen = (int(($tmpaa - 1) / 3) + 1) * 3;  # 1, 2, 3 -> 3
                                                        # 4, 5, 6 -> 6
                                                        # 7, 8, 9 -> 9
        }
        $tmpmax = $tmplen if ($tmpmax < $tmplen);
    }
    foreach $k (0..$#aaid) {
        $tmpaa = substr($aaseq[$k], $i, 1);
        if ($tmpaa !~ /\d/) {
            if ($tmpaa eq '-') {
                # if ($putcodon || (!$blockonly && !$nomismatch)) {
                if ($putcodon) {
                    $codonaln[$k] .= '-' x $tmpmax;
                    if ($aaidpos2mismatch{"$aaid[$k] $apos"}) {
                        $coloraln[$k] .= 'R' x $tmpmax;
                    } else {
                        $coloraln[$k] .= '-' x $tmpmax;
                    }
                }
            } elsif ($tmpaa =~ /[A-Z]/ || $tmpaa eq '*') {
                # if ($putcodon || (!$blockonly && !$nomismatch)) {
                if ($putcodon) {
                    $codonaln[$k] .= substr($codonseq[$k], $tmppos[$k], 3);
                    if ($aaidpos2mismatch{"$aaid[$k] $apos"}) {
                        $coloraln[$k] .= 'RRR';
                    } else {
                        $coloraln[$k] .= '---';
                    }
                }
                $tmppos[$k] += 3;
                # if ($putcodon || (!$blockonly && !$nomismatch)) {
                if ($putcodon) {
                    $codonaln[$k] .= '-' x ($tmpmax - 3) if ($tmpmax - 3 > 0);
                    if ($aaidpos2mismatch{"$aaid[$k] $apos"}) {
                        $coloraln[$k] .= 'R' x ($tmpmax - 3) if ($tmpmax - 3 > 0);
                    } else {
                        $coloraln[$k] .= '-' x ($tmpmax - 3) if ($tmpmax - 3 > 0);
                    }
                }
            }
        } elsif ($tmpaa =~ /\d/) {
            # if ($putcodon || (!$blockonly && !$nomismatch)) {
            if ($putcodon) {
                $codonaln[$k] .= substr($codonseq[$k], $tmppos[$k], $tmpaa);
                if ($aaidpos2mismatch{"$aaid[$k] $apos"}) {
                    $coloraln[$k] .= 'R' x $tmpaa;
                } else {
                    $coloraln[$k] .= '-' x $tmpaa;
                }
            }
            $tmppos[$k] += $tmpaa;
            # if ($putcodon || (!$blockonly &&!$nomismatch)) {
            if ($putcodon) {
                $codonaln[$k] .= '-' x ($tmpmax - $tmpaa);
                if ($aaidpos2mismatch{"$aaid[$k] $apos"}) {
                    $coloraln[$k] .= 'R' x ($tmpmax - $tmpaa);
                } else {
                    $coloraln[$k] .= '-' x ($tmpmax - $tmpaa);
                }
            }
        }
    }
    if (!$blockonly) {
        # if ($putcodon && substr($gblockseq, $i, 1) eq "#") {
        #     $maskseq .= '#' x $tmpmax;
        # } else {
        #     $maskseq .= ' ' x $tmpmax;
        # }
        if ($putcodon) {
            if (substr($gblockseq, $i, 1) eq "#") {
                $maskseq .= "#" x $tmpmax;
            } else {
                $maskseq .= " " x $tmpmax;
            }
        }
    }
}


#-----------#
#  -nogap?
#-----------#

$alilen = length($codonaln[0]);

if ($nogap) {
    $tmppos = 0;
    undef(@nogapaln);
    undef(@nogapcoloraln);
    undef($nogapmaskseq);
    while ($tmppos < $alilen) {
        $outok = 1;
        foreach $i (0..$#codonaln) {
            $tmpcodon = substr($codonaln[$i], $tmppos, 3);
            $outok = 0 if ($tmpcodon =~ /-/);
            $outok = 0 if ($tmpcodon =~ /(((U|T)A(A|G|R))|((T|U)GA))/);
        }
        if ($outok) {
            foreach $i (0..$#codonaln) {
                $tmpcodon = substr($codonaln[$i], $tmppos, 3);
                $nogapaln[$i] .= $tmpcodon;
                $tmpcolorcodon = substr($coloraln[$i], $tmppos, 3);
                $nogapcoloraln[$i] .= $tmpcolorcodon;
            }
            $nogapmaskseq .= substr($maskseq, $tmppos, 3);
        }
        $tmppos += 3;
    }
    @codonaln = @nogapaln;
    @coloraln = @nogapcoloraln;
    $maskseq = $nogapmaskseq;
}


#----------#
#  Output
#----------#

$maxn = 0;
foreach $i (0..$#aaid) {
    $maxn = length($aaid[$i]) if ($maxn < length($aaid[$i]));
}
$maxn = 10 if ($maxn < 10);

$alilen = length($codonaln[0]);

# foreach $i (0..$#aaid) {
#     1 while $codonaln[$i] =~ s/(.{60})(.+)/$1\n$2/;
#     1 while $coloraln[$i] =~ s/(.{60})(.+)/$1\n$2/;
# }
# 1 while $maskseq =~ s/(.{60})(.+)/$1\n$2/;

undef(%aaid2codonarr);
undef(%aaid2colorarr);
foreach $i (0..$#aaid) {
    @{$aaid2codonarr{$aaid[$i]}} = & my1while($codonaln[$i], 60);
    @{$aaid2colorarr{$aaid[$i]}} = & my1while($coloraln[$i], 60);
}
@maskarr = & my1while($maskseq, 60);


if ($outform eq "clustal") {

    #-----------#
    #  clustal
    #-----------#

    print "CLUSTAL W multiple sequence alignment\n";
    print "\n";

    @output1 = @{$aaid2codonarr{$aaid[0]}};
    if ($html) {
        foreach $i (0..$#output1) {
            foreach $k (0..$#aaid) {
                printf "%-${maxn}s    ", $aaid[$k];
                $outf = ${$aaid2codonarr{$aaid[$k]}}[$i];
                $outr = ${$aaid2colorarr{$aaid[$k]}}[$i];
                $rlen = length($outf);
                if ($outr =~ /R/) {
                    foreach $l (0..$rlen - 1) {
                        $tmpnuc = substr($outf, $l, 1);
                        $tmpr = substr($outr, $l, 1);
                        if ($tmpr eq "R") {
                            print "<FONT color='red'>$tmpnuc</FONT>";
                        } else {
                            print "$tmpnuc";
                        }
                    }
                    print "\n";
                } else {
                    print "$outf\n";
                }
            }
            $outmask = $maskarr[$i];
            printf "%-${maxn}s    $outmask\n", ' ' if (!$blockonly && $gblockseq =~ /#/);
            print "\n";
        }
    } else {
        foreach $i (0..$#output1) {
            foreach $k (0..$#aaid) {
                $outf = ${$aaid2codonarr{$aaid[$k]}}[$i];
                printf "%-${maxn}s    $outf\n", $aaid[$k];
            }
            $outmask = $maskarr[$i];
            printf "%-${maxn}s    $outmask\n", ' ' if (!$blockonly && $gblockseq =~ /#/);
            print "\n";
        }
    }

} elsif ($outform eq "codon") {

    #-------#
    # codon
    #-------#

    # @outaa = @aaseq;

    undef(@outaa);

    $withn = 0;
    foreach $j (0..$#aaid) {
        $withn = 1 if ($aaseq[$j] =~ /\d/);
    }

    if ($withn) {
        $alnlen = length($aaseq[0]);
        foreach $i (0..$alnlen - 1) {
            $maxaan = 0;
            foreach $j (0..$#aaid) {
                $tmpaa = substr($aaseq[$j], $i, 1);
                if ($tmpaa =~ /\d/ && $tmpaa > $maxaan) {
                    $maxaan = $tmpaa;
                }
            }
            if ($maxaan >= 4) {
                $tmplen = int(($tmpaa - 1) / 3) + 1;  # 4, 5, 6 -> 2
                                                      # 7, 8, 9 -> 3
            } else {
                $tmplen = 1;
            }
            foreach $j (0..$#aaid) {
                $pushaa = '-' x $tmplen;
                $tmpaa = substr($aaseq[$j], $i, 1);
                substr($pushaa, 0, 1) = $tmpaa;
                $outaa[$j] .= $pushaa;
            }
        }
    } else {
        @outaa = @aaseq;
    }

    undef(%aaid2aaarr);
    foreach $i (0..$#aaid) {
        # 1 while $outaa[$i] =~ s/(.{20})(.+)/$1\n$2/;
        @{$aaid2aaarr{$aaid[$i]}} = & my1while($outaa[$i], 20);
    }

    @output1 = @{$aaid2codonarr{$aaid[0]}};
    if ($html) {
        foreach $i (0..$#output1) {
            foreach $k (0..$#aaid) {
                printf "%-${maxn}s    ", '';
                $outa = ${$aaid2aaarr{$aaid[$k]}}[$i];
                # 1 while $outa =~ s/(\S)(\S+)/$1   $2/;
                @frags = & my1while($outa, 1);
                $outa = join("   ", @frags);

                $outf = ${$aaid2codonarr{$aaid[$k]}}[$i];
                # 1 while $outf =~ s/(\S{3})(\S+)/$1 $2/;
                @frags = & my1while($outf, 3);
                $outf = join(" ", @frags);

                $outr = ${$aaid2colorarr{$aaid[$k]}}[$i];
                # 1 while $outr =~ s/(\S{3})(\S+)/$1 $2/;
                @frags = & my1while($outr, 3);
                $outr = join(" ", @frags);

                if ($outr =~ /R/) {
                    $rlen = length($outf);
                    foreach $l (0..$rlen - 1) {
                        $tmppep = substr($outa, $l, 1);
                        $tmpr = substr($outr, $l, 1);
                        if ($tmpr eq "R") {
                            print "<FONT color='red'>$tmppep</FONT>";
                        } else {
                            print "$tmppep";
                        }
                    }
                    print "\n";
                    printf "%-${maxn}s    ", $aaid[$k];
                    foreach $l (0..$rlen - 1) {
                        $tmpnuc = substr($outf, $l, 1);
                        $tmpr = substr($outr, $l, 1);
                        if ($tmpr eq "R") {
                            print "<FONT color='red'>$tmpnuc</FONT>";
                        } else {
                            print "$tmpnuc";
                        }
                    }
                    print "\n";
                } else {
                    print "$outa\n";
                    printf "%-${maxn}s    ", $aaid[$k];
                    print "$outf\n";
                }
            }
            $outmask = $maskarr[$i];
            $outmask =~ s/\s/-/g;
            # 1 while $outmask =~ s/(\S{3})(\S+)/$1 $2/;
            @frags = & my1while($outmask, 3);
            $outmask = join(" ", @frags);
            $outmask =~ s/-/ /g;
            printf "%-${maxn}s    $outmask\n", ' ' if (!$blockonly && $gblockseq =~ /#/);
            print "\n";
        }
    } else {
        foreach $i (0..$#output1) {
            foreach $k (0..$#aaid) {
                $outa = ${$aaid2aaarr{$aaid[$k]}}[$i];
                # 1 while $outa =~ s/(\S)(\S+)/$1   $2/;
                @frags = & my1while($outa, 1);
                $outa = join("   ", @frags);
                printf "%-${maxn}s    $outa\n", '';

                $outf = ${$aaid2codonarr{$aaid[$k]}}[$i];
                # 1 while $outf =~ s/(\S{3})(\S+)/$1 $2/;
                @frags = & my1while($outf, 3);
                $outf = join(" ", @frags);
                printf "%-${maxn}s    $outf\n", $aaid[$k];
            }
            $outmask = $maskarr[$i];
            $outmask =~ s/\s/-/g;
            # 1 while $outmask =~ s/(\S{3})(\S+)/$1 $2/;
            @frags = & my1while($outmask, 3);
            $outmask = join(" ", @frags);
            $outmask =~ s/-/ /g;
            printf "%-${maxn}s    $outmask\n", ' ' if (!$blockonly && $gblockseq =~ /#/);
            print "\n";
        }
    }
} elsif ($outform eq "paml") {

    #--------#
    #  paml
    #--------#

    $nseq = $#aaid + 1;
    printf " %3d %6d\n", $nseq, $alilen;
    foreach $i (0..$#aaid) {
        print  "$aaid[$i]\n";
        if ($html) {
            @outf = @{$aaid2codonarr{$aaid[$i]}};
            @outr = @{$aaid2colorarr{$aaid[$i]}};
            foreach $k (0..$#outf) {
                if ($outr[$k] =~ /R/) {
                    $lenf = length($outf[$k]);
                    foreach $l (0..$lenf - 1) {
                        $tmpnuc = substr($outf[$k], $l, 1);
                        $tmpr = substr($outr[$k], $l, 1);
                        if ($tmpr eq "R") {
                            print "<FONT color='red'>$tmpnuc</FONT>";
                        } else {
                            print "$tmpnuc";
                        }
                    }
                    print "\n";
                } else {
                    print "$outf[$k]\n";
                }
            }
        } else {
            $outcodon = join("\n", @{$aaid2codonarr{$aaid[$i]}});
            print  "$outcodon\n";
        }
    }
} elsif ($outform eq "fasta") {

    #---------#
    #  fasta
    #---------#

    foreach $i (0..$#aaid) {
        print  ">$aaid[$i]\n";
        if ($html) {
            @outf = @{$aaid2codonarr{$aaid[$i]}};
            @outr = @{$aaid2colorarr{$aaid[$i]}};
            foreach $k (0..$#outf) {
                if ($outr[$k] =~ /R/) {
                    $lenf = length($outf[$k]);
                    foreach $l (0..$lenf - 1) {
                        $tmpnuc = substr($outf[$k], $l, 1);
                        $tmpr = substr($outr[$k], $l, 1);
                        if ($tmpr eq "R") {
                            print "<FONT color='red'>$tmpnuc</FONT>";
                        } else {
                            print "$tmpnuc";
                        }
                    }
                    print "\n";
                } else {
                    print "$outf[$k]\n";
                }
            }
        } else {
            $outcodon = join("\n", @{$aaid2codonarr{$aaid[$i]}});
            print  "$outcodon\n";
        }
    }
}

if ($html) {
    print "</pre>\n";
}


#-----------------------------------------------------------------------

sub pn2codon {
    #    pn2codon v6
    #
    #    input:   $pep    protein sequence
    #                         termination -> "_" or "*";
    #                         frameshift  -> digit
    #                         "-" or "."  -> gap
    #             $nuc    DNA or RNA sequence (lower/upper case letters)
    #
    #             $codontable  (corresponds to codon tables used in GenBank
    #                1  Universal code
    #                2  Vertebrate mitochondrial code
    #                3  Yeast mitochondrial code
    #                4  Mold, Protozoan, and Coelenterate Mitochondrial code
    #                   and Mycoplasma/Spiroplasma code
    #                5  Invertebrate mitochondrial
    #                6  Ciliate, Dasycladacean and Hexamita nuclear code
    #                9  Echinoderm and Flatworm mitochondrial code
    #               10  Euplotid nuclear code
    #               11  Bacterial, archaeal and plant plastid code
    #               12  Alternative yeast nuclear code
    #               13  Ascidian mitochondrial code
    #               14  Alternative flatworm mitochondrial code
    #               15  Blepharisma nuclear code
    #               16  Chlorophycean mitochondrial code
    #               21  Trematode mitochondrial code
    #               22  Scenedesmus obliquus mitochondrial code
    #               23  Thraustochytrium mitochondrial code
    #
    #    return:  hash
    #                $retval{'codonseq'}: codon seq (w/o gap)
    #                $retval{'message'}:  error/warning messame (array)
    #                $retval{'result'}:   1: OK, 2: mismatch, -1: no match found
    #
    #
    #                                      05/05/2002    Mikita Suyama
    #                                      12/07/2004
    #
    #    v2                                22/04/2005
    #    - reverse translation table (%p2c) was replaced
    #    - if there is no exact match, try fragment anchors
    #    - return value is changed to hash
    #
    #    v3                                26/03/2006
    #    - codon table for vertebrate mitochondria (vmitochondria) is added
    #
    #    v4
    #    - from: 10 residue window (incl '-')
    #      to:   10 residue window (excl '-')
    #                                      19/06/2006
    #    v5
    #    - bug fix:
    #        just before 'does not correspond'
    #        /$p2c{$tmpaa}/  ->   /$p2c{$tmpaa}/i
    #                                           ^
    #    - the first Met can be (A|G|C|R)TG
    #                                      01/08/2006
    #    v5.1
    #    - but fix:
    #        if ($tmpcodon !~ /((A|C|G|R)TG)/i) {
    #                                        ^
    #                                      2009/06/22
    #    v6
    #    - start Met -> aa-code "B" in %p2c
    #    - NCBI codon tables
    #                                      2011/11/13
    #      


    local($pep, $nuc, $codontable) = @_;


    local(%p2c);
    local($peplen, $qcodon, $codon);
    local($tmpaa, $message, $modpep, @qcodon, @fncodon, $wholecodon);
    local($i, $j, $anclen, @anchor, $peppos, $tmpcodon, $codonpos);

    local(%retval);


    if ($codontable == 1) {
        #-----------#
        # Universal Code (transl_table=1)
        #-----------#
        %p2c = (
            "B" => "((U|T|C|Y|A)(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "(((U|T)A(A|G|R))|((T|U)GA))",
            "*" => "(((U|T)A(A|G|R))|((T|U)GA))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    } elsif ($codontable == 2) {
        #--------------------------#
        # Vertebrate Mitochondrial Code (transl_table=2)
        #--------------------------#
        %p2c = (
            "B" => "((A(U|T).)|G(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "(CG.)",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y))",
            "_" => "(((U|T)A(A|G|R))|(AG(A|G|R)))",
            "*" => "(((U|T)A(A|G|R))|(AG(A|G|R)))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)(A|G|R))",
            "W" => "((U|T)G(A|G|R))",
            "X" => "...",
        );
    } elsif ($codontable == 3) {
        #--------------------------#
        # Yeast Mitochondrial Code (transl_table=3)
        #--------------------------#
        %p2c = (
            "B" => "(A(U|T)(A|G|R))",
            "L" => "((U|T)(U|T)(A|G|R))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "((AC.)|(C(U|T).))",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y))",
            "_" => "((U|T)A(A|G|R))",
            "*" => "((U|T)A(A|G|R))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)(A|G|R))",
            "W" => "((U|T)G(A|G|R))",
            "X" => "...",
        );
    } elsif ($codontable == 4) {
        #-----------#
        # Mold, Protozoan, and Coelenterate Mitochondrial Code
        # and Mycoplasma/Spiroplasma Code (transl_table=4)
        #-----------#
        %p2c = (
            "B" => "((A(U|T).)|((U|T)(U|T)(A|G|R))|(C(U|T)G)|(G(U|T)G))",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "((U|T)A(A|G|R))",
            "*" => "((U|T)A(A|G|R))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)G(A|G|R))",
            "X" => "...",
        );
    } elsif ($codontable == 5) {
        #-----------#
        # Invertebrate Mitochondrial Code (transl_table=5)
        #-----------#
        %p2c = (
            "B" => "((A(U|T).)|((U|T|A|G|R)(U|T)G))",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "(CG.)",
            "S" => "(((U|T)C.)|(AG.))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y))",
            "_" => "((U|T)A(A|G|R))",
            "*" => "((U|T)A(A|G|R))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)(A|G|R))",
            "W" => "((U|T)G(A|G|R))",
            "X" => "...",
        );
    } elsif ($codontable == 6) {
        #-----------#
        # Ciliate, Dasycladacean and Hexamita Nuclear Code (transl_table=6)
        #-----------#
        %p2c = (
            "B" => "(A(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "((T|U)GA)",
            "*" => "((T|U)GA)",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "((CA(A|G|R))|((U|T)A(A|G|R)))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    } elsif ($codontable == 9) {
        #-----------#
        # Echinoderm and Flatworm Mitochondrial Code (transl_table=9)
        #-----------#
        %p2c = (
            "B" => "((A|G|R)(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "(CG.)",
            "S" => "(((U|T)C.)|(AG.))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "((U|T)A(A|G|R))",
            "*" => "((U|T)A(A|G|R))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AAG)",
            "N" => "(AA(U|T|C|Y|A))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)G(A|G|R))",
            "X" => "...",
        );
    } elsif ($codontable == 10) {
        #-----------#
        # Euplotid Nuclear Code (transl_table=10)
        #-----------#
        %p2c = (
            "B" => "(A(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "(((U|T)A(A|G|R))|((T|U)GA))",
            "*" => "(((U|T)A(A|G|R))|((T|U)GA))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    } elsif ($codontable == 11) {
        #-----------#
        # Bacterial, Archaeal and Plant Plastid Code (transl_table=11)
        #-----------#
        %p2c = (
            "B" => "((A(U|T)G)|(.(U|T)G))",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "(((U|T)A(A|G|R))|((T|U)GA))",
            "*" => "(((U|T)A(A|G|R))|((T|U)GA))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    } elsif ($codontable == 12) {
        #-----------#
        # Alternative Yeast Nuclear Code (transl_table=12)
        #-----------#
        %p2c = (
            "B" => "((A|C)(U|T)G)",
            "L" => "((C(U|T)(U|T|C|Y|A))|((U|T)(U|T)(A|G|R)))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y))|(C(U|T)G))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "(((U|T)A(A|G|R))|((T|U)GA))",
            "*" => "(((U|T)A(A|G|R))|((T|U)GA))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    } elsif ($codontable == 13) {
        #-----------#
        # Ascidian Mitochondrial Code (transl_table=13)
        #-----------#
        %p2c = (
            "B" => "((T|U|A|G|R)(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "(CG.)",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "((GG.)|(AG(A|G|R)))",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y))",
            "_" => "((U|T)A(A|G|R))",
            "*" => "((U|T)A(A|G|R))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)(A|G|R))",
            "W" => "((U|T)G(A|G|R))",
            "X" => "...",
        );
    } elsif ($codontable == 14) {
        #-----------#
        # Alternative Flatworm Mitochondrial Code (transl_table=14)
        #-----------#
        %p2c = (
            "B" => "(A(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "(CG.)",
            "S" => "(((U|T)C.)|(AG.))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "((U|T)AG)",
            "*" => "((U|T)AG)",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AAG)",
            "N" => "(AA(U|T|C|Y|A))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y|A))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)G(G|A|R))",
            "X" => "...",
        );
    } elsif ($codontable == 15) {
        #-----------#
        # Blepharisma Nuclear Code (transl_table=15)
        #-----------#
        %p2c = (
            "B" => "(A(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "(((U|T)AA)|((T|U)GA))",
            "*" => "(((U|T)AA)|((T|U)GA))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "((CA(A|G|R))|((U|T)AG))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    } elsif ($codontable == 16) {
        #-----------#
        # Chlorophycean MitochondrialCode (transl_table=16)
        #-----------#
        %p2c = (
            "B" => "(A(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R))|((U|T)AG))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "(((U|T)AA)|((T|U)GA))",
            "*" => "(((U|T)AA)|((T|U)GA))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    } elsif ($codontable == 21) {
        #-----------#
        # Trematode Mitochondrial Code (transl_table=21)
        #-----------#
        %p2c = (
            "B" => "((A|G|R)(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "(CG.)",
            "S" => "(((U|T)C.)|(AG.))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y))",
            "_" => "((U|T)A(A|G|R))",
            "*" => "((U|T)A(A|G|R))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AAG)",
            "N" => "(AA(U|T|C|Y|A))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)(A|G|R))",
            "W" => "((U|T)G(A|G|R))",
            "X" => "...",
        );
    } elsif ($codontable == 22) {
        #-----------#
        # Scenedesmus obliquus mitochondrial Code (transl_table=22)
        #-----------#
        %p2c = (
            "B" => "(A(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R))|((T|U)AG))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C(U|T|C|Y|G))|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "((U|T)(C|A|G|R)A)",
            "*" => "((U|T)(C|A|G|R)A)",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    } elsif ($codontable == 23) {
        #-----------#
        # Thraustochytrium Mitochondrial Code (transl_table=23)
        #-----------#
        %p2c = (
            "B" => "(((A|G|R)(U|T)G)|(A(U|T)(U|T)))",
            "L" => "((C(U|T).)|((U|T)(U|T)G))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "(((U|T)A(A|G|R))|((T|U)GA)|((T|U)(T|U)A))",
            "*" => "(((U|T)A(A|G|R))|((T|U)GA)|((T|U)(T|U)A))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    }


    #---------------------------------------------------------------#
    # make codon sequence, $qcodon,  with all possible combinations
    #---------------------------------------------------------------#

    $peplen = length($pep);
    undef($qcodon);
    foreach $i (0..$peplen - 1) {
        $peppos = $i + 1;
        $tmpaa = substr($pep, $i, 1);
        if ($tmpaa =~ /[ACDEFGHIKLMNPQRSTVWY_\*XU]/) {
            if ($qcodon !~ /\w/ && substr($pep, $i, 1) eq "M") {

                #---------------#
                # the first Met
                #---------------#

                $qcodon .= $p2c{"B"};
            } else {
                $qcodon .= $p2c{substr($pep, $i, 1)};
            }
        } elsif ($tmpaa =~ /\d/) {
            $qcodon .= "." x $tmpaa;
        } elsif ($tmpaa =~ /[-\.]/) {
            # nothing to do
        } else {
            $message = "pepAlnPos $peppos: $tmpaa unknown AA type. Taken as 'X'";
            push(@{$retval{'message'}}, $message);
            $qcodon .= $p2c{'X'};
        }
    }
    # print "$qcodon\n";


    #-----------------------------#
    # does $nuc contain $qcodon ?
    #-----------------------------#

    if ($nuc =~ /$qcodon/i) {
        $codon = $&;

        $retval{'codonseq'} = $codon;
        $retval{'result'} = 1;

    } else {
        #-------------------#
        # make 10 aa anchor
        #-------------------#

#        undef(@{$retval{'message'}});

#        $modpep = $pep;
#        1 while $modpep =~ s/(.{10})(.{10,})/$1\n$2/;
#        @anchor = split(/\n/, $modpep);

        undef(@preanchor);
        undef($tmpanc);
        $nanc = 0;
        foreach $i (0..$peplen - 1) {
            $tmpaa = substr($pep, $i, 1);
            $tmpanc .= $tmpaa;
            ++$nanc if ($tmpaa !~ /-/);
            if ($nanc == 10 || $i == $peplen - 1) {
                push(@preanchor, $tmpanc);
                undef($tmpanc);
                $nanc = 0;
            }
        }
        undef(@anchor);
        $lastanchorlen = length($preanchor[-1]);
        if ($lastanchorlen < 10) {
            foreach $i (0..$#preanchor - 1) {
                if ($i < $#preanchor - 1) {
                    push(@anchor, $preanchor[$i]);
                } else {
                    push(@anchor, "$preanchor[$i]$preanchor[$i + 1]");
                }
            }
        } else {
            @anchor = @preanchor;
        }

        undef($wholecodon);
        foreach $i (0..$#anchor) {
            # print "    $anchor[$i]\n";
            $anclen = length($anchor[$i]);
            undef(@qcodon);
            undef(@fncodon);
            foreach $j (0..$anclen - 1) {
                $peppos = $i * 10 + $j + 1;
                $tmpaa = substr($anchor[$i], $j, 1);
                if ($tmpaa =~ /[ACDEFGHIKLMNPQRSTVWY_\*XU]/) {
                    if ($i == 0 && $qcodon[$i] !~ /\w/ && $tmpaa eq "M") {

                        #----------------------------------#
                        # the first Met can be AGT|GTG|CTG
                        #----------------------------------#

                        $qcodon[$i] .= "((A|C|G|R)TG)";
                    } else {
                        $qcodon[$i] .= $p2c{$tmpaa};
                    }
                    $fncodon[$i] .= $p2c{'X'};
                } elsif ($tmpaa =~ /\d/) {
                    $qcodon[$i] .= "." x $tmpaa;
                    $fncodon[$i] .= "." x $tmpaa;
                } elsif ($tmpaa =~ /[-\.]/) {
                    # nothing to do
                } else {
                    #del $message = "pepAlnPos $peppos: $tmpaa unknown AA type. Replaced by 'X'";
                    #del push(@{$retval{'message'}}, $message);
                    $qcodon[$i] .= $p2c{'X'};
                    $fncodon[$i] .= $p2c{'X'};
                }
            }
            if ($nuc =~ /$qcodon[$i]/i) {
                $wholecodon .= $qcodon[$i];
            } else {
                $wholecodon .= $fncodon[$i];
            }
        }

        if ($nuc =~ /$wholecodon/i) {
            $codon = $&;
            $codonpos = 0;
            $tmpnaa = 0;
            foreach $i (0..$peplen - 1) {
                $peppos = $i + 1;
                $tmpaa = substr($pep, $i, 1);
                undef($tmpcodon);
                if ($tmpaa !~ /\d/ && $tmpaa !~ /-/) {
                    ++$tmpnaa;
                    $tmpcodon = substr($codon, $codonpos, 3);
                    $codonpos += 3;
                    if ($tmpnaa == 1 && $tmpaa eq "M") {
                        if ($tmpcodon !~ /((A|C|G|R)TG)/i) {
                            $message = "pepAlnPos $peppos: $tmpaa does not correspond to $tmpcodon";
                            push(@{$retval{'message'}}, $message);
                        }
                    } elsif ($tmpcodon !~ /$p2c{$tmpaa}/i) {
                        $message = "pepAlnPos $peppos: $tmpaa does not correspond to $tmpcodon";
                        push(@{$retval{'message'}}, $message);
                    }
                } elsif ($tmpaa =~ /\d/i) {
                    $tmpcodon = substr($codon, $codonpos, $tmpaa);
                    $codonpos += $tmpaa;
                }
                # print "$tmpaa    $tmpcodon\n";
            }
            $codon;

            $retval{'codonseq'} = $codon;
            $retval{'result'} = 2;

        } else {

            $retval{'result'} = -1;

        }

    }

    return(%retval);
}


#-----------------------------------------------------------------------

sub my1while {
    local($seq, $wlen) = @_;

    local($seqlen, $tmppos, @frags);

    $seqlen = length($seq);
    $tmppos = 1;


    while ($tmppos <= $seqlen) {
        $tmpfrag = substr($seq, $tmppos - 1, $wlen);
        push(@frags, $tmpfrag);
        $tmppos += $wlen;
    }

    return(@frags);

}


#---------------------------------------------------------

sub showhelp {
print STDERR<<EOF;

pal2nal.pl  (v14)

Usage:  pal2nal.pl  pep.aln  nuc.fasta  [nuc.fasta...]  [options]


    pep.aln:    protein alignment either in CLUSTAL or FASTA format

    nuc.fasta:  DNA sequences (single multi-fasta or separated files)

    Options:  -h            Show help 

              -output (clustal|paml|fasta|codon)
                            Output format; default = clustal

              -blockonly    Show only user specified blocks
                            '#' under CLUSTAL alignment (see example)

              -nogap        remove columns with gaps and inframe stop codons

              -nomismatch   remove mismatched codons (mismatch between
                            pep and cDNA) from the output

              -codontable  N
                    1  Universal code (default)
                    2  Vertebrate mitochondrial code
                    3  Yeast mitochondrial code
                    4  Mold, Protozoan, and Coelenterate Mitochondrial code
                       and Mycoplasma/Spiroplasma code
                    5  Invertebrate mitochondrial
                    6  Ciliate, Dasycladacean and Hexamita nuclear code
                    9  Echinoderm and Flatworm mitochondrial code
                   10  Euplotid nuclear code
                   11  Bacterial, archaeal and plant plastid code
                   12  Alternative yeast nuclear code
                   13  Ascidian mitochondrial code
                   14  Alternative flatworm mitochondrial code
                   15  Blepharisma nuclear code
                   16  Chlorophycean mitochondrial code
                   21  Trematode mitochondrial code
                   22  Scenedesmus obliquus mitochondrial code
                   23  Thraustochytrium mitochondrial code


              -html         HTML output (only for the web server)

              -nostderr     No STDERR messages (only for the web server)


    - sequence order in pep.aln and nuc.fasta should be the same.

    - IDs in pep.aln are used in the output.

EOF
}
