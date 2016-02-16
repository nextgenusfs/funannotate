#!/usr/bin/perl

##########################################################################################
#	  This file is part of Proteinortho.
#	  (C) 2009/2010 Marcus Lechner
# 
#	  Proteinortho is free software; you can redistribute it and/or modify
#	  it under the terms of the GNU General Public License as published
#	  by the Free Software Foundation; either version 2, or (at your
#	  option) any later version.
#
#	  Proteinortho is distributed in the hope that it will be useful, but
#	  WITHOUT ANY WARRANTY; without even the implied warranty of
#	  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#	  General Public License for more details.
#
#	  You should have received a copy of the GNU General Public License
#	  along with Proteinortho; see the file COPYING.  If not, write to the
#	  Free Software Foundation, Inc., 59 Temple Place - Suite 330,
#	  Boston, MA 02111-1307, USA.	
##########################################################################################

##########################################################################################
# About
##########################################################################################
# Proteinortho
# input fasta files with proteins
# output matrix with orthologous proteins
# 
# @author Marcus Lechner
# @email lechner@staff.uni-marburg.de
# @company University of Maruburg
# @date 2014-07-07
#
##########################################################################################

##########################################################################################
# Imports
##########################################################################################
use strict;
use warnings "all";
use File::Basename;
use threads;
use threads::shared;
use Thread::Queue;

##########################################################################################
# Variables
##########################################################################################
our $version = "5.11";
our $step = 0;		# 0/1/2/3	-> do all / only apply step 1 / only apply step 2 / only apply step 3
our $verbose = 1;	# 0/1		-> don't / be verbose
our $debug = 0;		# 0/1		-> don't / show debug data
#our $reflexiv = 0;	# 0/1		-> check sets against themselves
our $synteny = 0;	# 0/1		-> Apply synteny algorithm
our $neighbourjoin = 0;	# 0/1		-> Merge neighbours
our $duplication = 2;	# 0-9 not 1	-> Repeats for duplication extension
our $cs = 3;	# int		-> cs-value
our $alpha = 0.5;	# Alpha value for ffadj_mcs.py
our $connectivity = 0.1;# Algebr. connectivity threshold
our $cpus   = 0;	# 0 = autodetect
our $evalue = "1e-05";
our $coverage = 50;	# Percent coverage threshold for two proteins
our $identity = 25;	# Percent identity threshold for two proteins
our $blastmode = "blastp+";
#our $tmpdir = "./";	# Dir for tmp-files
our $sim = 0.95;
our $report = 3;
our $startat = 0;
our $stopat = -1;
our $keep = 0;
our $force = 0;
our $selfblast = 0;
our $twilight = 0;
our $singles = 0;
our $clean = 0;
our $blastOptions = "";
our $nograph = 0;
our $desc = 0;

# Internal
our $blastversion = "unknown";	# Auto-detected blast version
our $blastpath = "";
our $makedb = "";		# makedb command
our $blast = "";		# blast command
our $jobque = Thread::Queue->new();	# Jobs todo
our $jobs_done:shared = 0;		# Counter
our $jobs_todo = 0;			# Sum of jobs
our $project = "myproject";		# Project name
our $graph_lock :shared;
our $syn_lock :shared;
our $all_jobs_submitted :shared = 0;
our $po_path = &get_po_path();		# Determine local path
our $last_i = -1;
our $last_j = -1;
our $run_id = "";
our %gene_counter;		# Holds the number of genes for each data file (for sorting)
our $threads_per_process :shared = 1;	# Number of subthreads for blast

##########################################################################################
# Parameters
##########################################################################################
our @files = ();
foreach my $option (@ARGV) {
  if ($option =~ m/^--?step=(0|1|2|3)$/)  		{ $step = $1;   }
  elsif ($option =~ m/^--?verbose$/) 			{ $verbose = 1;  }
  elsif ($option =~ m/^--?verbose=(0|1)$/) 		{ $verbose = $1;  }
  elsif ($option =~ m/^--?debug$/)     			{ $debug = 1;  }
  elsif ($option =~ m/^--?debug=(0|1)$/)     		{ $debug = $1;  }
  elsif ($option =~ m/^--?p=(.*)$/)     		{ $blastmode = $1; }
  elsif ($option =~ m/^--?e=(.*)$/)  			{ $evalue = $1; }
  elsif ($option =~ m/^--?cpus=(\d*)$/)    		{ $cpus = $1; }
  elsif ($option =~ m/^--?cpus=auto$/)	     		{ $cpus = 0; }
  elsif ($option =~ m/^--?alpha=([0-9\.]+)$/)    	{ $alpha = $1; }
  elsif ($option =~ m/^--?report=([0-9]+)$/)   		{ $report = $1; }  
  elsif ($option =~ m/^--?conn=([0-9\.]+)$/)   		{ $connectivity = $1; }
  elsif ($option =~ m/^--?cov=([0-9]+)$/)    		{ $coverage = $1; }
  elsif ($option =~ m/^--?blastpath=(.+)$/)  		{ $blastpath  = $1."/"; }
  elsif ($option =~ m/^--?identity=([0-9]+)$/)    	{ $identity = $1; }
  elsif ($option =~ m/^--?identity=twilight$/)   	{ $twilight = 1; }
  elsif ($option =~ m/^--?sim=([0-9\.]+)$/)    		{ $sim = $1; }
  elsif ($option =~ m/^--?startat=([0-9]+)$/)    	{ $startat = $1; }
  elsif ($option =~ m/^--?stopat=([0-9]+)$/)    	{ $stopat = $1; }
  elsif ($option =~ m/^--?selfblast$/) 			{ $selfblast = 1; }  
  elsif ($option =~ m/^--?selfblast=(0|1)$/) 		{ $selfblast = $1; }  
  elsif ($option =~ m/^--?singles$/) 			{ $singles = 1; }  
  elsif ($option =~ m/^--?singles=(0|1)$/) 		{ $singles = $1; }  
  elsif ($option =~ m/^--?poff$/) 			{ $synteny = 1; }  
  elsif ($option =~ m/^--?synteny$/) 			{ $synteny = 1; }  
  elsif ($option =~ m/^--?synteny=(0|1)$/) 		{ $synteny = $1; }
  elsif ($option =~ m/^--?dups=0$/)			{ $duplication = 0; } 
  elsif ($option =~ m/^--?dups=([1-8])$/)		{ $duplication = $1+1;}
  elsif ($option =~ m/^--?neighbourjoin$/)		{ $neighbourjoin = 1; }  
  elsif ($option =~ m/^--?neighbourjoin=(0|1)$/)	{ $neighbourjoin = $1; }
  elsif ($option =~ m/^--?cs=([0-9]+)$/)		{ $cs = $1; }  
  elsif ($option =~ m/^--?keep/)			{ $keep = 1; }  
  elsif ($option =~ m/^--?force/)			{ $force = 1; }  
  elsif ($option =~ m/^--?clean/)			{ $clean = 1; }  
  elsif ($option =~ m/^--?nograph/)			{ $nograph = 1; }  
  elsif ($option =~ m/^--?graph/)			{ $nograph = 0; }  
  elsif ($option =~ m/^--?desc/)			{ $desc = 1; }  
  elsif ($option =~ m/^--?project=(.*)$/)		{ $project = $1; }
  elsif ($option =~ m/^--?blastParameters=(.*)$/i)	{ $blastOptions = $1;}
  elsif ($option !~ /^-/)				{ push(@files,$option); }
  else  {&print_usage(); die "Invalid command line option: \'$option\'!\n"; }
}

##########################################################################################
# Check parameters
##########################################################################################
if ($startat != 0 || $stopat != -1) {
	if ($step != 2) {
		&Error("Parameters -startat and -stopat only work for step 2!");
	}
	$run_id = "_$startat-$stopat";
}

our $simgraph = "$project.blast-graph$run_id";		# Output file graph
our $syngraph = "$project.ffadj-graph$run_id";		# Output file synteny

our $rm_simgraph = "$project.removed_blast-graph";		# Output remove graph
our $rm_syngraph = "$project.removed_ffadj-graph";		# Output remove graph

our $csimgraph = "$project.proteinortho-graph";		# Output file graph
our $csyngraph = "$project.poff-graph";				# Output file synteny
our $simtable = "$project.proteinortho";			# Output file graph
our $syntable = "$project.poff";					# Output file synteny
our $desctable = "$project.descriptions";			# Output file seq descriptions

##########################################################################################
# Run
##########################################################################################
&print_header;		# Show Proteinortho Header
&auto_cpus;		# Set number of CPUs
&check_blast;		# Check blast version

# Always do
&check_files;		# Check files, count sequences
@files = ();
foreach my $file (sort { $gene_counter{$b} <=> $gene_counter{$a} } keys %gene_counter) {push(@files,$file);}	# Biggest first


# Step 1, check files and make indices
if ($step == 0 || $step == 1) {
	print STDERR "\n**Step 1**\n";
	&generate_indices;	# Generate index files for blast
	if ($desc) {
		&write_descriptions;	# Write sequence description file
	}
}

# Step 2, run blast and synteny algorithm
if ($step == 0 || $step == 2) {
	print STDERR "\n**Step 2**\n";
	&init_graph_output;	# Initiate Output file(s)
	&run_blast;		# Run blasts
}

# Step 3, spacial clustering
if ($step == 0 || $step == 3) {
	print STDERR "\n**Step 3**\n";
	&cluster;							# form clusters and write outputs
	if ($clean) {&clean;}						# remove blast indices
	print STDERR "\nAll finished.\n";
}


##########################################################################################
# Funktions
##########################################################################################
sub clean {
	print STDERR "Removing temporary files...\n";
	if ($blastmode =~ /blastp/) {
		foreach my $file (@files) {
			unlink("$file.phr");
			unlink("$file.pin");
			unlink("$file.psq");
		}
	}
	else {
		foreach my $file (@files) {
			unlink("$file.nhr");
			unlink("$file.nin");
			unlink("$file.nsq");
		}
	}

	if ($nograph) {
		unlink($simgraph);
		if ($synteny) {unlink($syngraph);}
	}
}

sub cluster {
	print STDERR "Clustering by similarity (Proteinortho mode)\n";
	system("$po_path/proteinortho5_clustering -verbose $verbose -conn $connectivity -rmgraph '$rm_simgraph' $simgraph* >'$simtable'");
	if ($singles) {
		print STDERR "Adding singles...\n";
		my $fastas = "'".join("' '",@files)."'";
		system("$po_path/proteinortho5_singletons.pl $fastas <'$simtable' >>'$simtable'");
	}
        print STDERR "-> Output written to $simtable\n";
	unless ($nograph) {
		print STDERR "Writing graph...\n";
#		system("$po_path/proteinortho5_clean_edges -e '$rm_simgraph' $simgraph* >'$csimgraph'");
		system("$po_path/proteinortho5_clean_edges2.pl '$rm_simgraph' $simgraph* >'$csimgraph'");
		unless ($keep) {unlink($rm_simgraph);}
		print STDERR "-> Output written to $csimgraph\n";
	}

	if ($synteny) {
		print STDERR "\nClustering by gene-order (POFF mode)\n";
		system("$po_path/proteinortho5_clustering -verbose $verbose -conn $connectivity -rmgraph '$rm_syngraph' $syngraph* >'$syntable'");
		if ($singles) {
		print STDERR "Adding singles...\n";
			my $fastas = "'".join("' '",@files)."'";
			system("$po_path/proteinortho5_singletons.pl $fastas <'$syntable' >>'$syntable'");
		}
		print STDERR "-> Output written to $syntable\n";
		unless ($nograph) {
			print STDERR "Writing graph...\n";
#			system("$po_path/proteinortho5_clean_edges -e '$rm_syngraph' $syngraph* >'$csyngraph'");
			system("$po_path/proteinortho5_clean_edges2.pl '$rm_syngraph' $syngraph* >'$csyngraph'");
			unless ($keep) {unlink($rm_syngraph);}
			print STDERR "-> Output written to $csyngraph\n";
		}
	}
}

sub print_header {
	print STDERR "Proteinortho with PoFF version $version - An orthology detection tool\n",
	             "******************************************************************\n";
}

sub print_usage {
print STDERR << "JUS";

Usage: proteinortho5.pl [OPTIONS] FASTA1 FASTA2 [FASTA...]
Options: -e=          E-value for blast [default: 1e-05]
         -p=          blast program {blastn|blastp|blastn+|blastp+}
                      [default: blastp+]
         -project=    prefix for all result file names [default: myproject]
         -synteny     activate PoFF extension to separate similar sequences
                      by contextual adjacencies (requires .gff for each .fasta)
         -dups=       PoFF: number of reiterations for adjacencies heuristic,
                      to determine duplicated regions (default: 0)
         -cs=         PoFF: Size of a maximum common substring (MCS) for
                      adjacency matches (default: 3)
         -alpha=      PoFF: weight of adjacencies vs. sequence similarity
                      (default: 0.5)
         -desc        write description files (for NCBI FASTA input only)
         -keep        stores temporary blast results for reuse
         -force       forces recalculation of blast results in any case
         -cpus=       number of processors to use [default: auto]
         -selfblast   apply selfblast, detects paralogs without orthologs
         -singles     report singleton genes without any hit
         -identity=   min. percent identity of best blast hits [default: 25]
         -cov=        min. coverage of best blast alignments in % [default: 50]
         -conn=       min. algebraic connectivity [default: 0.1]
         -sim=        min. similarity for additional hits (0..1) [default: 0.95]
         -step=       1 -> generate indices
                      2 -> run blast (and ff-adj, if -synteny is set)
                      3 -> clustering
                      0 -> all (default)
         -blastpath=  path to your local blast (if not installed globally)
         -verbose     keeps you informed about the progress
         -clean       remove all unnecessary files after processing
         -graph       generate .graph files (pairwise orthology relations)
         -debug       gives detailed information for bug tracking

         More specific blast parameters can be defined by
         -blastParameters='[parameters]' (e.g. -blastParameters='-seg no')

         In case jobs should be distributed onto several machines, use
         -startat=    File number to start with (default: 0)
         -stopat=     File number to end with (default: -1)
 
JUS
}

sub init_graph_output {
#	if (-e $graph) {
#		&Error("Graph output file '$graph' already exists.");
#	}
	open(GRAPH,">$simgraph") || die("Could not open graph '$simgraph': $!");
	print GRAPH "# file_a\tfile_b\n# a\tb\tevalue_ab\tbitscore_ab\tevalue_ba\tbitscore_ba\n";
	close(GRAPH);

	unless ($synteny) {return;}
#	if (-e $syn) {
#		&Error("Synteny Graph output file '$syn' already exists.");
#	}
	open(SYN,">$syngraph") || die("Could not open graph '$syngraph': $!");
	print SYN "# file_a\tfile_b\n# a\tb\tevalue_ab\tbitscore_ab\tevalue_ba\tbitscore_ba\tsame_strand\tsimscore\n";
	close(SYN);
}

sub set_threads_per_process {
	lock($jobs_done);
	my $willdo = ($jobs_todo-$jobs_done+$_[0]);

	if ($debug) {	
		print STDERR "\nTODO: $jobs_todo DONE: $jobs_done Running: $_[0] -> $willdo\n";
	}

	if ($willdo < 1) {return;}

	my $optimal = int($cpus/$willdo);
	lock($threads_per_process);
	if ($optimal > $threads_per_process) {
		$threads_per_process = $optimal;
		if ($debug) {
			print STDERR "\nBlast subthreads was set to $threads_per_process\n";
		}
	}
}

sub run_blast {
	# Jobs todo
	$jobs_todo = 0;
	for (my $k = $startat; $k < scalar(@files); $k++) {$jobs_todo += scalar(@files)-1-$k+$selfblast;if ($stopat > -1 && $stopat <= $k) {last;}}
	&set_threads_per_process(0);	# Check if we can apply more threads, nothing is running so far
	&print_blast_stats();

	# Spawn worker threads
	for (my $i = 0; $i < $cpus; $i++) {threads->create('workerthread');}

	# For each file against each other file
	for (my $i = $startat; $i < scalar(@files); $i++) {
		for (my $j = $i+1-$selfblast; $j < scalar(@files); $j++) {
			# Wait for queque to get empty (with some buffer)
			while ($jobque->pending() > 256 + 2*$cpus) {
				sleep(1);
			}
			# Syncronize with other processes
			$jobque->enqueue("$i $j");
		}
		# Start-Stop
		if ($stopat > -1 && $stopat <= $i) {last;}
	}

	# Tell all threads they can stop
	{lock($all_jobs_submitted); $all_jobs_submitted = 1;}

	# Wait until all jobs are done
	foreach (threads->list()) 	{$_->join();}
	&print_blast_stats();		print STDERR "\n";
	print STDERR "-> Output written to $simgraph\n";
}

sub workerthread {
	my $thread_id = threads->tid();
	my $temp_file = "$project-$run_id-$thread_id";

	# Clean up, just to be safe
	unlink("$temp_file.tmp");
	unlink("$temp_file.log");
	unlink("$temp_file.matching");

	while (1) {
		my ($i, $j);
		while (1) {
			# Fetch new jobs
			my $job = $jobque->dequeue_nb();
			# If there is nothing
			unless (defined($job)) {
				# Check if more jobs need to be done
	
				{
				lock($jobs_done);							# 	Productive
				if ($jobs_done >= $jobs_todo) {						# 	Productive
#				lock($all_jobs_submitted);						# 	DEBUGGING
#				if ($all_jobs_submitted) {						# 	DEBUGGING
					if ($debug) {print STDERR "Thread $thread_id\tis leaving\n";}
					return;
				}}
				# If so, wait
					if ($debug) {print STDERR "Thread $thread_id\tis seeking work ($jobs_done / $jobs_todo)\n";}
				sleep(1);			
			}
			else {
				# Parse job
				($i, $j) = split(" ",$job);
				# Break the fetch loop
				last;
			}
		}

		my $file_i = $files[$i];
		my $file_j = $files[$j];
		my $short_file_i = $file_i;	$short_file_i =~ s/^.*\///;
		my $short_file_j = $file_j;	$short_file_j =~ s/^.*\///;

		# Work
		&set_threads_per_process(scalar(threads->list()));
		my $result_ij = &blast($file_i,$file_j);
		
		my $result_ji;
		if ($file_i eq $file_j) {
			# One run is enough (selfblast)
			$result_ji = $result_ij;
		}
		else {
			$result_ji = &blast($file_j,$file_i);
		}

		my %lengths;
		if ($file_i eq $file_j) {
			# Get lengths once is enough (selfblast)
			%lengths = %{&get_gene_lengths($file_i)};
		}
		else {
			%lengths = %{&get_gene_lengths($file_i,$file_j)};
		}
		my %reciprocal_matches = %{&match(\%lengths,$result_ij,$result_ji)};

		# Remove secondary hits if better exist	(test here instead of later)
		%reciprocal_matches = %{&adaptive_best_blast_matches(\%reciprocal_matches)};

		if ($synteny) {
			my ($ordered_matches, $track_pointer, $close_copies_pointer) = &synteny_matches(\%reciprocal_matches,$file_i,$file_j);
			open(PREGRAPH,">>$temp_file.tmp") || die("Could not open temp file '$temp_file.tmp': $!");
			print PREGRAPH $ordered_matches;
			close(PREGRAPH);
			my $cmd = "$po_path/ffadj_mcs.py $temp_file.tmp $alpha";
			if ($duplication) {
				$cmd .= " --repeat-matching $duplication --min-cs-size $cs";
			}
			if ($debug) {print STDERR "$cmd\n";}
			my $synt_stats = qx($cmd);
			chomp $synt_stats;
			$synt_stats =~ s/#.+\n//;
			
			# Reverse mapping of full gene ids, two seperate maps in case of overlapps in short ids
			my %full_id_map_i;				
			my %full_id_map_j;
			foreach (sort keys %reciprocal_matches) {
				my ($a, $b) = split(" ",$_);
				my $aa = &convertNCBI($a);
				my $bb = &convertNCBI($b);
				if ($aa ne $a) {
					if ($debug) {print STDERR "j_map: $aa -> $a\n";}
					$full_id_map_j{$aa} = $a;}
				if ($bb ne $b) {
					if ($debug) {print STDERR "i_map: $bb -> $b\n";}
					$full_id_map_i{$bb} = $b;}
			}

			# Reverse mapping of gene position to short id
			my %track = %{$track_pointer};
			my %close = %{$close_copies_pointer};
			# Generate hash for synteny hits
			my %synteny;
			unless (-s "$temp_file.matching") {
				print STDERR "Error: Failed to run $po_path/ffadj_mcs.py for\n$file_i vs $file_j\nMoving source to $temp_file.err for debugging\nI will continue, but results may be insufficient.\n";
				system("mv $temp_file.tmp $temp_file.err");
				next;
			}
			open(OSYNGRAPH,"<$temp_file.matching") || die("Could not open temp file $temp_file.matching: $!'");
			while(<OSYNGRAPH>) {
					chomp;
					my ($i, $j, $score) = split(/\s+/,$_,3);
					if (!defined($score) || $i =~ /[^0-9]/ || $i == 0 || length($i) > 10) {next;}
					unless (defined($track{$file_i.$i})) {
						print STDERR "Could not find i: ".$file_i.$i."\n";	next;
					}
					unless (defined($track{$file_j.$j})) {
						print STDERR "Could not find j: ".$file_j.$j."\n";	next;
					}
					# Remap to full ID
					my $a = $track{$file_i.$i};
					if (defined($full_id_map_i{$a})) {$a = $full_id_map_i{$a};}
					my $b = $track{$file_j.$j};
					if (defined($full_id_map_j{$b})) {$b = $full_id_map_j{$b};}
					# Store
					$synteny{"$b $a"} = $score;

					# Close copies
					if ($neighbourjoin && defined($close{$i})) {
						my @partners = split(',',$close{$i});
						foreach (@partners) {
							my $c = $track{$file_i.$_};
							if (defined($full_id_map_i{$c})) {$c = $full_id_map_i{$c};}
							# Store
							$synteny{"$b $c"} = $score;
							if ($debug) {print STDERR "Storing addional proximity edge $a & $b -> $c\n";}
						}
					}
				}
				close(OSYNGRAPH);
				unlink("$temp_file.tmp");
				unlink("$temp_file.log");
				unlink("$temp_file.matching");

			{
			lock($syn_lock);
			
			open(SYN,">>$syngraph") || die("Could not open file '$syngraph': $!");
			print SYN "# $short_file_j\t$short_file_i\n";
			print SYN "# Scores: $synt_stats\n";
			foreach (sort keys %reciprocal_matches) {
				if (!defined($synteny{$_})) {if ($debug) {print STDERR "FAIL: $_\n";} next;}		# Not reported by gene-order algo
				my $line = "$_ ".$reciprocal_matches{$_}." $synteny{$_}";
				$line =~ s/ /\t/g;
				print SYN "$line\n";
			}
			close(SYN);
			}

		}

		{
			lock($graph_lock);
			open(GRAPH,">>$simgraph") || die("Could not open file '$simgraph': $!");
			print GRAPH "# $short_file_j\t$short_file_i\n";
			foreach (sort keys %reciprocal_matches) {
				my $line = "$_ ".$reciprocal_matches{$_};
				$line =~ s/ /\t/g;
				print GRAPH "$line\n";
			}
			close(GRAPH);
		}
		# Count
		{
			lock($jobs_done);
			$jobs_done++;
		}
		&print_blast_stats();	# Needs jobs_done to be free
	}
}

sub identitybylength {
        # Accoding to the formula of Rost, 1999 (Twilight-zone)
	if($_[0] <= 11)		{return 100;}
	if($_[0] <= 450)	{return 480*$_[0]**(-0.32*(1+exp(-$_[0]/1000)));}
	return 19.5;
}

sub adaptive_best_blast_matches {
	my %reciprocal_matches = %{(shift)};

	if ($debug) {
		print STDERR "\nStart with ";
		print STDERR scalar(keys %reciprocal_matches);
		print STDERR " edges\n";
	}

	my %best;
	my %best_gen;
	# Gather best
	foreach (keys %reciprocal_matches) {
		my ($a,$b) = split(" ",$_);
		my ($evalue_ab,$bitscore_ab,$evalue_ba,$bitscore_ba) = split(" ",$reciprocal_matches{$_});
		if (!defined($best{$a}) || $best{$a} < $bitscore_ab) {
			$best{$a} = $bitscore_ab;
			$best_gen{$a} = $b;
		}
		if (!defined($best{$b}) || $best{$b} < $bitscore_ba) {
			$best{$b} = $bitscore_ba;
			$best_gen{$b} = $a;
		}
	}

	if ($debug) {
		my %gene_num;
		# Count gene number
		foreach (keys %reciprocal_matches) {
			my ($a,$b) = split(" ",$_);
			$gene_num{$a}++;
			$gene_num{$b}++;
		}
		print STDERR "Number of genes: ".scalar(keys %gene_num)."\n";

	#	foreach (keys %best) {
	#		print STDERR "Best score for $_ is $best{$_} ($best_gen{$_})\n";
	#	}
	}

	# Remove worse
	foreach (keys %reciprocal_matches) {
		my ($a,$b) = split(" ",$_);
		my ($evalue_ab,$bitscore_ab,$evalue_ba,$bitscore_ba) = split(" ",$reciprocal_matches{$_});
		if 	($best{$a}*$sim > $bitscore_ab) {delete $reciprocal_matches{$_}; if ($debug) {my $v = $bitscore_ab/$best{$a}; print STDERR "Removed $_ because:\t$best{$a} vs $bitscore_ab ($v)\n";}}
		elsif 	($best{$b}*$sim > $bitscore_ba) {delete $reciprocal_matches{$_}; if ($debug) {my $v = $bitscore_ba/$best{$b}; print STDERR "Removed $_ because:\t$best{$b} vs $bitscore_ba ($v)\n";}}
	}

	if ($debug) {
	print STDERR "\nEnd with ";
	print STDERR scalar(keys %reciprocal_matches);
	print STDERR " edges\n";

	my %gene_num;
	# Count gene number
	foreach (keys %reciprocal_matches) {
		my ($a,$b) = split(" ",$_);
		$gene_num{$a}++;
		$gene_num{$b}++;
	}
	print STDERR "Number of genes: ".scalar(keys %gene_num)."\n";}

	return \%reciprocal_matches;
}

sub synteny_matches {
	my %reciprocal_matches = %{(shift)};
	my $file_i = shift; 
	my $file_j = shift;

	# Get order for both species (same hash as ids are non overlapping)
	my %order;
	my %track;
	for my $file ($file_i, $file_j) {
		# Get Coordinates for all genes
		my %coords = %{&read_details($file)};
		my $counter = 0;
		# Number them according to their order
		foreach my $id (sort 
			{
				my @a = split("\t",$coords{$a});
				my @b = split("\t",$coords{$b});
	
#				#chr strand pos
#				if ($a[0] ne $b[0]) {return $a[0] cmp $b[0];}
#				if ($a[1] ne $b[1]) {return $a[1] cmp $b[1];}
#				return $a[2] <=> $b[2];

				#chr pos
				if ($a[0] ne $b[0]) {return $a[0] cmp $b[0];}
				return $a[2] <=> $b[2];
			}	(keys(%coords))) {
			my @v = split("\t",$coords{$id});
			$order{$id} = ++$counter."\t".$v[1];	# Store strand info
			$track{$file.$counter} = $id;		# Reverse Mapping
		}
	}

	my $output = "";

	my @multis;	# array that contains all multi-edges
	# Convert reciprocal matches to ffadj input
	foreach (keys %reciprocal_matches) {
		my @values = split(" ",$reciprocal_matches{$_});
		my ($a, $b) = split(" ",$_);
		unless (defined($order{$a})) {$a = &convertNCBI($a);}
		unless (defined($order{$b})) {$b = &convertNCBI($b);}
		my @a = split(" ",$order{$a});
		my @b = split(" ",$order{$b});
		unless (defined($a[0])) {die();}

		unless (defined($multis[$a[0]])) {$multis[$a[0]] = $b[0];}
		else 				 {$multis[$a[0]] .= ','.$b[0];}

		$output .= "$b[0]\t$a[0]\t";			# Positions
		if ($a[1] eq $b[1]) 	{$output .= "1\t";}	# Same strand?
		else 			{$output .= "-1\t";}
		my $score = (&edgeweight($values[0])+&edgeweight($values[2]))/2;	# Score made from e-values
		$output .= $score."\n";	
	}

	# Check multis
	my %close_copies;
	if ($neighbourjoin) {
		for (my $i = 1; $i < scalar(@multis); $i++) {
			unless (defined($multis[$i])) {next;}
			my @partners = sort { $a <=> $b } split(',',$multis[$i]);
			if (scalar(@partners) <= 1) {next;}
			my $dist_limit = 2;	# How far can tandem copies be away from each other? (0/1 = off, 2 = immediate, ...
			my $last = 999999999999999;
			foreach my $new (@partners) {
				if (abs($last-$new) < $dist_limit) {
					if (!defined($close_copies{$last}))  	{$close_copies{$last} = $new;}
					else 					{$close_copies{$last} .= ','.$new;}
					$close_copies{$new} = $last;		# The list is sortet, so we are here for the frist time
				}
				$last = $new;
			}
		}
	}

	return ($output, \%track, \%close_copies);
}

# Read the length of all genes in a given set of fasta files
sub get_gene_lengths {
	my %lengths;
	while (my $file = shift @_) {
		my $last_gene = "";
		open(FASTA,"<$file") || &Error("Could not open '$file': $!");
		while (<FASTA>) {
			chomp;
			if ($_ =~ />/) {
				$_ =~ s/^>//;	$_ =~ s/\s.*//;
				$lengths{$_} = 0;
				$last_gene = $_;
			}
			else {
				$lengths{$last_gene} += length($_);
			}
		}
		close(FASTA);
	}

	return \%lengths;
}

sub print_blast_stats {
	if (!$verbose) {return;}
	{
		if ($jobs_todo == 0) {die("Nothing to do. This should not happen!");}
		lock($jobs_done);
		my $percent = int($jobs_done/$jobs_todo*10000)/100;
		print STDERR "\r                                                                               ";		
		print STDERR "\rRunning blast analysis: $percent% ($jobs_done/$jobs_todo)";
	}
}

sub match {
	my %length = %{(shift)};
	my @i = @{(shift)};
	my @j = @{(shift)};
	
	my %legal_i = &get_legal_matches(\%length,@i);
	my %legal_j = &get_legal_matches(\%length,@j);
	
	return &get_reciprocal_matches(\%legal_i,\%legal_j);
}

sub get_reciprocal_matches {
	my %i = %{(shift)};
	my %j = %{(shift)}; 
	my %reciprocal;

	foreach (keys %i) {
		my ($i, $j) = split(" ",$_);

		# Reciprocal hit?
		if (!defined($j{"$j $i"})) {next;}

		# Merge both
		$reciprocal{$_} = $i{$_}." ".$j{"$j $i"};	# evalue_ij, bitscore_ij, evalue_ji, bitscore_ji
	}

	return \%reciprocal;
}

sub get_legal_matches {
	my %length = %{(shift)};

	my %result;
	foreach (@_) {
		my ($query_id,$subject_id,$local_identity,$alignment_length,$mismatches,$openings,$query_start,$query_end,$subject_start,$subject_end,$local_evalue,$local_bitscore) = split(/\s+/);
	
		# Bug tracking
		unless (defined($length{$query_id})) 	{print STDERR "ERROR: Query gene ID '$query_id' is present in blast output but was not present in FASTA input! Skipping line.\n"; next;}
		unless (defined($length{$subject_id})) 	{print STDERR "ERROR: Subject gene ID '$subject_id' is present in blast output but was not present in FASTA input! Skipping line.\n"; next;}

		## Check for criteria
		# Well formatted
		if (!defined($local_bitscore)) 							{next;}
		# Percent identity
		if (!$twilight && $local_identity < $identity) 					{next;}
		if ( $twilight && $local_identity < &identitybylength($alignment_length))	{next;} 
		# Min. length
		if ($alignment_length < $length{$query_id}*($coverage/100)+0.5) 		{next;}
		if ($alignment_length < $length{$subject_id}*($coverage/100)+0.5) 		{next;}
		# It hit itself (only during selfblast)
		if ($query_id eq $subject_id) 							{next;}

		# Similar hits? Take the better one
		if (defined($result{"$query_id $subject_id"})) {
			my ($remote_evalue, $remote_bitscore) = split(" ",$result{"$query_id $subject_id"});
			if ($local_evalue > $remote_evalue) {next;}			
			if ($local_bitscore < $remote_bitscore) {next;}
		}

		# Store data
		$result{"$query_id $subject_id"} = "$local_evalue $local_bitscore";
	}

	return %result;
}

# Auto set the number of CPUs
sub auto_cpus {
	if ($cpus == 0) {
		my $cpu_x = 0;
		# Linux
		if (-e "/proc/cpuinfo") {
			$cpu_x = qx(grep processor /proc/cpuinfo | wc -l);
		}
		# Try Mac
		else {
			$cpu_x = qx(system_profiler | grep CPUs:);
		}
		$cpu_x =~ s/[^0-9]//g;
		if (length($cpu_x) == 0 || $cpu_x == 0) {print STDERR "failed! Use 1 core only\n";$cpu_x = 1;}
		print STDERR "Detected $cpu_x available CPU threads, ";
		$cpus = $cpu_x;
	}
	else {
		print STDERR "Using $cpus CPU threads, ";
	}
}

sub generate_indices {
	print STDERR "Generating indices\n";
	foreach my $file (@files) {
		if ($file =~ /\s/) {die("File name '$file' contains whitespaces. This might lead to undesired effects. Please change file name!\n");}
		if ($verbose) {print STDERR "Building database for '$file'\t(".$gene_counter{$file}." sequences)\n";}
		system("$makedb '$file' >/dev/null");
	}
	unlink('formatdb.log');
}

sub blast {
	my $command = "";
	if 	($blastmode eq "blastp" || $blastmode eq "blastn") 	{lock($threads_per_process); $command = $blastpath."blastall -a $threads_per_process -d '$_[0]' -i '$_[1]' -p $blastmode -m8 -e $evalue $blastOptions";}
	elsif	($blastmode eq "blastp+") 				{lock($threads_per_process); $command = $blastpath."blastp -num_threads $threads_per_process -db '$_[0]' -query '$_[1]' -evalue $evalue -outfmt 6 $blastOptions";}
	elsif	($blastmode eq "blastn+") 				{lock($threads_per_process); $command = $blastpath."blastn -num_threads $threads_per_process -db '$_[0]' -query '$_[1]' -evalue $evalue -outfmt 6 $blastOptions";}
	else	{die("This should not happen!");}

	my $a = $_[0];
	my $b = $_[1];
	$a =~ s/^.*\///;
	$b =~ s/^.*\///;
	my $bla = "$a.vs.$b.bla";

	# File does not exists yet or I am forced to rewrite it
	if (!(-s $bla) || $force) {
		if ($debug) {print STDERR "$command >$bla\n";}
		system("$command | sort -k11,11g >'$bla'");		# run blast and presort (speeds up best alignment search but is NOT mandatory)
		my @data = &readFile($bla);
		unless ($keep) {unlink($bla);}				# delete tmp file
		return \@data;
	}
	# Otherwise, use existing data
	if ($verbose) {print STDERR "\nNote: '$bla' exists, using pre-calculated data\n";}
	my @data = &readFile($bla);
	return \@data;
}

sub readFile {
	open(FILE,"<$_[0]") || die("Error, could not open file $_[0]: $!");
	my @data = <FILE>;
	close(FILE);
	chomp @data;
	return @data;
}

sub check_blast {
	if ($blastmode eq "blastp+" || $blastmode eq "blastn+") {
		my $tmp = $blastmode;
		$tmp =~ s/\+//g;
		my $cmd = $blastpath."$tmp -h";
		my $out = qx($cmd);
		if ($out =~ /DESCRIPTION.*?\n\s*(.+)\n/) {
			my @version = split(/\s+/,$1);
			$blastversion = pop @version;

			# Commands
			if 	($blastmode eq "blastp+") {$makedb = $blastpath."makeblastdb -dbtype prot -in";}
			elsif 	($blastmode eq "blastn+") {$makedb = $blastpath."makeblastdb -dbtype nucl -in";}
			else	{die("This should not happen!");}
  
			print STDERR "Detected NCBI BLAST version $blastversion\n";
			return;
		}
		&Error("Failed to detect '$blastmode'! Tried to call '$tmp'.");
	}
	elsif ($blastmode eq "blastp" || $blastmode eq "blastn") {
		my $cmd = $blastpath."blastall";
		my @blastv = qx($cmd);
		foreach (@blastv) {
			chomp;
			if ($_ =~ /blastall.+?([^\s]+)/) {
				$blastversion = $1;
				if 	($blastmode eq "blastp") {$makedb = $blastpath."formatdb -p T -o F -i";}
				elsif 	($blastmode eq "blastn") {$makedb = $blastpath."formatdb -p F -o F -i";}
				else	{die("This should not happen!");}

				print STDERR "Detected NCBI BLAST version $blastversion\n";
				return;
			}
		}
		&Error("Failed to detect '$blastmode'! Tried to call 'blastall'.");
	}

	&Error("Blast mode '$blastmode' is not supported. Feel free to ask the author to add it.");
}

# Check plausibility of files
sub check_files {
	my %ids;
	if (scalar(@files) < 2 && $step != 3)		{&print_usage; &Error("I need at least two files to compare something!");}
	print STDERR "Checking input files\n";
	foreach my $file (@files) {
		&read_details($file,\%ids);
	}
}

sub convertNCBI {
	my $long_id = shift;
	$long_id =~ s/\|$//g;					
	my @tmp = split(/\|/,$long_id);	# take the last column for NCBI format like patterns (e.g. gi|158333234|ref|YP_001514406.1|)
	return pop(@tmp);
}

sub read_details {
	my %ids;
	my %genes;
	my $file = shift;
	my $test = 0;
	if (defined($_[0])) {$test = 1; %ids = %{(shift)};}	# if no ID Hash is give, we do not want to test but to fetch the gff data

	if (!-e $file)		{&Error("File '$file' not found!");}
	open(FASTA,"<$file") || &Error("Could not open '$file': $!");
	while (<FASTA>) {
		if ($_ =~ />/) {
			$gene_counter{$file}++;
			chomp;
			$_ =~ s/^>//;
			$_ =~ s/\s.*//;
			if ($test) {
				if (defined($ids{$_}))	{&Error("Gene ID '$_' is defined at least twice:\n$ids{$_}\n$file");}
				$ids{$_} = $file;
			}
			if ($synteny) {
				my $short_id = &convertNCBI($_);
				$genes{$short_id} = 1;
			}
		}
	}
	close(FASTA);

	unless ($synteny) {return;}

	my %coordinates;
	if ($verbose && $test) {print STDERR "$file\t".scalar(keys %genes)." genes\n";}
	my $gff = &gff4fasta($file);
	open(GFF,"<$gff") || &Error("Could not open '$gff': $!");
	while (<GFF>) {
		if ($_ =~ /^#/) {next;}
		# e.g. NC_009925.1	RefSeq	CDS	9275	10096	.	-	0	ID=cds8;Name=YP_001514414.1;Parent=gene9;Dbxref=Genbank:YP_001514414.1,GeneID:5678848;gbkey=CDS;product=signal peptide peptidase SppA;protein_id=YP_001514414.1;transl_table=11
		my @col = split(/\t+/,$_);
		if ($col[2] ne "CDS") {next;}
		if ($col[8] =~ /Name=([^;]+)/i && defined($genes{$1})) {
			delete $genes{$1};
#			if (!$test) {$coordinates{$1} = "$col[0]\t$col[6]\t$col[3]";}	# store
			if (!$test && $col[6] eq "+") {$coordinates{$1} = "$col[0]\t$col[6]\t$col[3]";}	# store
			if (!$test && $col[6] eq "-") {$coordinates{$1} = "$col[0]\t$col[6]\t$col[4]";}	# store
		}
		elsif ($col[8] =~ /ID=([^;]+)/i && defined($genes{$1})) {
			delete $genes{$1};
#			if (!$test) {$coordinates{$1} = "$col[0]\t$col[6]\t$col[3]";}	# store
			if (!$test && $col[6] eq "+") {$coordinates{$1} = "$col[0]\t$col[6]\t$col[3]";}	# store
			if (!$test && $col[6] eq "-") {$coordinates{$1} = "$col[0]\t$col[6]\t$col[4]";}	# store
		}
	}
	close(GFF);

	if (scalar(keys %genes)) {
		my @tmp = keys %genes;
		&Error("No coordinate found for these gene(s): ".join(",",@tmp)."\nusing '$gff' and '$file'");
	}

	if (!$test) {return \%coordinates;}		# store
}

sub Error {
	print STDERR "Error: ".$_[0]."\n";
	exit 0;
}

# Remove .fasta/.faa etc. and change it to .gff
sub gff4fasta {
	my $gff = shift;
	$gff =~ s/\.[^.]+$/.gff/;
	return $gff;
}

sub get_po_path {
	my @tmppath = fileparse($0); # path to the C++-part of this program
	return $tmppath[1];
}

sub edgeweight {
	# 1e-10 = 0.15, 1e-20 = 0.3, 1e-40 = 0.6, 1e-66+ = 1.0
	if ($_[0] == 0) {return 1;}
	my $x = -1*&log10($_[0])/100*1.5;
	if ($x > 1) {return 1;}
	if ($x <= 0) {return 0.0001;}
	return $x;
}

sub log10 {
	return log($_[0])/log(10);
}

sub write_descriptions {
	print STDERR "Writing sequence descriptions\n";
	open DESC, '>', $desctable;
	foreach my $file (@files) {
		if ($verbose) {print STDERR "Extracting descriptions from '$file'\t(".$gene_counter{$file}." entries)\n";}
		open FASTA, '<', $file;
		while (<FASTA>) {
			chomp;
			if (m/^>(\S+)(\s+(.*))?$/) {
				print DESC $1, "\t", ($3 || "unannotated sequence"), "\n";
			}
		}
	}
	print STDERR "-> Output written to $desctable\n";

}
