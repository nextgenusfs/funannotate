RunIprScan 1.1.0


LICENSE

This application free for academic and commercial use.


AUTHOR

Michael Thon <mike@michaelrthon.com>
http://www.michaelrthon.com


PREREQUISITES

A Java runtime environment (JRE).
RunIprScan has been tested on Ubuntu Linux and Mac OSX.  It should also work on
Windows, but I haven't tested it. I assume that users already are familiar with
the unix/linux command line. If you want a graphical user interface, checkout
Geneious (http://www.geneious.com) and install the InterProScan plugin.


INTRODUCTION

RunIprScan is a command line utility for batch submission of protein sequences
to the InterProScan server at the European Bioinformatics Institute. EBI has a
few example scripts on their website that accomplish the same task.  The main
difference is that RunIprScan submits sequences batchwise, up to 25 at a time
(the maximum allowed by EBI) which improves throughput. I've used RunIprScan to
submit tens of thousands of proteins to the server, and have let it run for
several days at a time.

The results are saved in xml format (one xml file per sequence) in an output
directory and there is an option to create a separate tab delimited file of the
Gene Ontology terms that are mapped to the InterPro terms. I also plan to add an
exporter for the InterPro terms.

Keep in mind that when using RunIprScan, your sequences are sent to the web
server as EBI, and as such, you are subject to their terms of service.


INSTALLATION

Unzip the distribution file and place the resulting folder somewhere on your
computer. Copy the runiprscan shell script into a bin directory (for example
~/bin) and edit it so the path points to the RunIprScan.jar file. the shell
script won't work on Windows, but you can execute the jar file directly but you
can also execute the RunIprScan.jar file directly, with something like:
java -jar RunIprScan.jar.


USE

A typical command line looks like this:

  runiprscan -v -i 3seqs.fasta -m mike@michaelrthon.com -o out

where 3seqs.fasta is a file that contains your protein sequences in FASTA
format. The output directory (out in the example above) should be created before
running the application. Please provide a valid email address on the command
line - if your sequences produce errors on the server then EBI will
automatically send you error reports.

Don't run more than one instance of RunIprScan at a time.  Doing so will exceed
EBI's limit of 25 sequences at a time, and will likely get you banned from their
service.

Occasionally, a sequence will produce an error on the server and no xml file is
produced.  If this happens, you can re-run RunIprScan using the same input file
and output directory, and it will check the output directory for missing xml
files, and resubmit them to the server. Usually, they will be processed
successfully the second time. I usually run the same runpirscan command 2 to 3
times to make sure all the proteins get processed successfully. If a protein
sequence consistently produces an error, there is probably a problem with the
format of the sequence.

Scanning 10,000 proteins with RunIprScan can take 2-3 days, and the processing
time can vary a lot, depending on the server load at EBI.  However, unless you
have a big compute cluster, RunIprScan is probably faster than a local
installation of the InterProScan application.  I have run InterProScan locally
on an 8 cpu server and it was much slower than RunIprScan.

Included with RunIprScan are three python scripts for processing the xml files
and extracting information from them.

ipr2go.py  extract GO terms into a three column tab separated output

ipr_table.py Parse one or more directories of xml files and output a table listing
             the number of times each InterPro term is encountered

extract_ipr.py Given a directory of XML files and an InterPro term, this script
               will extract all protein sequences that are annoteted with the term.


CHANGELOG

v1.1
Removed the -p option and replaced it with a python script

v1.0.1
