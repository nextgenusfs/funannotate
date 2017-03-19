#!/bin/bash

#setup shell script for funannotate databases
RED='\033[0;91m'
NC='\033[0m'

#likely run some checks here
command -v hmmpress >/dev/null 2>&1 || { echo "Funannotate setup requires HMMer 3.1 but it's not in PATH.  Aborting." >&2; exit 1; }
command -v wget >/dev/null 2>&1 || { echo "Funannotate setup requires wget but it's not in PATH.  Aborting." >&2; exit 1; }
command -v makeblastdb >/dev/null 2>&1 || { echo "Funannotate setup requires BLAST+ but it's not in PATH.  Aborting." >&2; exit 1; }

#try to be smart about finding previous databases if a user updates via homebrew.
#parse input commands to script
for i in "$@"
do
case $i in
    -d=*|--database=*)
    outputdir="${i#*=}"
    shift # past argument=value
    ;;
    -m=*|--mode=*)
    mode="${i#*=}"
    shift # past argument=value
    ;;
    *)
            # unknown option
    ;;
esac
done
echo "OutputDir    = ${outputdir}"
echo "Script mode  = ${mode}"

dir=$(pwd)
if [ -z "$outputdir" ]; then
    #check if softlink is already present, if so get target
    if [ -e 'DB' ]; then
        echo "Found DB symlinked folder"
        if [[ "$OSTYPE" == "darwin"* ]]; then
            outputdir=$(readlink DB)
        else
            outputdir=$(readlink -f DB)
        fi
        echo -n "DB directory set to ($outputdir), continue [y/n]: "
        read question1
        if [ $question1 == 'n' ]; then
            echo -n "Enter path to DB directory: "
            read dbname
            outputdir=$dbname
        fi
    else
        #no softlink found, check if in libexec folder which means homebrew install, try to look for previous version
        if [[ $dir == *"libexec"* ]]; then
            echo "HomeBrew installation detected, looking for any previous versions"
            pre_vers=$(ls ../../ | sort | tail -2 | head -1)
            curr_vers=$(ls ../../ | sort | tail -1 | head -1)
            if [ "$pre_vers" == "$curr_vers" ]; then
                echo "This is the first HomeBrew install detected."
                outputdir='/usr/local/share/funannotate'
                echo -n "Default DB directory set to ($outputdir), continue [y/n]: "
                read question1
                if [ $question1 == 'n' ]; then
                    echo -n "Enter path to DB directory: "
                    read dbname
                    outputdir=$dbname
                fi
            else
                if [[ "$OSTYPE" == "darwin"* ]]; then
                    outputdir=$(readlink ../../$pre_vers/libexec/DB)
                else
                    outputdir=$(readlink -f ../../$pre_vers/libexec/DB)
                fi
                echo "Symlink found to $outputdir, setting up DB"
            fi    
        else
            echo "HomeBrew installation not detected, specify DB installation directory"
            outputdir="$HOME/funannotate/"
            echo -n "Default DB directory set to ($outputdir), continue [y/n]: "
            read question1
            if [ $question1 == 'n' ]; then
                echo -n "Enter path to DB directory: "
                read dbname
                outputdir=$dbname
            fi
        fi
    fi
fi


if [ "$mode" == 'all' ]; then
	db='pass'
	dep='pass'
elif [ "$mode" == 'db' ]; then
    db='pass'
    dep='fail'
elif [ "$mode" == 'dep' ]; then
    db='fail'
    dep='pass'
else
    echo -e "    To download databases and check dependencies:   ./setup.sh
    To just download databases:  ./setup.sh db
    To just check dependencies:  ./setup.sh dep"
	exit
fi

if [ "$db" = 'pass' ]; then
    #start downloading databases
    cd $dir
    echo "Creating folder to store DB files: $outputdir"
    mkdir -p $outputdir
    ln -s $outputdir DB
    cd $outputdir

    #Do MEROPS first as need to download manually, wait for download to apear before moving on
    echo "-----------------------------------------------"
    echo "Okay, starting downloading of databases...."
    echo "-----------------------------------------------"

    #check if Merops is already downloaded
    if [ ! -f merops_formatted.fa ]; then
        wget -c --tries=0 --read-timeout=20 --output-document=merops_scan.lib https://uwmadison.box.com/shared/static/fvx5wt5ghhv5991mjytbz1lbyigoex7h.lib
        #echo "You need to manually download the MEROPS protease database as it requires a log in"
        #echo "download: merops_scan.lib   from here: https://merops.sanger.ac.uk/download/"
        #echo "then move the file into /usr/local/share/funannotate, once the file is in the folder the script will proceed."
        #until [ -f merops_scan.lib ]
        #do
        #     sleep 5
        #done
        tr -d '\r' < merops_scan.lib | sed 's/ - /#/g' | while read line; do set -- "$line"; IFS="#"; declare -a Array=($*); if [[ "${Array[0]}" == ">"* ]]; then echo ${Array[0]} ${Array[2]}; else echo $line; fi; done > merops_formatted.fa
        makeblastdb -in merops_formatted.fa -input_type fasta -dbtype prot -title MEROPS -parse_seqids -out MEROPS
    else
        echo "MEROPS DB found, skipping download"
        echo "-----------------------------------------------"
    fi

    #get uniprot and format database
    if [ ! -f uniprot_sprot.fasta ]; then
        echo "Now downloading/formatting UniProt DB"
        wget -c --tries=0 --read-timeout=20 ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
        gunzip uniprot_sprot.fasta.gz
        makeblastdb -in uniprot_sprot.fasta -input_type fasta -dbtype prot -title uniprot -parse_seqids -out uniprot
        echo "-----------------------------------------------"
    else
        echo "UniProt DB found, skipping download. To update delete uniprot_sprot.fasta."
        echo "-----------------------------------------------"
    fi

    #get PFAM database and associated mapping file
    if [ ! -f Pfam-A.hmm ]; then
        echo "Now downloading/formatting PFam-A DB"
        wget -c --tries=0 --read-timeout=20 ftp://ftp.ebi.ac.uk/pub/databases/Pfam//current_release/Pfam-A.hmm.gz
        gunzip Pfam-A.hmm.gz
        hmmpress Pfam-A.hmm
        echo "-----------------------------------------------"
    else
        echo "Pfam-A DB found, skipping download. To update delete Pfam-A.hmm"
        echo "-----------------------------------------------"
    fi

    #get pFAM mapping tsv vile
    if [ ! -f Pfam-A.clans.tsv ]; then
        echo "Now downloading PFAM mapping file"
        wget -c --tries=0 --read-timeout=20 ftp://ftp.ebi.ac.uk/pub/databases/Pfam//current_release/Pfam-A.clans.tsv.gz
        gunzip Pfam-A.clans.tsv.gz
        echo "-----------------------------------------------"
    else
        echo "PFAM mapping found, skipping download. To update delete Pfam-A.clans.tsv"
        echo "-----------------------------------------------"
    fi

    #get dbCAN database
    if [ ! -f dbCAN.hmm ]; then
        echo "Now downloading/formatting dbCAN CAZyme DB"
        wget -c --tries=0 --read-timeout=20 http://csbl.bmb.uga.edu/dbCAN/download/dbCAN-fam-HMMs.txt
        wget -c --tries=0 --read-timeout=20 http://csbl.bmb.uga.edu/dbCAN/download/FamInfo.txt
        mv FamInfo.txt dbCAN.info.txt
        sed 's/\.hmm$//g' dbCAN-fam-HMMs.txt > dbCAN.hmm
        hmmpress dbCAN.hmm
        echo "-----------------------------------------------"
    else
        echo "dbCAN DB found, skipping download"
        echo "-----------------------------------------------"
    fi

    #download Eggnog
    #if [ ! -f fuNOG_4.5.hmm ]; then
    #    echo "Now downloading/formatting EggNog 4.5 DB"
    #    wget -c --tries=0 --read-timeout=20 http://eggnogdb.embl.de/download/eggnog_4.5/data/fuNOG/fuNOG.hmm.tar.gz
    #    wget -c --tries=0 --read-timeout=20 http://eggnogdb.embl.de/download/eggnog_4.5/data/fuNOG/fuNOG.annotations.tsv.gz
    #    gunzip fuNOG.annotations.tsv.gz
    #    tar -zxf fuNOG.hmm.tar.gz
    #    find fuNOG_hmm/ -maxdepth 1 -type f -name '*.hmm' -exec cat '{}' \; > fuNOG_4.5.hmm
    #    hmmpress fuNOG_4.5.hmm
    #    rm fuNOG.hmm.tar.gz
    #    rm -R fuNOG_hmm/
    #    echo "-----------------------------------------------"
    #else
    #    echo "EggNog 4.5 DB found, skipping download"
    #    echo "-----------------------------------------------"
    #fi
    
    #get BUSCO and fungi models
    rm -R fungi
    if [ ! -d fungi ]; then
        echo "Downloading BUSCO fungi models"
        wget -c --tries=0 --read-timeout=20 http://busco.ezlab.org/v2/datasets/fungi_odb9.tar.gz
        tar -zxf fungi_odb9.tar.gz
        mv fungi_odb9 fungi
        echo "-----------------------------------------------"
    else
        echo "BUSCO fungi DB found, skipping download"
        echo "-----------------------------------------------"
    fi
    rm -R outgroups
    if [ ! -d outgroups ]; then
        echo "Downloading BUSCO outgroups"
        wget -c --tries=0 --read-timeout=20 --content-disposition https://uwmadison.box.com/shared/static/4pl3ngptpjjfs1cu4se6g27ei0wptsdt.gz
        tar -zxf busco_outgroups.tar.gz
        echo "-----------------------------------------------"
    else
        echo "BUSCO outgroups found, skipping download"
        echo "-----------------------------------------------"
    fi
    
    if [ ! -f funannotate.repeat.proteins.fa ]; then
        echo "Downloading Repeat Protein DB"
        wget -c --tries=0 --read-timeout=20 --content-disposition https://uwmadison.box.com/shared/static/vcftxq6yuzc3u1nykiahxcqzk3jlvyzx.gz
        tar -zxf funannotate.repeat.proteins.fa.tar.gz
        makeblastdb -in funannotate.repeat.proteins.fa -input_type fasta -dbtype prot -title REPEATS -out REPEATS
        echo "-----------------------------------------------"
    else
        echo "Funannotate repeat protein DB, skipping download"
        echo "-----------------------------------------------"
    fi
            
    if [ ! -f go.obo ]; then
        echo "Downloading Gene Ontology"
        wget -c --tries=0 --read-timeout=20 http://geneontology.org/ontology/go.obo
        echo "-----------------------------------------------"
    else
        echo "Gene Ontology already exists, skipping download. To update delete go.obo."
        echo "-----------------------------------------------"
    fi

    #download MiBIG database for getting best SM hit from curated database.
    if [ ! -f MIBiG_prot_seqs.fa ]; then
        echo "Downloading MIBiG protein fasta files"
        wget -c --tries=0 --read-timeout=20 http://mibig.secondarymetabolites.org/MIBiG_prot_seqs_1.3.fasta
        mv MIBiG_prot_seqs_1.3.fasta MIBiG_prot_seqs.fa
        makeblastdb -in MIBiG_prot_seqs.fa -input_type fasta -dbtype prot -title MIBiG -out MIBiG
        echo "-----------------------------------------------"
    else
        echo "MIBiG database already exists, skipping download"
        echo "-----------------------------------------------"
    fi

    #download InterProScan xml mapping file
    if [ ! -f interpro.xml ]; then
        echo "Downloading InterPro mapping xml file"
        wget -c --tries=0 --read-timeout=20 ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz
        gunzip interpro.xml.gz
        echo "-----------------------------------------------"
    else
        echo "InterPro mapping file already exists, skipping download. To update delete interpro.xml"
        echo "-----------------------------------------------"
    fi
fi

if [ "$dep" = 'pass' ]; then
    #okay now check dependencies and report which are installed and which are not
    echo "Checking External Dependencies...."
    echo "-----------------------------------------------"

    #setup some programs and look for dependencies

    #make sure in funannotate directory
    cd $dir

    check='pass'
    for i in {blastp,hmmsearch,hmmscan,augustus,'gmes_petap.pl',mummer,nucmer,show-coords,exonerate,gmap,blat,RepeatModeler,RepeatMasker,pslCDnaFilter,bedtools,bamtools,'gag.py',tbl2asn,'braker.pl',funannotate,mafft,trimal,raxmlHPC-PTHREADS,tRNAscan-SE,'rmOutToGFF3.pl','proteinortho5.pl'}; do
        var=$(command -v $i)
        if [ "$var" ]; then
            echo "$i installed.........$var"
        else
            echo -e "${RED}ERROR:${NC}  $i is not installed"
            check='fail'
        fi
    done
    echo "-----------------------------------------------"
    echo "Checking Python modules...."
    echo "-----------------------------------------------"
    python -c 'import pkgutil; print("biopython installed"if pkgutil.find_loader("Bio") else "ERROR: biopython not installed")'
    python -c 'import pkgutil; print("psutil installed"if pkgutil.find_loader("psutil") else "ERROR: psutil not installed")'
    python -c 'import pkgutil; print("natsort installed"if pkgutil.find_loader("natsort") else "ERROR: natsort not installed")'
    python -c 'import pkgutil; print("goatools installed"if pkgutil.find_loader("goatools") else "ERROR: goatools not installed")'
    python -c 'import pkgutil; print("sklearn installed"if pkgutil.find_loader("sklearn") else "ERROR: sklearn not installed")'
    python -c 'import pkgutil; print("matplotlib installed"if pkgutil.find_loader("matplotlib") else "ERROR: matplotlib not installed")'
    python -c 'import pkgutil; print("seaborn installed"if pkgutil.find_loader("seaborn") else "ERROR: seaborn not installed")'
    python -c 'import pkgutil; print("numpy installed"if pkgutil.find_loader("numpy") else "ERROR: numpy not installed")'
    python -c 'import pkgutil; print("pandas installed"if pkgutil.find_loader("pandas") else "ERROR: pandas not installed")'

    echo "-----------------------------------------------"
    echo "Checking Perl modules...."
    echo "-----------------------------------------------"

    module_exists() {
    perl -e 'use '$1 2>/dev/null; }

    for i in {Bio::SeqIO,Pod::Usage,File::Basename,threads,threads::shared,Thread::Queue,Getopt::Long,FindBin,File::Spec,File::Path,Data::Dumper,YAML,Carp,Hash::Merge,Logger::Simple,Parallel::ForkManager}; do
        module_exists $i && echo "$i installed" || echo -e "${RED}ERROR:${NC} $i not installed"
    done

    if ! [ "$EVM_HOME" ]; then
        echo -e "${RED}ERROR:${NC}  EVM_HOME variable has not been set
            example: export EVM_HOME=/usr/local/opt/evidencemodeler"
        check='fail'
    fi

    if ! [ "$AUGUSTUS_CONFIG_PATH" ]; then
        echo -e "${RED}ERROR:${NC}  AUGUSTUS_CONFIG_PATH variable has not been set
            example: export AUGUSTUS_CONFIG_PATH=/usr/local/opt/augustus/libexec/config"
        check='fail'
    fi

    if [ "$check" = 'fail' ];then
        echo -e "\nAt least one external dependency needs to be installed before running funannotate\n"
        exit
    else
        echo "-----------------------------------------------"
        echo -e "Script complete, funannotate is ready to roll!\n"
        exit
    fi
fi
