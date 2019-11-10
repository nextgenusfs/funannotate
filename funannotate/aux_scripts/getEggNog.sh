#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage: getEggNog.sh fuNOG directory"
    exit
fi


EGGNOG=$1
wget -c --tries=0 --read-timeout=20 http://eggnogdb.embl.de/download/eggnog_4.5/data/$EGGNOG/$EGGNOG.hmm.tar.gz
wget -c --tries=0 --read-timeout=20 http://eggnogdb.embl.de/download/eggnog_4.5/data/$EGGNOG/$EGGNOG.annotations.tsv.gz
gunzip $EGGNOG.annotations.tsv.gz
tar -zxf $EGGNOG.hmm.tar.gz
find $EGGNOG\_hmm/ -maxdepth 1 -type f -name '*.hmm' -exec cat '{}' \; > $EGGNOG\_4.5.hmm
hmmpress $EGGNOG\_4.5.hmm
rm $EGGNOG.hmm.tar.gz
rm -R $EGGNOG\_hmm/
for i in $EGGNOG\*; do
    mv $i $2/
done
echo "Done, $EGGNOG DB is now ready to use"