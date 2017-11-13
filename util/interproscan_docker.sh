#!/bin/bash

#Run Docker interproscan
if [ -z "$2" ]; then
    echo "Usage: interproscan_docker.sh -i=input_proteins.fa -c=CPUS"
    exit
fi


#parse input commands to script
for i in "$@"
do
case $i in
    -i=*|--input=*)
    input="${i#*=}"
    shift # past argument=value
    ;;
    -c=*|--cpus=*)
    cpus="${i#*=}"
    shift # past argument=value
    ;;
    *)
            # unknown option
    ;;
esac
done
echo "InputProteins    = ${input}"
echo "Num CPUS  = ${cpus}"

#download the interproscan properties file
curl https://raw.githubusercontent.com/ebi-pf-team/interproscan/5.22-61.0/core/jms-implementation/support-mini-x86-32/interproscan.properties > interproscan.properties.tmp
export var=${cpus}; perl -plne 's/^(number.of.embedded.workers).*/$1=1/;s/^(maxnumber.of.embedded.workers).*/$1=$ENV{var}/' interproscan.properties.tmp > interproscan.properties
rm interproscan.properties.tmp

#now setup the docker run
docker run -u $UID:$GROUPS --rm \
    -v `pwd`:/dir \
    -v `pwd`:/in \
    -v `pwd`/interproscan.properties:/interproscan-5.22-61.0/interproscan.properties \
    blaxterlab/interproscan:latest \
    interproscan.sh -i /in/$input -d /dir -dp -f XML -goterms -pa