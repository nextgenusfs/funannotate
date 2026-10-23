#!/usr/bin/env bash
# Run the tools of the funannotate-base image for real, not just --version: Trinity on its sample data, a MySQL
# server (the conda mysqld, and the apt mariadbd) with DBD::mysql talking to it (PASA's --pasa_db mysql), and the
# Perl/Python/R libraries. Run inside the image as the runtime user:
#   docker run --rm -v "$PWD/install_scripts:/smoke:ro" <image> bash /smoke/image_smoke_tests.sh
# Used by .github/workflows/image-prune-test.yml after prune_runtime_env.sh. Every check prints "ok:" or "FAIL:";
# exit status 1 if any failed (all checks still run).

set -uo pipefail

: "${PIXI_ENV_PATH:=/pixi/.pixi/envs/base}"
work="$(mktemp -d)"
fail=0
ok() { echo "ok:   $1"; }
bad() { echo "FAIL: $1"; fail=1; }
check() { # check "description" command...
    local desc="$1"
    shift
    if "$@" > "${work}/last.log" 2>&1; then ok "${desc}"; else bad "${desc}"; sed 's/^/        /' "${work}/last.log" | tail -8; fi
}
has() { grep -q "$1" "${work}/last.log"; }

echo "== tool versions"
check "java" java -version
check "salmon 2.x" bash -c 'salmon --no-version-check --version | grep -qE "^salmon 2\."'
check "blastn" blastn -version
check "rmblastn" rmblastn -version
check "hisat2" hisat2 --version
check "bowtie2" bowtie2 --version
check "samtools" samtools --version
check "augustus" bash -c 'augustus --version 2>&1 | grep -q AUGUSTUS'
check "RepeatMasker" RepeatMasker -v
check "mmseqs" mmseqs version
check "diamond" diamond version
check "Trinity --version" Trinity --version
check "evidence_modeler present" test -x "${PIXI_ENV_PATH}/bin/evidence_modeler"
check "Perl EVM scripts present" test -f "${EVM_HOME:-${PIXI_ENV_PATH}/opt/evm}/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl"
check "PASA launcher present" test -f "${PASAHOME:-${PIXI_ENV_PATH}/opt/pasa/src}/Launch_PASA_pipeline.pl"
check "PASA launcher loads its Perl libraries" bash -c 'perl "${PASAHOME:-'"${PIXI_ENV_PATH}"'/opt/pasa/src}/Launch_PASA_pipeline.pl" -h 2>&1 | grep -qi "pasa"'

echo "== more tools funannotate calls"
present() { check "$1 on PATH" command -v "$1"; }
for t in gmap fasta36 blat pblat snap glimmerhmm tbl2asn tRNAscan-SE CodingQuarry RepeatModeler kallisto stringtie fastp pigz tantan bedtools mafft trimal iqtree proteinortho pslCDnaFilter; do
    present "$t"
done
check "minimap2" minimap2 --version
check "miniprot" miniprot --version
check "exonerate" bash -c 'exonerate --version 2>&1 | grep -q "exonerate version"'
check "hmmsearch" hmmsearch -h
check "eggnog-mapper" emapper.py --version
echo "   (TransDecoder files: $(ls "${PIXI_ENV_PATH}"/bin | grep -i transdecoder | tr '\n' ' ') / $(ls -d "${PIXI_ENV_PATH}"/opt/transdecoder* 2> /dev/null | tr '\n' ' '))"
echo "   (TransDecoder.LongOrfs: $(ls -l "${PIXI_ENV_PATH}/bin/TransDecoder.LongOrfs" 2>&1 | cut -c1-200); running as: $(id -un); first line: $(head -c 120 "${PIXI_ENV_PATH}/bin/TransDecoder.LongOrfs" 2>&1 | head -1))"
check "TransDecoder.LongOrfs on a made-up transcript" bash -c '
    d=$(mktemp -d) && cd "$d" &&
    printf ">t1\nGGGATGGCTAGCATCGATCGGCTAGCTAGCATCGATCGATCGGCTAGCATCGATCGGATCGATCGGCTAGCTAGCATCGATCGATCGGCTAGCATCGATCGGATCGATCGGCTAGCTAGCATCGATCGATTAAGGG\n" > t.fa &&
    TransDecoder.LongOrfs -t t.fa -O out -m 30 > log 2>&1 && test -f out/longest_orfs.pep || { cat log; exit 1; }'
check "R: seqlogo and ggplot2 (what TransDecoder needs)" Rscript -e 'library(seqLogo); library(ggplot2)'

echo "== libraries"
check "perl: DBI DBD::mysql DBD::SQLite DB_File" perl -MDBI -MDBD::mysql -MDBD::SQLite -MDB_File -e 1
check "python: numpy pandas scipy Bio sklearn matplotlib" python -c 'import numpy, pandas, scipy, Bio, sklearn, matplotlib'
check "R" Rscript -e 'invisible(library(stats))'

echo "== MySQL (PASA --pasa_db mysql): server, client and DBD::mysql"
mysql_roundtrip() { # mysql_roundtrip name server-binary init-cmd...
    local name="$1" server="$2" d="${work}/${1}"
    shift 2
    mkdir -p "${d}"
    "$@" > "${d}/init.log" 2>&1 || { bad "${name}: initialize the data directory"; tail -5 "${d}/init.log"; return; }
    "${server}" --no-defaults --datadir="${d}/data" --socket="${d}/s.sock" --pid-file="${d}/pid" --skip-networking \
        > "${d}/server.log" 2>&1 &
    local pid=$! i
    for i in $(seq 1 60); do [ -S "${d}/s.sock" ] && break; sleep 1; done
    if [ ! -S "${d}/s.sock" ]; then bad "${name}: server did not start"; tail -8 "${d}/server.log"; kill "${pid}" 2> /dev/null; return; fi
    ok "${name}: server started"
    check "${name}: SQL through the mysql client" bash -c "mysql --no-defaults -S '${d}/s.sock' -u root -e 'create database t; use t; create table x (i int); insert into x values (1); select * from x'"
    check "${name}: DBD::mysql connects" perl -MDBI -e '$d = DBI->connect("DBI:mysql:mysql_socket=$ARGV[0];user=root", "", "", {RaiseError => 1}); print $d->selectrow_array("select 1"), "\n"' "${d}/s.sock"
    kill "${pid}" 2> /dev/null
    wait "${pid}" 2> /dev/null
}
mysql_roundtrip conda-mysqld "$(command -v mysqld)" "$(command -v mysqld)" --no-defaults --initialize-insecure --datadir="${work}/conda-mysqld/data"
if command -v mariadbd > /dev/null; then
    mysql_roundtrip apt-mariadbd "$(command -v mariadbd)" "$(command -v mariadb-install-db)" --no-defaults --datadir="${work}/apt-mariadbd/data" --auth-root-authentication-method=normal
else
    bad "apt mariadbd is not installed"
fi

echo "== Trinity assembly on its sample data (Inchworm, Chrysalis, Butterfly)"
sample="${TRINITY_HOME:-${PIXI_ENV_PATH}/opt/trinityrnaseq}/sample_data/test_Trinity_Assembly"
left="$(ls "${sample}"/reads.left.fq* 2> /dev/null | head -1)"
right="$(ls "${sample}"/reads.right.fq* 2> /dev/null | head -1)"
if [ -z "${left}" ] || [ -z "${right}" ]; then
    bad "Trinity sample reads not found in ${sample}"
else
    if timeout 1500 Trinity --seqType fq --left "${left}" --right "${right}" --SS_lib_type RF --max_memory 2G --CPU 2 \
        --output "${work}/trinity_test_out" > "${work}/trinity.log" 2>&1 && [ -s "${work}/trinity_test_out.Trinity.fasta" -o -s "${work}/trinity_test_out/Trinity.fasta" ]; then
        ok "Trinity assembled $(grep -c '^>' "${work}"/trinity_test_out*Trinity.fasta "${work}/trinity_test_out/Trinity.fasta" 2> /dev/null | head -1)"
    else
        bad "Trinity assembly on the sample data"
        tail -25 "${work}/trinity.log"
    fi
fi

rm -rf "${work}"
echo
if [ "${fail}" = "0" ]; then echo "all image smoke tests passed"; else echo "some image smoke tests FAILED"; fi
exit "${fail}"
