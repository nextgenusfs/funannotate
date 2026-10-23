#!/bin/bash
# Sourced by pixi on every activation of the osx-arm64 environment. Keep it fast and
# idempotent, and do not `exit` (the script is sourced).
#
# The conda perl-db_file package for osx-arm64 crashes on "use DB_File" (it does not link
# libdb), and its copy in site_perl hides perl's working core copy. Put perl's core library
# first so that Trinity's insilico_read_normalization.pl runs. PASA's hook libraries are added
# as on Linux (see docs/conda.rst).
if [ -n "${CONDA_PREFIX:-}" ]; then
    _fa_core=$(ls -d "${CONDA_PREFIX}"/lib/perl5/*/core_perl 2> /dev/null | head -1)
    for _fa_dir in "${PASAHOME:-${CONDA_PREFIX}/opt/pasa/src}/PerlLib" "${PASAHOME:-${CONDA_PREFIX}/opt/pasa/src}/SAMPLE_HOOKS" "${_fa_core}"; do
        if [ -n "${_fa_dir}" ]; then
            case ":${PERL5LIB:-}:" in
                *":${_fa_dir}:"*) ;;
                *) export PERL5LIB="${_fa_dir}${PERL5LIB:+:${PERL5LIB}}" ;;
            esac
        fi
    done
    unset _fa_core _fa_dir
fi
