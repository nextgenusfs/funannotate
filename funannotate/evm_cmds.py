"""Command builders for the EVM recombine and convert steps.

Kept out of aux_scripts/funannotate-runEVM.py so they can be tested. The Rust
EVM tools gained --skip_empty, --threads and --concat_out after 1.9.0, and older
builds reject unknown options, so every new option is added only when the
installed binary lists it in --help.
"""

import functools
import os
import shutil
import subprocess


@functools.lru_cache(maxsize=None)
def tool_supports(tool, flag):
    """True if `tool --help` mentions `flag`. False if it cannot be run."""
    try:
        proc = subprocess.run(
            [tool, "--help"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            timeout=60,
        )
    except (OSError, subprocess.SubprocessError):
        return False
    return flag in proc.stdout


def build_recombine_cmd(
    engine, recombine_bin, partitions_file, output_name, threads=0, supports=tool_supports
):
    if engine == "rust":
        cmd = [recombine_bin, "--partitions", os.path.basename(partitions_file), "-O", output_name]
        if threads > 0 and supports(recombine_bin, "--threads"):
            cmd += ["--threads", str(threads)]
        return cmd
    return [
        "perl",
        recombine_bin,
        "--partitions",
        os.path.basename(partitions_file),
        "--output_file_name",
        output_name,
    ]


def build_convert_cmd(
    engine,
    convert_bin,
    partitions_file,
    output_name,
    fasta,
    threads=0,
    concat_out=None,
    supports=tool_supports,
):
    """Return the convert command, or None if the Rust binary is not on PATH.

    `concat_out` is a file name relative to the working directory. It is passed
    as --concat_out only when the binary supports it; check with
    `"--concat_out" in cmd` to see whether it was used.
    """
    if engine == "rust":
        if not shutil.which(convert_bin):
            return None
        cmd = [convert_bin, "--partitions", os.path.basename(partitions_file), "-O", output_name]
        # no empty evm.out.gff3 for each empty evm.out
        if supports(convert_bin, "--skip_empty"):
            cmd.append("--skip_empty")
        if threads > 0 and supports(convert_bin, "--threads"):
            cmd += ["--threads", str(threads)]
        if concat_out and supports(convert_bin, "--concat_out"):
            cmd += ["--concat_out", concat_out]
        return cmd
    return [
        "perl",
        convert_bin,
        "--partitions",
        os.path.basename(partitions_file),
        "--output",
        output_name,
        "--genome",
        os.path.abspath(fasta),
    ]


def collect_gff3(tmpdir, out_path, concat_path=None, gff3_name="evm.out.gff3"):
    """Write all per-contig GFF3 to `out_path`.

    Uses `concat_path` (already in partitions-listing order) when the converter
    wrote it. Otherwise walks `tmpdir` and copies every `gff3_name`, as before.
    Returns "concat" or "walk".
    """
    if concat_path and os.path.isfile(concat_path):
        shutil.move(concat_path, out_path)
        return "concat"
    with open(out_path, "w") as out:
        for root, _dirs, files in os.walk(tmpdir):
            for name in files:
                if name == gff3_name:
                    with open(os.path.join(root, name), "r") as handle:
                        shutil.copyfileobj(handle, out)
    return "walk"
