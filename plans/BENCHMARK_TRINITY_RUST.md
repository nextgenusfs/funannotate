# Benchmarking Trinity Rust Optimization

## Overview

This document describes how to benchmark funannotate's Trinity workflow performance with and without Rust-optimized utilities.

## Quick Start

### Run Full Benchmark

```bash
# Make sure pixi environment is activated
pixi shell
# or
eval "$(pixi shell-hook -s bash)"

# Run benchmark with test dataset
./scripts/benchmark_trinity_rust.sh \
    ~/projects/funannotate/trinity_example/Cordyceps_militaris \
    ~/benchmark_results
```

### What This Does

The benchmark script:

1. **Validates** input dataset (genome + paired-end RNA-seq reads)
2. **Checks** dependencies (funannotate, Trinity, Rust tools, etc.)
3. **Runs** funannotate train with **Rust Trinity tools enabled**
   - Measures: Elapsed time, peak memory, CPU usage, transcript count
4. **Runs** same workflow with **Rust tools disabled** (uses Perl versions)
   - Uses same dataset, same parameters for fair comparison
5. **Generates** detailed comparison report with metrics

### Output Structure

```
~/benchmark_results/
├── logs/
│   ├── rust_enabled.log      # Detailed log from Rust run
│   └── perl_only.log         # Detailed log from Perl run
├── results/
│   ├── rust_enabled.metrics  # /usr/bin/time -v output (Rust)
│   ├── perl_only.metrics     # /usr/bin/time -v output (Perl)
│   └── BENCHMARK_REPORT.md   # Formatted comparison report
├── work_rust/                # Trinity output with Rust tools
│   └── output/
│       └── transcripts.fasta
└── work_perl/                # Trinity output with Perl tools
    └── output/
        └── transcripts.fasta
```

## Benchmark Details

### Test Dataset

The benchmark uses RNA-seq data from *Cordyceps militaris*:

```
Genome:       Cordyceps_militaris_ATCC_34164.fasta (33 MB)
Left reads:   Cordyceps_militaris_norm_R1.fastq.gz (250 MB)
Right reads:  Cordyceps_militaris_norm_R2.fastq.gz (263 MB)
Total:        ~546 MB (good size for meaningful benchmark)
```

### Benchmark Parameters

- **funannotate train** workflow
- **4 CPUs** (limited to avoid system saturation)
- **16GB memory** limit
- **--jaccard_clip** enabled for dense genomes
- **Two independent runs** (Rust vs. Perl)

### What Gets Measured

**From `/usr/bin/time -v`:**
- Elapsed time (wall clock)
- User CPU time
- System CPU time
- Maximum resident set size (memory)
- Page faults (major/minor)
- Context switches
- I/O statistics

**From Trinity output:**
- Number of transcripts assembled
- Total bases assembled
- Transcript statistics

## Understanding Results

### Key Metrics

| Metric | What It Means |
|--------|--------------|
| **Elapsed** | Total wall-clock time (what you experience) |
| **User CPU** | Time spent running code |
| **Sys CPU** | Time spent in system calls/I/O |
| **Max RSS** | Peak memory usage in KB |
| **Minor faults** | Page swaps from RAM (normal) |
| **Major faults** | Page swaps from disk (slow!) |

### Performance Indicators

**Speedup = Perl Time / Rust Time**
- 1.0x = No difference
- 2.0x = Rust is twice as fast
- 5.0x+ = Significant improvement

### Expected Improvements

Based on Trinity internals, expect **5-50x speedup** for:
- `sam_to_read_coords`: SAM parsing (5-20x faster)
- `fragment_coverage_writer`: Coverage generation (10-50x faster)
- `define_coverage_partitions`: Partition definition (3-15x faster)
- `extract_reads_per_partition`: Read extraction (10-40x faster)

## Advanced Usage

### Custom Dataset

```bash
./scripts/benchmark_trinity_rust.sh \
    /path/to/your/genome_dir \
    ~/my_benchmark_results
```

Requirements:
- One `.fasta` or `.fa` file (genome)
- At least one `*R1*.fastq.gz` file (left reads)
- At least one `*R2*.fastq.gz` file (right reads)

### Analyzing Results

```bash
# View detailed metrics
cat ~/benchmark_results/results/rust_enabled.metrics
cat ~/benchmark_results/results/perl_only.metrics

# View logs for troubleshooting
tail -f ~/benchmark_results/logs/rust_enabled.log
tail -f ~/benchmark_results/logs/perl_only.log

# View generated report
cat ~/benchmark_results/results/BENCHMARK_REPORT.md
```

### Manual Comparison

If you want to run individual tests:

```bash
# Run with Rust tools (normal PATH)
funannotate train -i genome.fasta -o output_rust \
    --left_norm R1.fastq.gz --right_norm R2.fastq.gz \
    -s "Species" --cpus 4 --memory 16G

# Run with Rust tools hidden (use Perl)
# Temporarily rename Rust binaries
CONDA_PREFIX=$(python -c "import sys; print(sys.prefix)")
for tool in sam_to_read_coords extract_reads_per_partition \
            fragment_coverage_writer define_coverage_partitions; do
    [ -f "$CONDA_PREFIX/bin/$tool" ] && \
        mv "$CONDA_PREFIX/bin/$tool" "$CONDA_PREFIX/bin/${tool}.disabled"
done

funannotate train -i genome.fasta -o output_perl \
    --left_norm R1.fastq.gz --right_norm R2.fastq.gz \
    -s "Species" --cpus 4 --memory 16G

# Restore Rust tools
for tool in sam_to_read_coords extract_reads_per_partition \
            fragment_coverage_writer define_coverage_partitions; do
    [ -f "$CONDA_PREFIX/bin/${tool}.disabled" ] && \
        mv "$CONDA_PREFIX/bin/${tool}.disabled" "$CONDA_PREFIX/bin/$tool"
done
```

## Troubleshooting

### "Some dependencies are missing"

Make sure pixi environment is activated:
```bash
pixi shell
# or
eval "$(pixi shell-hook -s bash)"
```

### "Rust tools not found"

Ensure Trinity Rust was built during pixi activation:
```bash
pixi install  # rebuilds environment and runs build scripts
which sam_to_read_coords  # should show path
```

### "funannotate train failed"

Check the log file:
```bash
tail -100 ~/benchmark_results/logs/rust_enabled.log
```

Common issues:
- `FUNANNOTATE_DB` not set
- Insufficient disk space in output directory
- Memory limits too low for large datasets

### Transcript Count Differs

This is **normal**. Trinity has non-deterministic assembly steps, so:
- Different CPU counts may produce slightly different results
- Rust vs Perl may assemble reads differently (both valid)
- Focus on **elapsed time** for performance comparison, not transcript count

## Performance Tuning

### Adjust for Faster Results

```bash
# Use fewer CPUs/memory for quicker tests (less accurate):
./scripts/benchmark_trinity_rust.sh ... --cpus 2 --memory 8G

# Or modify the script directly (line ~200):
# Change: --cpus 4 --memory 16G
# To:     --cpus 2 --memory 8G
```

### Adjust for Better Accuracy

```bash
# Use more CPUs/memory for thorough benchmarks:
# Modify script to use --cpus 8 --memory 32G

# Or run custom test:
funannotate train -i genome.fasta -o output \
    --left_norm R1.fastq.gz --right_norm R2.fastq.gz \
    -s "Species" --cpus 8 --memory 32G --jaccard_clip
```

## CI/CD Integration

### Run As Part of Tests

```bash
# In a CI pipeline
set -e
pixi shell
./scripts/benchmark_trinity_rust.sh \
    ~/projects/funannotate/trinity_example/Cordyceps_militaris \
    /tmp/benchmark_results

# Check if Rust was faster than Perl
grep "Speedup" /tmp/benchmark_results/results/BENCHMARK_REPORT.md | grep -E "[2-9]\..*x|^[5-9].*x"
```

## Related Documentation

- `TRINITY_RUST_INTEGRATION.md` - Technical overview of Rust tools
- `TRINITY_TOOLS_REPLACEMENT.md` - Detailed tool replacement analysis
- `../trinityrnaseq/BENCHMARK_SUMMARY.md` - Upstream Trinity benchmarks

## References

- Trinity source: https://github.com/hyphaltip/trinityrnaseq
- Benchmark reference data: `~/projects/funannotate/trinity_example/`
