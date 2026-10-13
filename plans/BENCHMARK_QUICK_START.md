# Trinity Rust Optimization - Quick Benchmark Guide

## 30-Second Setup

```bash
# 1. Enter pixi environment
pixi shell

# 2. Run benchmark (will take 1-3 hours depending on hardware)
./scripts/benchmark_trinity_rust.sh \
    ~/projects/funannotate/trinity_example/Cordyceps_militaris \
    ~/bench_results

# 3. View results
cat ~/bench_results/results/BENCHMARK_REPORT.md
```

## Manual Testing (Quick)

### Check which version is active
```bash
./scripts/toggle_trinity_rust.sh status
```

### Switch to Perl version (uses slower scripts)
```bash
./scripts/toggle_trinity_rust.sh disable
# Then run: funannotate train ...
```

### Switch back to Rust (optimized)
```bash
./scripts/toggle_trinity_rust.sh enable
# Then run: funannotate train ...
```

## Understanding the Output

### Report Example

```
# Elapsed Time
Rust:  45 minutes
Perl:  120 minutes
Speedup: 2.67x faster with Rust
```

### Expected Performance

Based on component analysis:

| Tool | Expected Speedup |
|------|------------------|
| sam_to_read_coords | 5-20x |
| fragment_coverage_writer | 10-50x |
| define_coverage_partitions | 3-15x |
| extract_reads_per_partition | 10-40x |
| **Overall workflow** | **5-15x** |

Note: Overall is less than individual tools because:
- Other pipeline stages (Inchworm, Chrysalis, Butterfly) still use Perl/C++
- Not all time is spent in optimized utilities
- I/O overhead from temporary files

## What Gets Tested

The benchmark runs `funannotate train` (RNA-seq training pipeline) which includes:

1. **HISAT2 alignment** → BAM file
2. **Trinity genome-guided assembly** (this is where Rust tools help):
   - `sam_to_read_coords` - Parse BAM file coordinates
   - `fragment_coverage_writer` - Generate coverage from fragments
   - `define_coverage_partitions` - Partition genome by coverage
   - `extract_reads_per_partition` - Extract reads for each partition
   - Inchworm k-mer assembly (per partition)
   - Chrysalis graph clustering
   - Butterfly contig reconstruction
3. **PASA transcript alignment** to refined genome

## File Structure

```
~/bench_results/
├── logs/                      # Detailed execution logs
│   ├── rust_enabled.log
│   └── perl_only.log
├── results/                   # Metrics and report
│   ├── rust_enabled.metrics
│   ├── perl_only.metrics
│   └── BENCHMARK_REPORT.md   # Main result file
├── work_rust/                 # Output with Rust tools
│   └── output/transcripts.fasta
└── work_perl/                 # Output with Perl tools
    └── output/transcripts.fasta
```

## Troubleshooting

### "Command not found" errors

```bash
# Ensure pixi environment is active
pixi shell
which funannotate  # should return a path
which sam_to_read_coords  # Rust tool should exist
```

### Benchmark takes too long

The benchmark uses 4 CPUs and 16GB memory to be conservative. For faster tests:

1. Edit `scripts/benchmark_trinity_rust.sh` (around line 200)
2. Change `--cpus 4 --memory 16G` to `--cpus 2 --memory 8G`
3. Results will be less accurate but run faster

### One run failed

This is okay! Check the log:
```bash
tail -50 ~/bench_results/logs/rust_enabled.log  # For Rust run
tail -50 ~/bench_results/logs/perl_only.log     # For Perl run
```

Common issues:
- Disk space: Trinity needs ~10GB free
- Memory: 16GB may not be enough for large datasets
- Dependencies: `FUNANNOTATE_DB` environment variable

### Transcripts differ between runs

**This is normal!** Trinity has random assembly steps, so:
- Rust vs Perl may produce different results
- Different CPUs may produce different results
- Both are valid; focus on timing, not transcript count

## Advanced Workflows

### Run just one version manually

```bash
# With Rust (fast):
time funannotate train -i genome.fasta -o output_rust \
    --left_norm R1.fastq.gz --right_norm R2.fastq.gz \
    -s "Species" --cpus 4 --memory 16G

# Without Rust (slow):
./scripts/toggle_trinity_rust.sh disable
time funannotate train -i genome.fasta -o output_perl \
    --left_norm R1.fastq.gz --right_norm R2.fastq.gz \
    -s "Species" --cpus 4 --memory 16G
./scripts/toggle_trinity_rust.sh enable
```

### Benchmark different datasets

```bash
# Your dataset
./scripts/benchmark_trinity_rust.sh \
    /path/to/your/data \
    ~/results_mydata

# Multiple datasets
for dataset in dataset1 dataset2 dataset3; do
    ./scripts/benchmark_trinity_rust.sh \
        /data/$dataset \
        ~/results/$dataset
done
```

## Key Scripts

| Script | Purpose |
|--------|---------|
| `scripts/benchmark_trinity_rust.sh` | Full automated benchmark |
| `scripts/toggle_trinity_rust.sh` | Quick manual testing |
| `scripts/pixi_install_trinity_rust.sh` | Build Rust tools (auto on pixi activation) |

## Related Documentation

- **Full guide**: `BENCHMARK_TRINITY_RUST.md`
- **Integration details**: `TRINITY_RUST_INTEGRATION.md`
- **Tools replaced**: `TRINITY_TOOLS_REPLACEMENT.md`

## Questions?

Check the full documentation:
```bash
less BENCHMARK_TRINITY_RUST.md        # Complete guide
less TRINITY_RUST_INTEGRATION.md      # Technical details
less TRINITY_TOOLS_REPLACEMENT.md     # What was changed
```

## Key Takeaway

Rust Trinity tools are **drop-in replacements** - no code changes needed, automatic speedup via pixi activation.

Benchmark shows you exactly how much faster! 🚀
