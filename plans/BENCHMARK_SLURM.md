# Running Trinity Rust Benchmark on HPCC (SLURM)

## Quick Start

```bash
# Navigate to funannotate directory
cd /bigdata/stajichlab/jstajich/projects/funannotate/funannotate-live

# Create log directory (one-time)
mkdir -p slurm_logs

# Submit benchmark job
sbatch benchmark_trinity_rust.sbatch

# Monitor progress
squeue -u $USER
tail -f slurm_logs/trinity_rust_bench_<job_id>.out

# View results (when done)
cat ~/bench_results/results/BENCHMARK_REPORT.md
```

## SLURM Configuration

The sbatch script is configured with:

```bash
--partition=epyc           # Longer jobs (up to 7 days)
--cpus-per-task=4          # 4 CPU cores
--mem=32GB                 # 32 GB memory (allows for overhead)
--time=0-08:00:00          # 8 hour time limit
--account=stajichlab       # Lab account
```

### Why These Settings?

- **epyc partition**: Needed for potentially 2-4 hour runtime (exceeds 2 hour `short` partition limit)
- **4 CPUs**: Matches benchmark configuration; more CPUs won't help much (Trinity is single-threaded)
- **32GB memory**: 16GB for Trinity + buffer for overhead
- **8 hours**: Safe margin for occasional longer runs

## Adjusting Resources

### For Faster Tests (less accurate)

Edit `benchmark_trinity_rust.sbatch`:

```bash
# Change:
#SBATCH --time=0-08:00:00    # to:
#SBATCH --time=0-02:00:00    # 2 hours (can use 'short' partition)

#SBATCH --partition=epyc     # to:
#SBATCH --partition=short    # Shorter queue, faster start

#SBATCH --mem=32GB           # to:
#SBATCH --mem=16GB           # Less memory overhead
```

### For Larger Datasets

```bash
#SBATCH --cpus-per-task=8    # Use more CPUs
#SBATCH --mem=64GB           # More memory
#SBATCH --time=0-12:00:00    # Longer runtime
```

## Available SLURM Partitions on UCR HPCC

| Partition | Max Time | Use Case |
|-----------|----------|----------|
| `short` | 2 hours | Quick tests |
| `short_gpu` | 2 hours | GPU compute |
| `stajichlab` | 14 days | Lab nodes (priority) |
| `epyc` | 7 days | Long-running compute |
| `highmem` | 7 days | Memory-intensive jobs |

**This benchmark fits best on**: `stajichlab` or `epyc`

## Monitoring Your Job

### Check job status
```bash
# Your jobs
squeue -u $USER

# Specific job
squeue -j <job_id>

# With more details
squeue -u $USER -o "%.18i %.9P %.20j %.8u %.8T %.10M %.9l %.6D %R"
```

### Watch output in real-time
```bash
tail -f slurm_logs/trinity_rust_bench_<job_id>.out

# Or with 'follow' mode (auto-refresh)
watch -n 2 'tail -20 slurm_logs/trinity_rust_bench_*.out'
```

### Check resource usage (while running)
```bash
# If job is running
sstat -j <job_id> --format=AveCPU,AveVMSize,MaxVMSize,MaxRSS
```

## Job Output Files

After submission, you'll have:

```
slurm_logs/
├── trinity_rust_bench_<job_id>.out    # Standard output
└── trinity_rust_bench_<job_id>.err    # Standard error

~/bench_results/
├── logs/
│   ├── rust_enabled.log      # Rust version execution log
│   └── perl_only.log         # Perl version execution log
├── results/
│   ├── rust_enabled.metrics  # Timing/memory metrics (Rust)
│   ├── perl_only.metrics     # Timing/memory metrics (Perl)
│   └── BENCHMARK_REPORT.md   # Main report with comparison
├── work_rust/                # Trinity output (with Rust)
└── work_perl/                # Trinity output (with Perl)
```

## Viewing Results

### After job completes

```bash
# View the main report
cat ~/bench_results/results/BENCHMARK_REPORT.md

# View detailed metrics
cat ~/bench_results/results/rust_enabled.metrics
cat ~/bench_results/results/perl_only.metrics

# View execution logs
tail -100 ~/bench_results/logs/rust_enabled.log
tail -100 ~/bench_results/logs/perl_only.log

# Check generated transcripts
ls -lh ~/bench_results/work_rust/output/transcripts.fasta
ls -lh ~/bench_results/work_perl/output/transcripts.fasta
```

## Common Issues & Solutions

### Job fails with "pixi not found"

The system may not have pixi in the default PATH. Check if it's available:

```bash
# Try loading a module (if available)
module avail pixi
module load pixi  # if available

# Or check where pixi is installed
which pixi
find ~ -name "pixi" -type f 2>/dev/null | head -1
```

If pixi is installed, update the sbatch to load it:
```bash
#SBATCH --partition=epyc
# Add this line:
module load pixi  # if available, or adjust to correct module name
```

### Job times out (hits time limit)

Increase time in sbatch:
```bash
#SBATCH --time=0-12:00:00    # Increase to 12 hours
```

Or adjust benchmark to use fewer resources:
```bash
# Edit scripts/benchmark_trinity_rust.sh
# Around line 200, change:
# --cpus 4 --memory 16G
# to:
# --cpus 2 --memory 8G
```

### Insufficient memory error

Increase memory allocation:
```bash
#SBATCH --mem=48GB    # Increase from 32GB
```

### Dataset not found error

Ensure dataset path is correct:
```bash
# Verify dataset exists
ls ~/projects/funannotate/trinity_example/Cordyceps_militaris/

# If it doesn't, update the sbatch script:
# DATASET_DIR="path/to/your/dataset"
```

## Email Notifications

The sbatch script sends emails on job completion or failure. Update your email:

```bash
#SBATCH --mail-user=your.email@ucr.edu
```

Or disable notifications:
```bash
# Remove or comment out these lines:
# #SBATCH --mail-type=END,FAIL
# #SBATCH --mail-user=...
```

## Running Multiple Benchmarks

You can submit multiple jobs to compare different scenarios:

```bash
# Benchmark 1: Default (4 CPUs, 16GB)
sbatch benchmark_trinity_rust.sbatch

# Benchmark 2: More resources (8 CPUs, 32GB)
# Create a copy with different settings
cp benchmark_trinity_rust.sbatch benchmark_trinity_rust_8cpu.sbatch
# Edit: --cpus-per-task=8, --mem=64GB
sbatch benchmark_trinity_rust_8cpu.sbatch

# Check all your jobs
squeue -u $USER
```

## Performance Expectations

Based on test dataset (33 MB genome, 250+263 MB RNA-seq):

| Phase | Estimated Time |
|-------|-----------------|
| Pixi setup | 2-5 min |
| HISAT2 alignment | 10-15 min |
| Trinity genome-guided | 60-180 min |
| PASA refinement | 20-40 min |
| **Total** | **90-240 min** |

With 2x speedup from Rust:
- **Rust version**: ~60-120 minutes
- **Perl version**: ~120-240 minutes

## Advanced: Running on Different Partition

To use lab's dedicated partition (if available):

```bash
#SBATCH --partition=stajichlab    # Instead of epyc
#SBATCH --account=stajichlab
```

This gives priority access and may start faster.

## Related Documentation

- `BENCHMARK_QUICK_START.md` - Quick local testing guide
- `BENCHMARK_TRINITY_RUST.md` - Full benchmarking documentation
- [UCR HPCC Wiki](http://cluster.hpcc.ucr.edu/) - SLURM documentation

## Questions?

```bash
# Check SLURM help
man sbatch
man squeue

# Check job details
sinfo -n <node_name>
slurm-stat  # If available on your system
```
