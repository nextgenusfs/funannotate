# Troubleshooting Rust Tool Installations

## Common Error: "destination path already exists and is not an empty directory"

### Symptom

When running `pixi install` or submitting an sbatch job, you see errors like:

```
fatal: destination path '/bigdata/stajichlab/jstajich/projects/funannotate/funannotate-live/.pixi/envs/default/opt/evm-rust/src' 
already exists and is not an empty directory.
```

### Cause

A previous installation attempt failed or was interrupted, leaving incomplete/corrupt files in the installation directory. When `pixi install` tries to rebuild, the installation scripts attempt to clone fresh copies but fail because the directories already exist.

### Solution

#### Option 1: Automatic Cleanup (Recommended)

Use the cleanup utility script:

```bash
# Check what would be cleaned (dry-run)
./scripts/cleanup_rust_builds.sh

# Actually remove incomplete installations
./scripts/cleanup_rust_builds.sh --force

# Rebuild from scratch
pixi install
```

#### Option 2: Manual Cleanup

If you're in a pixi environment:

```bash
pixi shell

# Remove incomplete EVM installation
rm -rf ${CONDA_PREFIX}/opt/evm-rust/src

# Remove incomplete PASA installation  
rm -rf ${CONDA_PREFIX}/opt/pasa-rust-3.0/src

# Clean up Trinity build artifacts (optional)
cd /path/to/trinityrnaseq/rust_bio_utils
cargo clean

# Rebuild
pixi install
```

If you're NOT in a pixi environment:

```bash
# Navigate to funannotate directory
cd /bigdata/stajichlab/jstajich/projects/funannotate/funannotate-live

# Activate pixi shell
pixi shell

# Clean up and rebuild (as above)
```

#### Option 3: Force Rebuild in sbatch

Update `benchmark_trinity_rust.sbatch` to automatically clean before rebuilding:

```bash
# Uncomment these lines in the PIXI ENVIRONMENT ACTIVATION section:
rm -rf .pixi/envs/default/opt/evm-rust/src
rm -rf .pixi/envs/default/opt/pasa-rust-3.0/src

# Then run:
sbatch benchmark_trinity_rust.sbatch
```

## Verification After Cleanup

After cleaning up, verify the installation succeeded:

```bash
pixi shell

# Check Rust Trinity tools
which sam_to_read_coords
which extract_reads_per_partition
which fragment_coverage_writer
which define_coverage_partitions

# Check Rust EVM
which evidence_modeler

# Check Rust PASA  
which pasa_rust
```

All should return paths in `${CONDA_PREFIX}/bin/`.

## If Problems Persist

### Check environment variables

```bash
pixi shell
echo $CONDA_PREFIX
echo $FUNANNOTATE_EVM_ENGINE
echo $PASAHOME
echo $TRINITY_HOME
```

### Verify pixi activation scripts ran

```bash
# Check if build scripts exist
ls -la scripts/pixi_install_*.sh

# Manually run them to see detailed output
bash scripts/pixi_install_rust_evm.sh
bash scripts/pixi_install_rust_pasa.sh  
bash scripts/pixi_install_trinity_rust.sh
```

### Check for disk space issues

```bash
# Make sure you have enough space
df -h /bigdata/
du -sh .pixi/

# If low on space, clean up:
pixi project clean
```

### Check git/network connectivity

If cloning from GitHub fails:

```bash
# Test git connectivity
git ls-remote https://github.com/hyphaltip/EVidenceModeler_rust.git

# Or use local checkouts instead:
# ../EVidenceModeler_rust
# ../PASA_rust  
# ../trinityrnaseq (already being used for Trinity)
```

## Related Issues

### "CONDA_PREFIX is not set"

```bash
# Solution: You're not in a pixi environment
pixi shell
# Then run cleanup or installation commands
```

### "pixi command not found"

```bash
# Solution: pixi not in PATH
# Either:
# 1. Add pixi to PATH (if installed)
# 2. Load pixi module: module load pixi
# 3. Use full path: /path/to/pixi shell
```

### "Permission denied" errors

```bash
# Make sure scripts are executable
chmod +x scripts/*.sh
chmod +x benchmark_trinity_rust.sbatch

# And you have write permissions to:
# - .pixi/ directory
# - ~/bench_results/
# - /tmp (or temporary directories)
```

## Clean Rebuild from Scratch

If you want to completely rebuild the pixi environment:

```bash
# 1. Clean up everything
pixi project clean
rm -rf .pixi

# 2. Rebuild fresh
pixi install

# 3. Verify
pixi shell
which sam_to_read_coords
```

## Debugging: Verbose Output

To see detailed installation output:

```bash
# Run pixi with verbose flag
pixi install -v

# Or manually run build scripts with debug output
bash -x scripts/pixi_install_rust_evm.sh
bash -x scripts/pixi_install_rust_pasa.sh
bash -x scripts/pixi_install_trinity_rust.sh
```

## Scripts That Handle Cleanup

The following scripts have been updated to handle failed installations gracefully:

- `scripts/pixi_install_rust_evm.sh`: Removes incomplete EVM source before cloning
- `scripts/pixi_install_rust_pasa.sh`: Removes incomplete PASA source before cloning
- `scripts/pixi_install_trinity_rust.sh`: Runs `cargo clean` before building
- `scripts/cleanup_rust_builds.sh`: Utility to manually clean installations
- `benchmark_trinity_rust.sbatch`: Commented-out cleanup lines (uncomment if needed)

## Quick Reference

```bash
# Check what needs cleaning
./scripts/cleanup_rust_builds.sh

# Clean everything
./scripts/cleanup_rust_builds.sh --force

# Rebuild
pixi install

# Verify
pixi shell
which sam_to_read_coords evidence_modeler pasa_rust
```

## Still Having Issues?

1. Check `TROUBLESHOOTING_INSTALLATIONS.md` (this file)
2. Look at detailed logs: `slurm_logs/*.err` or `pixi install -v`
3. Review `BENCHMARK_SLURM.md` for SLURM-specific issues
4. Check `BENCHMARK_TRINITY_RUST.md` for Trinity-specific issues

## Report a Bug

If you encounter an issue not covered here:

1. Run cleanup: `./scripts/cleanup_rust_builds.sh --force`
2. Rebuild with verbose output: `pixi install -v 2>&1 | tee install.log`
3. Check the log file for the actual error
4. Report with: `cat install.log` and `uname -a`
