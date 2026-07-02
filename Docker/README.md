# Funannotate Docker Images

Dockerfiles for building and running funannotate with Rust-optimized PASA and EVidenceModeler components.

## Available Images

### 1. Standard Image (Dockerfile)
**Size**: ~3-4GB  
**Use case**: Full funannotate annotation pipeline (recommended for most users)  
**Features**:
- All bioconda dependencies
- Rust-optimized PASApipeline (rust_optimize branch, transcript assembly)
- Rust-optimized EVidenceModeler (rust_optimize branch, consensus modeling)
- Trinity (RNA-seq assembly)
- Augustus, GeneMark, SNAP, GlimmerHMM (gene predictors)
- BLAST, DIAMOND, HMMER (protein searches)
- Minimap2, HISAT2, Samtools (sequence alignment)
- No pre-installed databases (you can mount or download)

**Build**:
```bash
cd Docker
./build.sh                           # Builds as funannotate:latest
./build.sh Dockerfile funannotate:1.9-rust
docker build -f Dockerfile -t my-registry/funannotate:latest .
```

**Run**:
```bash
# Get help
docker run -it funannotate:latest funannotate --help

# Annotation workflow (mount local data)
docker run -it -v $(pwd)/data:/work funannotate:latest \
  funannotate predict -i genome.fasta -o output/

# Run interactively
docker run -it -v $(pwd):/work funannotate:latest /bin/bash
```

### 2. Development Image (Dockerfile.dev)
**Size**: ~3-4GB  
**Use case**: Development with local PASApipeline and EVidenceModeler (rust_optimize branch) checkouts  
**Features**:
- Same as standard image, but built from local source
- Uses `../PASApipeline` and `../EVidenceModeler` from your filesystem
- Ideal for testing changes before committing
- Includes git for development

**Build** (requires local checkouts as siblings of `funannotate-live/`,
i.e. the *parent* of this repo is the Docker build context so `COPY
trinityrnaseq`/`PASApipeline`/`EVidenceModeler` in Dockerfile.dev can see
them):
```bash
# From the parent directory of funannotate-live, with sibling checkouts of
# trinityrnaseq, PASApipeline, EVidenceModeler (rust_optimize branch)
docker build -f funannotate-live/Dockerfile.dev -t funannotate:dev .
```

**Requirements** (all three, checked out to the `rust_optimize` branch, as
siblings of `funannotate-live/`):
- `../trinityrnaseq`
- `../PASApipeline`
- `../EVidenceModeler`

### 3. Slim with Databases (Dockerfile2)
**Size**: ~10-15GB (includes all BUSCO lineage databases)  
**Use case**: Production deployment with pre-installed databases  
**Features**:
- Extends from funannotate:latest
- Pre-installed BUSCO databases for all lineages
- Pre-installed InterPro, Pfam, dbCAN references
- Ready to run without additional setup
- Larger but faster for users

**Build**:
```bash
cd Docker
./build.sh Dockerfile2 funannotate:with-databases

# Or directly
docker build -f Dockerfile2 -t funannotate:with-databases .
```

**Run**:
```bash
docker run -it -v $(pwd):/work funannotate:with-databases \
  funannotate predict -i genome.fasta -o output/ \
  -s "Species name" -d /opt/databases
```

## Building Images

### Quick Build
```bash
cd funannotate-live/Docker
./build.sh                    # Build default (Dockerfile -> funannotate:latest)
```

### Custom Builds
```bash
# Build with custom tag
./build.sh Dockerfile funannotate:1.9.0-rust

# Build with custom registry
./build.sh Dockerfile myregistry.com/funannotate:latest

# Build with cache
./build.sh Dockerfile funannotate:latest --cache
```

### Build Script Options
```bash
./build.sh [DOCKERFILE] [TAG] [EXTRA_ARGS]

# Examples:
./build.sh                                    # Default: Dockerfile -> funannotate:latest
./build.sh Dockerfile                        # Same as default
./build.sh Dockerfile2 funannotate:db        # Build slim with tag
./build.sh Dockerfile funannotate:dev --cache  # Build with cache
```

## Dockerfile Breakdown

### Dockerfile (Standard)

**Stage 1: Build**
1. Start from `continuumio/miniconda3`
2. Install mamba, conda-pack
3. Create conda environment with all bioconda dependencies
4. Install funannotate Python package from GitHub
5. Clone and build EVidenceModeler (rust_optimize branch) from source
6. Clone and build PASApipeline (rust_optimize branch) from source (C++ + Rust components)
7. Pack environment with conda-pack

**Stage 2: Runtime**
1. Start from `debian:bullseye`
2. Copy conda environment from build stage
3. Install system packages (Augustus, Snap, Perl, SQLite)
4. Fix symlinks
5. Set environment variables
6. Create non-root user
7. Verify all tools are installed

**Multi-stage benefits**:
- Build artifacts not included in final image
- Smaller final image size
- Faster deployment from registry
- Reproducible builds

### Dockerfile.dev (Development)

Similar to standard Dockerfile but:
- Copies local `EVidenceModeler` and `PASApipeline` (rust_optimize branch) from host
- Uses local source instead of cloning from GitHub
- Includes git in runtime image

### Dockerfile2 (Database Image)

Extends the standard image with:
- Runs `funannotate setup` to download all BUSCO lineage databases
- Pre-installs InterPro, Pfam, dbCAN, and other reference databases
- Much larger (~10GB) but eliminates database download step at runtime

## Key Environment Variables

```bash
PATH="/venv/bin:/venv/opt/pasa-rust/bin:$PATH"
AUGUSTUS_CONFIG_PATH="/usr/share/augustus/config"
EVM_HOME="/venv/opt/evm-rust/src"
PASAHOME="/venv/opt/pasa-rust/src"
FUNANNOTATE_EVM_ENGINE="rust"           # Use Rust implementation
TRINITYHOME="/venv/bin"
TRINITY_HOME="/venv/bin"
QUARRY_PATH="/venv/opt/codingquarry-2.0/QuarryFiles"
ZOE="/usr/share/snap"
FUNANNOTATE_DB="/opt/databases"         # Database location (mount or download)
```

## Example Workflows

### Basic Annotation
```bash
docker run -it \
  -v $(pwd)/data:/work \
  funannotate:latest \
  funannotate predict \
    -i genome.fasta \
    -o output/ \
    -s "Species name"
```

### With RNA-seq Training
```bash
docker run -it \
  -v $(pwd)/data:/work \
  funannotate:latest \
  bash -c "
    funannotate train \
      -i genome.fasta \
      -o training_output/ \
      --left reads_R1.fastq.gz \
      --right reads_R2.fastq.gz \
      -s 'Species name' && \
    funannotate predict \
      -i genome.fasta \
      -o annotation_output/ \
      -s 'Species name'
  "
```

### Using Pre-installed Databases
```bash
docker run -it \
  -v $(pwd)/data:/work \
  funannotate:with-databases \
  funannotate annotate \
    -i annotation_output/ \
    -s "Species name" \
    -d /opt/databases
```

### Interactive Shell
```bash
docker run -it \
  -v $(pwd)/data:/work \
  --rm \
  funannotate:latest \
  /bin/bash
```

### Check Dependencies
```bash
docker run -it --rm funannotate:latest \
  funannotate check --show-versions
```

## Docker Compose Example

```yaml
version: '3.8'

services:
  funannotate:
    image: funannotate:latest
    volumes:
      - ./data:/work
      - ./output:/output
    environment:
      FUNANNOTATE_DB: /work/databases
    working_dir: /work
    command: >
      funannotate predict
      -i genome.fasta
      -o /output/annotations
      -s "Species name"
    user: "1000:1000"

  # Optional: separate service for training
  funannotate-train:
    image: funannotate:latest
    volumes:
      - ./data:/work
      - ./training:/training
    working_dir: /work
    command: >
      funannotate train
      -i genome.fasta
      -o /training/output
      --left reads_R1.fastq.gz
      --right reads_R2.fastq.gz
      -s "Species name"
```

## Troubleshooting

### Docker Build Fails

**Insufficient disk space**:
```bash
docker system prune -a  # Remove all unused images/containers
df -h                   # Check available space (need ~20GB)
```

**Network issues during clone**:
```bash
docker build --build-arg http_proxy=http://proxy:port ... -f Dockerfile .
```

**Build cache issues**:
```bash
docker build --no-cache -f Dockerfile -t funannotate:latest .
```

### Docker Run Issues

**Permission denied writing to volume**:
```bash
# Run with user ID
docker run -u $(id -u):$(id -g) -v $(pwd):/work funannotate:latest ...

# Or use sudo
sudo docker run -v $(pwd):/work funannotate:latest ...
```

**Out of memory**:
```bash
# Increase Docker desktop memory (Settings -> Resources)
# Or increase swap space on Linux
```

**Cannot find database**:
```bash
# Use database image
docker run -v $(pwd):/work funannotate:with-databases ...

# Or download databases to host
funannotate setup -d /path/to/db
docker run -v /path/to/db:/opt/databases funannotate:latest ...
```

## Publishing Images

To publish to Docker Hub or a private registry:

```bash
# Tag image
docker tag funannotate:latest myregistry/funannotate:1.9-rust
docker tag funannotate:latest myregistry/funannotate:latest

# Push to registry
docker push myregistry/funannotate:1.9-rust
docker push myregistry/funannotate:latest

# For Docker Hub
docker login docker.io
docker tag funannotate:latest username/funannotate:latest
docker push username/funannotate:latest
```

## Notes

- All images are built for **Linux x86_64 only**
- Build time: 20-40 minutes depending on system (includes Rust compilation)
- Final images are ~3-4GB (or ~10-15GB with databases)
- Multi-stage build keeps final image small
- Non-root user `funannotate` (UID 1000) is used by default for security
- Verify installations automatically during build (warnings are OK)

## References

- [funannotate repository](https://github.com/nextgenusfs/funannotate)
- [PASApipeline (rust_optimize branch) repository](https://github.com/hyphaltip/PASApipeline)
- [EVidenceModeler (rust_optimize branch) repository](https://github.com/hyphaltip/EVidenceModeler)
- [Docker documentation](https://docs.docker.com/)
- [Docker Compose documentation](https://docs.docker.com/compose/)
