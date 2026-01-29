# SMCSMC-msprime

A Python prototype reimplementing the SMCSMC particle filter for demographic inference using **msprime** and **tskit**.

## Goals

1. **Scalability**: Handle 1000s of genomes (vs. ~100 in current implementation)
2. **Performance**: Leverage msprime's efficient tree sequence format (100-1000× faster than SCRM)

## Installation

Using pixi (recommended):

```bash
cd new/smcsmc-msprime
pixi install
pixi shell
```

Dependencies:
- msprime >=1.3.4
- tskit >=1.0.0
- numpy >=2.4
- scipy >=1.17

## Architecture

```
Particle Filter
    ├── Particle (wraps tskit.TreeSequence)
    │   ├── weight (importance weight)
    │   ├── ts (genealogy/ARG)
    │   └── methods: extend(), likelihood()
    │
    ├── ParticleContainer
    │   ├── propagate particles along genome
    │   ├── calculate weights (likelihood)
    │   └── resample when ESS low
    │
    └── Inference
        ├── extract event counts
        └── EM/VB parameter updates
```

## Key Components (from DESIGN.md)

**Particle**: Wraps a `tskit.TreeSequence` with:
- Genealogical ARG (ancestral recombination graph)
- Importance weight
- Methods for extension and likelihood calculation

**Likelihood**: Felsenstein's pruning algorithm on tree sequence

**Resampling**: Systematic resampling when effective sample size drops

**Extension**: Simulate ARG forward using msprime

**Inference**: Extract coalescence/recombination events, update demographic parameters

## Quick Start (Coming Soon)

```python
from smcsmc import ParticleFilter, Demography

# Define demographic model
demography = Demography.constant_Ne(Ne=10000)

# Load data
data = load_vcf("data.vcf")

# Run particle filter
pf = ParticleFilter(
    num_particles=1000,
    demography=demography,
    data=data
)

# Infer parameters
results = pf.run()
print(f"Estimated Ne: {results.Ne}")
```