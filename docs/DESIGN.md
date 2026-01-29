# SMCSMC-msprime: Prototype Design Document

## Overview

This is a Python prototype reimplementing the SMCSMC particle filter for demographic inference using **msprime** and **tskit** instead of SCRM. The goal is to achieve:

1. **Simplicity**: Pure Python implementation with clean APIs
2. **Scalability**: Handle 1000s of genomes instead of ~100
3. **Performance**: Leverage msprime's efficient tree sequence format
4. **Maintainability**: Use actively maintained libraries with good documentation

## Background: Current Implementation

### Architecture

The current SMCSMC is a C++/Python hybrid:

- **C++ backend** (~7,600 lines):
  - Particle filter core
  - Extends SCRM's `Forest` class for genealogy representation
  - Custom arena memory allocation
  - Hand-optimized likelihood calculations

- **SCRM integration**:
  - Git submodule providing coalescent simulation
  - Pointer-based tree structures
  - ~100 sample limit in practice

### Core Algorithm

SMCSMC implements a **Sequential Monte Carlo (SMC) particle filter** where:

1. Each **particle** represents a genealogical ARG (Ancestral Recombination Graph)
2. Particles are **propagated** along the genome by simulating recombination/coalescence
3. Particles are **weighted** by likelihood at each variant site
4. **Resampling** occurs when effective sample size (ESS) drops too low
5. **Parameter estimation** via EM/Variational Bayes on event counts

### Key Components

**ForestState (Particle)**:
- Extends SCRM's `Forest` class
- Members: `posterior_weight_`, `pilot_weight_`, `multiplicity_`
- Methods:
  - `extend_ARG()`: Simulate genealogy forward to new position
  - `calculate_likelihood()`: Recursive likelihood on tree
  - Delayed importance sampling for numerical stability

**ParticleContainer**:
- Maintains N particles
- `update_weight_at_site()`: Calculate likelihoods
- `resample()`: Systematic resampling when ESS < threshold
- Tracks normalization factors

**CountModel**:
- Extracts coalescence/recombination event counts
- Updates demographic parameters (Ne, migration rates)
- Implements online EM algorithm

### Computational Bottlenecks

1. **Tree manipulation**: Pointer-based ARG structures don't scale
2. **Memory**: Custom allocators needed to prevent fragmentation
3. **SCRM limitations**: Struggles with >100 samples
4. **Likelihood**: Recursive tree traversals (mitigated by fast exp approximations)

## New Approach: msprime/tskit

### Why msprime?

**msprime** is 100-1000× faster than SCRM and handles millions of samples:

1. **Succinct Tree Sequences**: O(sequence_length + mutations) not O(samples²)
2. **Highly optimized C backend** with Python API
3. **Active maintenance** and ecosystem (tskit, tsinfer, tsdate)
4. **Native ARG support** with efficient tree iteration
5. **Scales to biobank-sized datasets** (UK Biobank, All of Us)

### Architecture

```
┌─────────────────────────────────────────────────────┐
│                 Python Frontend                      │
│  (CLI, data loading, model specification, results)   │
└─────────────────┬───────────────────────────────────┘
                  │
┌─────────────────▼───────────────────────────────────┐
│            Particle Filter Core (Python)             │
│  • ParticleContainer: manage N particles             │
│  • Resampling: systematic resampling                 │
│  • Weight calculation & ESS monitoring               │
└─────────────────┬───────────────────────────────────┘
                  │
┌─────────────────▼───────────────────────────────────┐
│          Particle Class (TreeSequence wrapper)       │
│  • weight: float                                     │
│  • ts: tskit.TreeSequence                           │
│  • extend(): simulate forward                        │
│  • likelihood(): compute on tree                     │
└─────────────────┬───────────────────────────────────┘
                  │
┌─────────────────▼───────────────────────────────────┐
│              msprime + tskit (C/Python)              │
│  • Coalescent simulation (msprime.sim_ancestry)      │
│  • Tree sequence storage (tskit)                     │
│  • Efficient tree iteration                          │
└──────────────────────────────────────────────────────┘
```

### Key Design Decisions

**1. Particle Representation**

**Current (SCRM)**:
```cpp
class ForestState : public Forest {
    double posterior_weight_;
    double pilot_weight_;
    int multiplicity_;
    // ... complex pointer-based tree structure
}
```

**New (msprime)**:
```python
class Particle:
    def __init__(self, ts: tskit.TreeSequence):
        self.ts = ts              # Tree sequence (ARG)
        self.weight = 1.0          # Posterior weight
        self.pilot_weight = 1.0    # Pilot weight for resampling
        self.multiplicity = 1      # Number of copies
```

**Advantages**:
- Immutable tree sequences (copy-on-write)
- Built-in serialization
- Efficient iteration over local trees
- Metadata support

**2. ARG Extension**

**Current approach**:
- Call SCRM's `sampleNextBase()` to extend genealogy
- Record recombination/coalescence events
- Apply importance weights for biased sampling

**New approach**:
```python
def extend_particle(particle, start_pos, end_pos, demography):
    """Extend particle's ARG from start_pos to end_pos."""

    # Use msprime to simulate new genealogy segment
    # with recombination starting from current state
    new_ts = msprime.sim_ancestry(
        initial_state=particle.ts,
        start_time=particle.ts.max_root_time,
        demography=demography,
        recombination_rate=...,
        sequence_length=end_pos - start_pos,
        # ... other params
    )

    # Calculate importance weight for biased sampling
    importance_weight = calculate_extension_weight(
        particle.ts, new_ts, start_pos, end_pos
    )

    particle.ts = new_ts
    particle.weight *= importance_weight
    return particle
```

**Challenges**:
- msprime expects to simulate from scratch, not extend existing ARG
- May need to use `initial_state` parameter or build custom extension
- Importance sampling weights need careful calculation

**3. Likelihood Calculation**

**Current approach**:
- Recursive function traversing SCRM's tree
- Fast exponential approximations for mutation probability
- Handles unphased data by summing over configurations

**New approach**:
```python
def calculate_likelihood(ts, haplotypes, mutation_rate):
    """
    Calculate P(data | genealogy) using tree sequence.

    Uses Felsenstein's pruning algorithm on each local tree.
    """
    log_likelihood = 0.0

    for tree in ts.trees():
        # For each variant in this tree's span
        for variant in ts.variants(left=tree.interval.left,
                                     right=tree.interval.right):
            # Felsenstein pruning from leaves to root
            prob = felsenstein_pruning(
                tree,
                variant.genotypes,
                mutation_rate,
                tree.interval.right - tree.interval.left
            )
            log_likelihood += np.log(prob)

    return np.exp(log_likelihood)
```

**Advantages**:
- tskit provides efficient tree iteration
- Can access branch lengths directly
- Natural handling of variant sites

**4. Resampling**

**Current approach**:
- Systematic resampling when ESS < threshold
- Particles have `multiplicity` field
- Copy constructor for particle duplication

**New approach**:
```python
def systematic_resampling(particles, weights):
    """
    Resample particles using systematic resampling.
    Returns new list of particles.
    """
    N = len(particles)
    weights = np.array(weights) / np.sum(weights)

    # Systematic resampling positions
    positions = (np.arange(N) + np.random.uniform()) / N
    cumsum = np.cumsum(weights)

    new_particles = []
    j = 0
    for pos in positions:
        while cumsum[j] < pos:
            j += 1
        # Copy particle (tree sequence is immutable, cheap to copy)
        new_particles.append(particles[j].copy())

    return new_particles
```

**Advantages**:
- Simpler than multiplicity tracking
- Tree sequences are cheap to copy (copy-on-write)
- Standard algorithm implementation

**5. Parameter Estimation**

**Current approach**:
- Extract event counts from each particle's `eventTrees`
- Sum across particles weighted by importance
- EM/VB updates for Ne, migration, recombination rates

**New approach**:
```python
def extract_event_counts(particles, weights):
    """Extract coalescence events from particle ensemble."""

    event_counts = defaultdict(float)  # {(epoch, event_type): count}

    for particle, weight in zip(particles, weights):
        for edge in particle.ts.edges():
            # Identify which epoch this coalescence occurred in
            time = particle.ts.node(edge.parent).time
            epoch = get_epoch(time, change_times)

            event_counts[(epoch, 'coal')] += weight

            # Similarly for recombination (track edges that split)
            # ...

    return event_counts
```

**Advantages**:
- Direct access to edge table
- Times and topologies readily available
- Can use tskit statistics functions

## Implementation Plan

### Phase 1: Core Components (Week 1)

**Files to create**:
1. `src/smcsmc/particle.py`: Particle class wrapping TreeSequence
2. `src/smcsmc/likelihood.py`: Felsenstein pruning on trees
3. `src/smcsmc/resampling.py`: Systematic resampling
4. `src/smcsmc/demography.py`: Demographic model specification

**Deliverable**: Can initialize particles, calculate likelihoods, resample

### Phase 2: Particle Filter Loop (Week 2)

**Files to create**:
5. `src/smcsmc/particle_filter.py`: Main SMC loop
6. `src/smcsmc/extension.py`: ARG extension logic

**Deliverable**: Can run particle filter along a chromosome

### Phase 3: Parameter Estimation (Week 3)

**Files to create**:
7. `src/smcsmc/inference.py`: Event counting and EM updates
8. `src/smcsmc/utils.py`: Helper functions

**Deliverable**: Can estimate Ne from simulated data

### Phase 4: Testing & Validation (Week 4)

**Files to create**:
9. `tests/test_*.py`: Unit tests for each module
10. `examples/simple_demo.py`: End-to-end example

**Deliverable**: Validated against known demographic scenarios

## Current Prototype Status

### Completed
- ✅ Pixi project setup with dependencies (msprime, tskit, numpy, scipy)
- ✅ Package structure created (`src/smcsmc/`)

### In Progress
- 🔄 Particle class implementation

### TODO
- ⬜ Likelihood calculation
- ⬜ Resampling logic
- ⬜ Particle filter loop
- ⬜ ARG extension
- ⬜ Parameter estimation
- ⬜ Demo script

## Challenges & Open Questions

### 1. ARG Extension Strategy

**Problem**: msprime simulates from scratch; we need to extend existing ARGs

**Options**:
- **A**: Use `initial_state` parameter (if it supports mid-chromosome starts)
- **B**: Extract final state, use as starting point for next segment
- **C**: Simulate entire chromosome, then "observe" it incrementally
- **D**: Modify msprime or use lower-level API

**Recommendation**: Start with **C** (simplest) for prototype, optimize later

### 2. Importance Sampling

**Problem**: Need to calculate importance weights for biased sampling

**Current approach**:
- Weight by recombination rate to focus samples on high-recombination regions
- Delayed application of weights for numerical stability

**New approach**:
- Calculate proposal vs. target density ratio
- May be simpler with msprime's explicit simulation

### 3. Unphased Data

**Problem**: Current implementation sums over haplotype configurations

**Options**:
- **A**: Require phased input (simplest)
- **B**: Implement configuration enumeration
- **C**: Use phase-aware likelihood (Li & Stephens style)

**Recommendation**: **A** for prototype, **C** for production

### 4. Memory Efficiency

**Problem**: Tree sequences can be large for many samples × long sequences

**Solutions**:
- Use tskit's simplification to remove ancient nodes
- Only store recent k coalescences
- Checkpoint to disk periodically

## Performance Expectations

### Current Implementation
- **Samples**: ~100 maximum practical
- **Sequence length**: Whole chromosome (100s Mb)
- **Particles**: 100-1000
- **Runtime**: Hours to days per chromosome

### Expected with msprime
- **Samples**: 1000-10000 (100× improvement)
- **Sequence length**: Whole genome
- **Particles**: 1000-10000
- **Runtime**: Minutes to hours (msprime is 100-1000× faster)

**Caveat**: Initial Python prototype will be slower than optimized C++.
After validation, can:
- Profile and optimize hotspots
- Use Numba JIT compilation
- Rewrite critical sections in Rust/C++
- Parallelize across particles

## Next Steps

1. Implement `Particle` class with basic operations
2. Implement likelihood calculation on tree sequences
3. Implement resampling
4. Create simple demo showing particle filter on simulated data
5. Validate against known demographic scenarios
6. Benchmark and profile
7. Optimize based on profiling results

## References

- **Current SMCSMC**: Henderson et al. (methodology papers)
- **msprime**: Kelleher et al. 2016, PLOS Comp Bio
- **tskit**: Kelleher et al. 2019, Genetics
- **SCRM**: Staab et al. 2015, Bioinformatics
- **SMC algorithms**: Doucet & Johansen 2009 (tutorial)
