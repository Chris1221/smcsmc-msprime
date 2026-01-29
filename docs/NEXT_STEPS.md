# Next Steps for Implementation

## Current Progress

### ✅ Completed
1. Analyzed original SMCSMC implementation (C++/SCRM)
2. Identified performance bottlenecks and scalability issues
3. Designed msprime-based architecture
4. Set up pixi project with dependencies
5. Created documentation (DESIGN.md, README.md)

### 🔄 In Progress
- Implementing core Particle class

### ⬜ TODO
- Likelihood calculation module
- Resampling logic
- Particle filter loop
- Demo script

## Directory Structure

```
smcsmc-msprime/
├── pyproject.toml          # Pixi configuration
├── README.md               # Project overview
├── DESIGN.md              # Detailed design document
├── NEXT_STEPS.md          # This file
├── src/
│   └── smcsmc/
│       ├── __init__.py    # Package init
│       ├── particle.py    # Particle class (TreeSequence wrapper)
│       ├── likelihood.py  # Felsenstein pruning
│       ├── resampling.py  # Systematic resampling
│       ├── particle_filter.py  # Main SMC loop
│       ├── extension.py   # ARG extension logic
│       ├── demography.py  # Demographic models
│       └── inference.py   # Parameter estimation
├── tests/
│   └── test_*.py          # Unit tests
└── examples/
    └── demo.py            # Simple demo script
```

## Implementation Order

### Phase 1: Core Components (Priority 1)

**1. Particle Class** (`src/smcsmc/particle.py`)
```python
class Particle:
    """Wrapper around tskit.TreeSequence for particle filtering."""

    def __init__(self, ts: tskit.TreeSequence):
        self.ts = ts
        self.weight = 1.0
        self.pilot_weight = 1.0
        self.multiplicity = 1

    def copy(self):
        """Create a copy of this particle."""
        pass

    def likelihood(self, genotypes, mutation_rate):
        """Calculate P(data | genealogy)."""
        pass
```

**Key methods needed**:
- `__init__`: Initialize with tree sequence
- `copy`: For resampling (cheap with immutable ts)
- `likelihood`: Calculate data likelihood given genealogy

**2. Likelihood Calculation** (`src/smcsmc/likelihood.py`)
```python
def felsenstein_pruning(tree, genotypes, mutation_rate, span):
    """
    Calculate likelihood using Felsenstein's pruning algorithm.

    Args:
        tree: tskit.Tree object
        genotypes: Array of observed genotypes at leaves
        mutation_rate: Per-base mutation rate
        span: Genomic span of this tree

    Returns:
        float: Likelihood P(genotypes | tree)
    """
    pass
```

**Key algorithm**: Recursively calculate likelihoods from leaves to root

**3. Resampling** (`src/smcsmc/resampling.py`)
```python
def systematic_resampling(particles, weights):
    """
    Systematic resampling of particles.

    Args:
        particles: List of Particle objects
        weights: Array of (normalized) weights

    Returns:
        List of resampled Particle objects
    """
    pass

def effective_sample_size(weights):
    """Calculate ESS = 1 / sum(w_i^2)."""
    pass
```

### Phase 2: Particle Filter (Priority 2)

**4. Particle Filter Loop** (`src/smcsmc/particle_filter.py`)
```python
class ParticleFilter:
    def __init__(self, num_particles, demography, data):
        self.num_particles = num_particles
        self.demography = demography
        self.data = data
        self.particles = []

    def initialize_particles(self):
        """Create initial particle ensemble."""
        pass

    def propagate(self, start_pos, end_pos):
        """Propagate particles from start_pos to end_pos."""
        pass

    def update_weights(self, position):
        """Update particle weights at variant site."""
        pass

    def resample_if_needed(self, threshold=0.5):
        """Resample if ESS < threshold * num_particles."""
        pass

    def run(self):
        """Run full particle filter along chromosome."""
        pass
```

**5. ARG Extension** (`src/smcsmc/extension.py`)
```python
def extend_particle(particle, start_pos, end_pos, demography, recomb_rate):
    """
    Extend particle's ARG from start_pos to end_pos.

    This is the tricky part - need to figure out how to use msprime
    to extend an existing tree sequence.

    Options:
    A. Simulate whole chromosome upfront, observe incrementally
    B. Use msprime initial_state parameter
    C. Custom extension logic

    Returns:
        Particle: Extended particle with updated weight
    """
    pass
```

### Phase 3: Demo & Testing (Priority 3)

**6. Demo Script** (`examples/demo.py`)
```python
import msprime
import numpy as np
from smcsmc import ParticleFilter, Demography

# Simulate data
ts = msprime.sim_ancestry(
    samples=10,
    population_size=10000,
    sequence_length=1e6,
    recombination_rate=1e-8,
    random_seed=42
)
ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=42)

# Extract genotypes
genotypes = ts.genotype_matrix()

# Run particle filter
pf = ParticleFilter(
    num_particles=100,
    demography=Demography.constant_Ne(Ne=10000),
    data=genotypes
)

results = pf.run()
print(f"Estimated Ne: {results.Ne}")
```

## Key Decisions Needed

### 1. ARG Extension Strategy
**Question**: How to extend tree sequences incrementally?

**Options**:
- **Simple (for prototype)**: Simulate entire chromosome, observe incrementally
  - Pros: Easy to implement
  - Cons: Not realistic (knows future), memory intensive

- **Realistic**: Use msprime's `initial_state` or custom extension
  - Pros: Correct importance sampling
  - Cons: More complex, may need msprime modifications

**Recommendation**: Start with simple approach, validate algorithm, then add realistic extension

### 2. Likelihood Model
**Question**: Infinite sites vs. finite sites mutation model?

**Current SMCSMC**: Uses fast exponential approximations

**Options**:
- **Infinite sites**: Assumes each mutation is unique (standard for coalescent)
- **Finite sites**: Allows recurrent mutation (more realistic but slower)

**Recommendation**: Start with infinite sites for simplicity

### 3. Phasing
**Question**: Handle phased or unphased data?

**Current SMCSMC**: Sums over haplotype configurations for unphased data

**Options**:
- **Phased only**: Simplest
- **Sum over configs**: Match current implementation
- **Statistical phasing**: Use Li & Stephens or similar

**Recommendation**: Phased only for prototype

## Testing Strategy

### Unit Tests
- `test_particle.py`: Test Particle class operations
- `test_likelihood.py`: Validate against known tree/genotype combos
- `test_resampling.py`: Check resampling produces correct distribution
- `test_particle_filter.py`: Integration test on simulated data

### Validation Tests
1. **Constant Ne**: Should recover known constant population size
2. **Exponential growth**: Should detect growth rate
3. **Bottleneck**: Should identify bottleneck timing and severity
4. **Compare to original**: Match SMCSMC results on small datasets

### Performance Benchmarks
- Time to process 100 samples × 1 Mb
- Memory usage scaling with samples
- Compare to original SMCSMC (should be comparable or better)

## Questions for Discussion

1. **Do we want to match original SMCSMC exactly** (for validation)?
   - If yes: Need same importance sampling scheme
   - If no: Can simplify some aspects

2. **What demographic models to support initially?**
   - Start with constant Ne?
   - Add piecewise constant Ne?
   - Migration?

3. **Performance targets?**
   - How many samples?
   - How long is acceptable runtime?

4. **Do we need the full EM algorithm** or just forward filter?
   - Just filtering: Easier to implement
   - Full EM: Needed for actual inference

## Resources

- **msprime docs**: https://tskit.dev/msprime/docs/stable/
- **tskit docs**: https://tskit.dev/tskit/docs/stable/
- **Felsenstein pruning**: Classic phylogenetics algorithm
- **SMC tutorial**: Doucet & Johansen 2009

## When Resuming Work

1. Read DESIGN.md for full context
2. Start with implementing `particle.py`
3. Then `likelihood.py` (most algorithmically interesting)
4. Then `resampling.py` (straightforward)
5. Build up to `particle_filter.py`
6. Create `demo.py` to test end-to-end

The key challenge is the ARG extension - everything else is relatively standard SMC.
