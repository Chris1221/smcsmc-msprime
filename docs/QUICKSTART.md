# Quick Start Guide

Get up and running with smcsmc-msprime in 5 minutes!

## Installation

```bash
cd smcsmc-msprime
pixi install
pixi shell
```

## Run the Demo

```bash
python examples/demo.py
```

Expected output: Particle filter runs successfully on simulated data and reports log likelihood!

## Your First Particle Filter

Create a file `my_analysis.py`:

```python
import numpy as np
import msprime
import smcsmc_msprime as smc

# 1. Simulate some data (or load your own)
print("Simulating data...")
ts = msprime.sim_ancestry(
    samples=20,
    population_size=10_000,
    sequence_length=100_000,
    recombination_rate=1e-8,
    random_seed=42,
    ploidy=1
)
ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=42)

# Extract genotypes and positions
genotypes = ts.genotype_matrix()
positions = np.array([site.position for site in ts.sites()])

print(f"Generated {ts.num_sites} variant sites")

# 2. Create a demographic model
demography = smc.constant_Ne(Ne=10_000)

# 3. Run the particle filter
print("Running particle filter...")
results = smc.run_particle_filter(
    genotypes=genotypes,
    positions=positions,
    num_particles=50,
    demography=demography,
    mutation_rate=1e-8,
    recombination_rate=1e-8,
    sequence_length=100_000,
    random_seed=123,
    verbose=True
)

# 4. View results
print(f"\nLog Likelihood: {results.log_likelihood:.2f}")
print(f"Mean ESS: {np.mean(results.ess_history):.1f}")
print(f"Resampling events: {len(results.resample_positions)}")

# Access final particles
print("\nFinal particles:")
for i, p in enumerate(results.final_particles[:3]):
    print(f"  Particle {i}: {p.num_trees} trees, weight={p.weight:.6f}")
```

Run it:
```bash
python my_analysis.py
```

## Different Demographic Models

### Constant Population Size
```python
demography = smc.constant_Ne(Ne=10_000)
```

### Two-Epoch Model (size change)
```python
demography = smc.two_epoch(
    Ne_present=10_000,    # Recent size
    Ne_ancestral=5_000,   # Ancient size
    change_time=1_000     # When it changed (generations ago)
)
```

### Bottleneck
```python
demography = smc.bottleneck(
    Ne_present=10_000,
    Ne_bottleneck=1_000,
    bottleneck_start=500,   # Started 500 gen ago
    bottleneck_end=1_000    # Ended 1000 gen ago
)
```

### Exponential Growth
```python
demography = smc.exponential(
    Ne_present=10_000,
    growth_rate=0.01  # Positive = growth, negative = decline
)
```

### Piecewise Constant
```python
demography = smc.DemographyHelper.piecewise_constant(
    sizes=[10_000, 5_000, 20_000],  # Most recent first
    change_times=[1_000, 5_000]     # Times of changes
)
```

## Adjusting Particle Filter Settings

```python
from smcsmc_msprime import ParticleFilterConfig

config = ParticleFilterConfig(
    num_particles=100,              # More particles = better, but slower
    resample_threshold=0.5,         # Resample when ESS < 0.5*N
    mutation_rate=1e-8,
    recombination_rate=1e-8,
    sequence_length=1e6,
    random_seed=42,
    verbose=True                    # Print progress
)

# Use the config
from smcsmc_msprime import ParticleFilter

pf = ParticleFilter(
    num_samples=20,
    genotypes=genotypes,
    positions=positions,
    demography=demography,
    config=config
)

results = pf.run()
```

## Working with Results

```python
# Log likelihood
print(f"Total log likelihood: {results.log_likelihood}")
print(f"Per-site log likelihood: {results.log_likelihood / len(positions)}")

# ESS history
import matplotlib.pyplot as plt
plt.plot(results.positions, results.ess_history)
plt.xlabel("Position")
plt.ylabel("Effective Sample Size")
plt.title("ESS Along the Genome")
plt.show()

# Resampling positions
print(f"Resampled at positions: {results.resample_positions}")

# Final particles
for i, particle in enumerate(results.final_particles):
    print(f"Particle {i}:")
    print(f"  Trees: {particle.num_trees}")
    print(f"  Weight: {particle.weight}")
    print(f"  Sequence length: {particle.sequence_length}")
```

## Loading Your Own Data

If you have genotype data:

```python
# Your genotype matrix: rows=sites, columns=samples
# Values should be 0 or 1 (biallelic SNPs)
genotypes = np.array([
    [0, 0, 1, 1, 0, 1, ...],  # Site 1
    [1, 0, 0, 1, 0, 0, ...],  # Site 2
    ...
])

# Physical positions of each site (in base pairs)
positions = np.array([1000, 2500, 5200, ...])

# Run the filter
results = smc.run_particle_filter(
    genotypes=genotypes,
    positions=positions,
    num_particles=100,
    demography=your_demographic_model,
    mutation_rate=1e-8,
    recombination_rate=1e-8,
    sequence_length=positions[-1] * 1.1,  # Slightly longer than last position
)
```

## Troubleshooting

### "All particles have zero weight"
- Your demographic model might be very unlikely given the data
- Try a different demographic model
- Check that mutation/recombination rates are reasonable

### Out of memory
- Reduce `num_particles`
- Process shorter sequences
- Use `sequence_length` that matches your actual data

### Slow performance
- Reduce `num_particles`
- Process fewer sites
- Consider downsampling sites with low information

## Next Steps

- Read [DESIGN.md](DESIGN.md) for architecture details
- Read [IMPLEMENTATION.md](IMPLEMENTATION.md) for what's been built
- Read [STATUS.md](STATUS.md) for current project status
- Check [examples/README.md](examples/README.md) for more examples
- Read [NEXT_STEPS.md](NEXT_STEPS.md) for future development

## Getting Help

- Check the docstrings: `help(smc.ParticleFilter)`
- Read the source code in `src/smcsmc_msprime/`
- Look at [examples/demo.py](examples/demo.py)

## Key Concepts

### Particle
A particle represents one possible genealogical history (ARG) for your samples. The particle filter maintains an ensemble of weighted particles.

### Tree Sequence
The underlying data structure from tskit that efficiently represents genealogies. Each particle contains a tree sequence.

### Likelihood
How well a particle's genealogy explains the observed genetic data. Calculated using Felsenstein's pruning algorithm.

### Resampling
When most particles have very low weight, we "resample" - duplicate high-weight particles and discard low-weight ones. This prevents particle degeneracy.

### ESS (Effective Sample Size)
Measures particle diversity. High ESS = diverse particles (good). Low ESS = most weight on few particles (bad, triggers resampling).

## Happy filtering! 🎉
