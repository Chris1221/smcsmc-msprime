# Examples

This directory contains example scripts demonstrating how to use the smcsmc_msprime package.

## Demo Script

### [demo.py](demo.py)

A complete end-to-end demonstration of the particle filter:

1. **Simulates genetic data** under a known demographic model using msprime
2. **Runs the particle filter** on the simulated data
3. **Reports results** including log likelihood, ESS statistics, and resampling events

**Usage:**
```bash
python examples/demo.py
```

**What it does:**
- Simulates 20 haploid samples over 100kb with constant Ne=10,000
- Generates ~70 variant sites under mutation rate 1e-8
- Runs 50-particle filter with systematic resampling
- Calculates log likelihood using Felsenstein's pruning algorithm
- Reports particle filter statistics (ESS, resampling positions, etc.)

**Expected output:**
```
============================================================
SMCSMC-MSPRIME DEMO
============================================================

Simulating data...
  Generated ~70 variant sites

Running particle filter...
  50 particles, 20 samples, 72 sites
  Log Likelihood: -3862.13
  Mean ESS: 48.15
  Resampling events: 3

Demo complete!
```

## Running the Examples

Make sure you're in the pixi environment:

```bash
cd smcsmc-msprime
pixi shell
python examples/demo.py
```

## Modifying Parameters

You can edit [demo.py](demo.py) to test different scenarios:

```python
results = run_demo(
    num_particles=100,        # Try more particles
    num_samples=50,           # Try more samples
    sequence_length=1_000_000,  # Try longer sequences
    population_size=5_000,    # Try different Ne
    mutation_rate=2e-8,       # Try different mutation rate
    recombination_rate=2e-8,  # Try different recombination rate
)
```

## Creating Your Own Examples

Basic template for using the particle filter:

```python
import numpy as np
import msprime
import smcsmc_msprime as smc

# 1. Load or simulate your data
genotypes = ...  # Shape: (num_sites, num_samples)
positions = ...  # Variant positions

# 2. Define demographic model
demography = smc.constant_Ne(Ne=10000)
# Or: demography = smc.bottleneck(Ne_present=10000, Ne_bottleneck=1000, ...)
# Or: demography = smc.two_epoch(Ne_present=10000, Ne_ancestral=5000, change_time=1000)

# 3. Run particle filter
results = smc.run_particle_filter(
    genotypes=genotypes,
    positions=positions,
    num_particles=100,
    demography=demography,
    mutation_rate=1e-8,
    recombination_rate=1e-8,
    sequence_length=1e6,
    random_seed=42,
    verbose=True
)

# 4. Analyze results
print(f"Log likelihood: {results.log_likelihood}")
print(f"Mean ESS: {np.mean(results.ess_history)}")
print(f"Resampling events: {len(results.resample_positions)}")

# Access final particles
for i, particle in enumerate(results.final_particles[:3]):
    print(f"Particle {i}: {particle.num_trees} trees, weight={particle.weight}")
```

## Future Examples

Planned examples to be added:

- [ ] **Parameter estimation**: Using EM algorithm to infer Ne
- [ ] **Model comparison**: Comparing different demographic models
- [ ] **Real data**: Loading and analyzing VCF files
- [ ] **Validation**: Comparing results with original SMCSMC
- [ ] **Visualization**: Plotting particle weights, ESS, demographic trajectories
- [ ] **Multiple populations**: Migration and population structure
- [ ] **Performance benchmarking**: Scaling tests with various parameters
