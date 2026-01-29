#!/usr/bin/env python3
"""
Simple demonstration of SMCSMC-msprime particle filter.

This script:
1. Simulates genetic data under a known demographic model
2. Runs the particle filter on the simulated data
3. Reports the likelihood and particle filter statistics

Usage:
    python examples/demo.py
"""

import numpy as np
import msprime
import sys

# Import smcsmc_msprime package
import smcsmc_msprime as smc


def simulate_data(
    num_samples=20,
    sequence_length=100_000,
    population_size=10_000,
    mutation_rate=1e-8,
    recombination_rate=1e-8,
    random_seed=42
):
    """
    Simulate genetic data under a simple demographic model.

    Args:
        num_samples: Number of haploid samples
        sequence_length: Length of sequence
        population_size: Effective population size
        mutation_rate: Mutation rate per base per generation
        recombination_rate: Recombination rate per base per generation
        random_seed: Random seed for reproducibility

    Returns:
        Tuple of (genotypes, positions, tree_sequence)
    """
    print("=" * 60)
    print("SIMULATING DATA")
    print("=" * 60)
    print(f"Samples: {num_samples}")
    print(f"Sequence length: {sequence_length:,} bp")
    print(f"Population size (Ne): {population_size:,}")
    print(f"Mutation rate: {mutation_rate:.2e}")
    print(f"Recombination rate: {recombination_rate:.2e}")
    print()

    # Simulate tree sequence
    print("Simulating ancestry...")
    ts = msprime.sim_ancestry(
        samples=num_samples,
        population_size=population_size,
        sequence_length=sequence_length,
        recombination_rate=recombination_rate,
        random_seed=random_seed,
        ploidy=1  # Haploid
    )

    print(f"  Generated {ts.num_trees} local trees")
    print(f"  Total nodes: {ts.num_nodes}")
    print(f"  Total edges: {ts.num_edges}")
    print()

    # Add mutations
    print("Adding mutations...")
    ts = msprime.sim_mutations(
        ts,
        rate=mutation_rate,
        random_seed=random_seed
    )

    print(f"  Generated {ts.num_sites} variant sites")
    print(f"  Total mutations: {ts.num_mutations}")
    print()

    # Extract genotypes and positions
    genotypes = ts.genotype_matrix()
    positions = np.array([site.position for site in ts.sites()])

    print(f"Genotype matrix shape: {genotypes.shape}")
    print(f"Position range: {positions[0]:.0f} - {positions[-1]:.0f}")
    print()

    return genotypes, positions, ts


def run_demo(
    num_particles=50,
    num_samples=20,
    sequence_length=100_000,
    population_size=10_000,
    mutation_rate=1e-8,
    recombination_rate=1e-8,
    random_seed=42
):
    """
    Run a complete demo of the particle filter.
    """
    print("\n")
    print("=" * 60)
    print("SMCSMC-MSPRIME DEMO")
    print("=" * 60)
    print()

    # Simulate data
    genotypes, positions, true_ts = simulate_data(
        num_samples=num_samples,
        sequence_length=sequence_length,
        population_size=population_size,
        mutation_rate=mutation_rate,
        recombination_rate=recombination_rate,
        random_seed=random_seed
    )

    # Create demographic model (matching the simulation)
    demography = smc.constant_Ne(Ne=population_size)

    print("=" * 60)
    print("RUNNING PARTICLE FILTER")
    print("=" * 60)
    print()

    # Run particle filter
    results = smc.run_particle_filter(
        genotypes=genotypes,
        positions=positions,
        num_particles=num_particles,
        demography=demography,
        mutation_rate=mutation_rate,
        recombination_rate=recombination_rate,
        sequence_length=sequence_length,
        random_seed=random_seed + 1,  # Different seed for particles
        verbose=True
    )

    # Display results
    print("\n")
    print("=" * 60)
    print("RESULTS")
    print("=" * 60)
    print(f"Log Likelihood: {results.log_likelihood:.2f}")
    print(f"Average log likelihood per site: {results.log_likelihood / len(positions):.4f}")
    print()
    print(f"Number of resampling events: {len(results.resample_positions)}")
    if len(results.resample_positions) > 0:
        print(f"Resample positions: {results.resample_positions[:5]}...")
    print()

    if len(results.ess_history) > 0:
        print(f"ESS statistics:")
        print(f"  Mean: {np.mean(results.ess_history):.2f}")
        print(f"  Min: {np.min(results.ess_history):.2f}")
        print(f"  Max: {np.max(results.ess_history):.2f}")
        print()

    # Final particle statistics
    print("Final particle ensemble:")
    weights = np.array([p.weight for p in results.final_particles])
    print(f"  Mean weight: {np.mean(weights):.6f}")
    print(f"  Std weight: {np.std(weights):.6f}")
    print(f"  Max weight: {np.max(weights):.6f}")
    print(f"  Min weight: {np.min(weights):.6f}")
    print()

    # Compare a few particles to true tree sequence
    print("Sample particle details:")
    for i in range(min(3, len(results.final_particles))):
        p = results.final_particles[i]
        print(f"  Particle {i}: {p.num_trees} trees, "
              f"weight={p.weight:.6f}, "
              f"{p.num_mutations} mutations")

    print()

    # PARAMETER INFERENCE
    print("\n")
    print("=" * 60)
    print("PARAMETER INFERENCE")
    print("=" * 60)
    print()
    print("Inferring population size from particle ensemble...")
    print()

    # Try different inference methods
    methods = ["tmrca", "pairwise", "events"]
    method_names = {
        "tmrca": "TMRCA-based",
        "pairwise": "Pairwise coalescence",
        "events": "All coalescence events"
    }

    for method in methods:
        print(f"Method: {method_names[method]}")
        inference_results = smc.infer_constant_Ne(
            particles=results.final_particles,
            method=method
        )

        comparison = smc.compare_to_true_parameters(
            inference_results,
            true_Ne=population_size
        )

        print(f"  Estimated Ne: {inference_results.Ne_estimate:,.0f} ± {inference_results.Ne_std:,.0f}")
        print(f"  True Ne: {population_size:,.0f}")
        print(f"  Relative error: {comparison['relative_error_percent']:.1f}%")

        if comparison['within_1_std']:
            print(f"  ✓ True value within 1σ")
        elif comparison['within_2_std']:
            print(f"  ✓ True value within 2σ")
        else:
            print(f"  ✗ True value outside 2σ")
        print()

    # Use TMRCA method for final summary
    print("Using TMRCA-based inference for final results:")
    final_inference = smc.infer_constant_Ne(
        particles=results.final_particles,
        method="tmrca"
    )

    print()
    print(smc.summarize_inference(final_inference, true_Ne=population_size))

    print()
    print("=" * 60)
    print("DEMO COMPLETE")
    print("=" * 60)
    print()

    return results, final_inference


def test_likelihood_calculation():
    """
    Test that the likelihood calculation works correctly.
    """
    print("\n")
    print("=" * 60)
    print("TESTING LIKELIHOOD CALCULATION")
    print("=" * 60)
    print()

    # Simulate small tree sequence
    ts = msprime.sim_ancestry(
        samples=10,
        population_size=10_000,
        sequence_length=10_000,
        recombination_rate=1e-8,
        random_seed=123,
        ploidy=1
    )
    ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=123)

    print(f"Simulated {ts.num_sites} sites")

    # Calculate likelihood using our method
    log_lik = smc.calculate_likelihood_at_variants(
        ts=ts,
        mutation_rate=1e-8,
        use_log=True
    )

    print(f"Log likelihood: {log_lik:.4f}")
    print(f"Likelihood: {np.exp(log_lik):.4e}")
    print()

    print("Likelihood test complete!")
    print()


if __name__ == "__main__":
    # Test likelihood calculation
    test_likelihood_calculation()

    # Run main demo
    results, inference = run_demo(
        num_particles=50000,
        num_samples=20,
        sequence_length=100_000,
        population_size=10_000,
        mutation_rate=1e-8,
        recombination_rate=1e-8,
        random_seed=42
    )

    print("Demo finished successfully!")
    print()
    print(f"Final Ne estimate: {inference.Ne_estimate:,.0f} (true: 10,000)")
