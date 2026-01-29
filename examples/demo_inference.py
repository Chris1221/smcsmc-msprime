#!/usr/bin/env python3
"""
Parameter inference demo with proper data scale.

This demonstrates that the particle filter can recover demographic parameters
with sufficient data (tens of thousands of sites) and particles (20k+).
"""

import numpy as np
import msprime
import smcsmc_msprime as smc
import time


def run_inference_demo(
    num_particles=20_000,
    num_samples=50,
    sequence_length=10_000_000,  # 10 Mb for many sites
    population_size=10_000,
    mutation_rate=2e-8,  # Higher rate for more sites
    recombination_rate=1e-8,
    random_seed=42
):
    """
    Run parameter inference demo with realistic data scale.
    """
    print("=" * 70)
    print("SMCSMC PARAMETER INFERENCE DEMO")
    print("=" * 70)
    print()
    print(f"Configuration:")
    print(f"  Samples: {num_samples}")
    print(f"  Sequence length: {sequence_length:,} bp")
    print(f"  True Ne: {population_size:,}")
    print(f"  Mutation rate: {mutation_rate:.2e}")
    print(f"  Recombination rate: {recombination_rate:.2e}")
    print(f"  Num particles: {num_particles:,}")
    print()

    # Simulate data
    print("=" * 70)
    print("STEP 1: SIMULATING DATA")
    print("=" * 70)
    print()

    start_time = time.time()

    print("Simulating ancestry...")
    ts = msprime.sim_ancestry(
        samples=num_samples,
        population_size=population_size,
        sequence_length=sequence_length,
        recombination_rate=recombination_rate,
        random_seed=random_seed,
        ploidy=1
    )

    print(f"  {ts.num_trees} local trees")
    print(f"  {ts.num_nodes} nodes, {ts.num_edges} edges")
    print()

    print("Adding mutations...")
    ts = msprime.sim_mutations(ts, rate=mutation_rate, random_seed=random_seed)

    print(f"  {ts.num_sites:,} variant sites")
    print(f"  {ts.num_mutations:,} mutations")
    print()

    genotypes = ts.genotype_matrix()
    positions = np.array([site.position for site in ts.sites()])

    sim_time = time.time() - start_time
    print(f"Simulation time: {sim_time:.1f} seconds")
    print()

    # Run particle filter
    print("=" * 70)
    print("STEP 2: RUNNING PARTICLE FILTER")
    print("=" * 70)
    print()

    demography = smc.constant_Ne(Ne=population_size)

    start_time = time.time()

    results = smc.run_particle_filter(
        genotypes=genotypes,
        positions=positions,
        num_particles=num_particles,
        demography=demography,
        mutation_rate=mutation_rate,
        recombination_rate=recombination_rate,
        sequence_length=sequence_length,
        random_seed=random_seed + 1,
        verbose=True
    )

    filter_time = time.time() - start_time
    print()
    print(f"Particle filter time: {filter_time:.1f} seconds")
    print(f"Time per site: {filter_time / len(positions) * 1000:.2f} ms")
    print()

    # Summary statistics
    print("=" * 70)
    print("STEP 3: PARTICLE FILTER SUMMARY")
    print("=" * 70)
    print()
    print(f"Log Likelihood: {results.log_likelihood:.2f}")
    print(f"Log likelihood per site: {results.log_likelihood / len(positions):.4f}")
    print()
    print(f"Resampling events: {len(results.resample_positions)}")
    print(f"Resampling rate: {len(results.resample_positions) / len(positions) * 100:.1f}%")
    print()
    print(f"ESS statistics:")
    print(f"  Mean: {np.mean(results.ess_history):,.1f}")
    print(f"  Min: {np.min(results.ess_history):,.1f}")
    print(f"  Max: {np.max(results.ess_history):,.1f}")
    print(f"  Median: {np.median(results.ess_history):,.1f}")
    print()

    # Parameter inference
    print("=" * 70)
    print("STEP 4: PARAMETER INFERENCE")
    print("=" * 70)
    print()

    print("Testing different inference methods...")
    print()

    methods = ["tmrca", "pairwise", "events"]
    method_names = {
        "tmrca": "TMRCA-based (recommended)",
        "pairwise": "Pairwise coalescence",
        "events": "All coalescence events"
    }

    results_by_method = {}

    for method in methods:
        print(f"Method: {method_names[method]}")
        start_time = time.time()

        inference = smc.infer_constant_Ne(
            particles=results.final_particles,
            method=method
        )

        inference_time = time.time() - start_time

        comparison = smc.compare_to_true_parameters(inference, population_size)

        print(f"  Estimated Ne: {inference.Ne_estimate:,.0f} ± {inference.Ne_std:,.0f}")
        print(f"  True Ne: {population_size:,}")
        print(f"  Error: {comparison['relative_error_percent']:.1f}%")
        print(f"  Inference time: {inference_time:.2f}s")

        if comparison['within_1_std']:
            print(f"  ✓ Excellent: True value within 1σ")
        elif comparison['within_2_std']:
            print(f"  ✓ Good: True value within 2σ")
        else:
            print(f"  ⚠ Poor: True value outside 2σ")

        results_by_method[method] = inference
        print()

    # Detailed summary for best method
    print("=" * 70)
    print("DETAILED RESULTS (TMRCA method)")
    print("=" * 70)
    print()

    best_inference = results_by_method["tmrca"]
    print(smc.summarize_inference(best_inference, true_Ne=population_size))
    print()

    # Final summary
    print("=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    print()
    print(f"True Ne:      {population_size:,}")
    print(f"Estimated Ne: {best_inference.Ne_estimate:,.0f} ± {best_inference.Ne_std:,.0f}")
    print()

    error_pct = abs((best_inference.Ne_estimate - population_size) / population_size * 100)
    print(f"Accuracy: {100 - error_pct:.1f}% (error: {error_pct:.1f}%)")
    print()

    comparison = smc.compare_to_true_parameters(best_inference, population_size)
    if comparison['within_1_std']:
        print("✓ Inference SUCCESSFUL: True parameter recovered within 1σ")
    elif comparison['within_2_std']:
        print("✓ Inference GOOD: True parameter recovered within 2σ")
    elif error_pct < 20:
        print("⚠ Inference ACCEPTABLE: Relative error < 20%")
    else:
        print("✗ Inference POOR: Large deviation from true value")

    print()
    print("=" * 70)
    print()

    return results, best_inference


if __name__ == "__main__":
    import sys

    # Default parameters
    params = {
        "num_particles": 20_000,
        "num_samples": 50,
        "sequence_length": 10_000_000,
        "population_size": 10_000,
        "mutation_rate": 2e-8,
        "recombination_rate": 1e-8,
        "random_seed": 42
    }

    # Parse command line arguments for quick adjustments
    if len(sys.argv) > 1:
        params["num_particles"] = int(sys.argv[1])
    if len(sys.argv) > 2:
        params["num_samples"] = int(sys.argv[2])
    if len(sys.argv) > 3:
        params["population_size"] = int(sys.argv[3])

    print()
    print("To customize: python demo_inference.py [num_particles] [num_samples] [true_Ne]")
    print(f"Using: {params['num_particles']:,} particles, {params['num_samples']} samples, Ne={params['population_size']:,}")
    print()

    results, inference = run_inference_demo(**params)

    print("Demo complete!")
