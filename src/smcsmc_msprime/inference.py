"""
Parameter inference from particle filter results.

This module implements demographic parameter estimation from the particle ensemble
using coalescence time distributions and event counting.
"""

import numpy as np
import tskit
from typing import List, Dict, Tuple, Optional
from collections import defaultdict
from dataclasses import dataclass, field

from .particle import Particle


@dataclass
class InferenceResults:
    """Results from parameter inference."""

    Ne_estimate: float = 0.0
    """Estimated effective population size"""

    Ne_std: float = 0.0
    """Standard deviation of Ne estimate"""

    coalescence_times: List[float] = field(default_factory=list)
    """Coalescence times extracted from particles"""

    weighted_coalescence_times: List[Tuple[float, float]] = field(default_factory=list)
    """(time, weight) pairs for weighted analysis"""

    num_coalescence_events: int = 0
    """Total number of coalescence events"""

    mean_coalescence_time: float = 0.0
    """Mean coalescence time"""

    tmrca: float = 0.0
    """Time to most recent common ancestor (TMRCA)"""


def extract_coalescence_events(
    particles: List[Particle],
    weights: Optional[np.ndarray] = None
) -> List[Tuple[float, float]]:
    """
    Extract coalescence events from particle ensemble.

    Returns (time, weight) pairs for each coalescence event across all particles.

    Args:
        particles: List of Particle objects
        weights: Optional weights (default: use particle.weight)

    Returns:
        List of (coalescence_time, weight) tuples
    """
    if weights is None:
        weights = np.array([p.weight for p in particles])

    # Normalize weights
    weights = weights / np.sum(weights)

    events = []

    for particle, weight in zip(particles, weights):
        ts = particle.ts

        # Extract coalescence times from the tree sequence
        # Each internal node represents a coalescence event
        for node in ts.nodes():
            # Skip sample nodes (leaves)
            if node.id in ts.samples():
                continue

            # This is an internal node = coalescence event
            coalescence_time = node.time

            # Weight this event by the particle weight
            events.append((coalescence_time, weight))

    return events


def extract_coalescence_times_per_tree(
    particles: List[Particle],
    weights: Optional[np.ndarray] = None
) -> List[Tuple[float, float, float]]:
    """
    Extract coalescence times from each tree in the tree sequence.

    Returns (position, time, weight) for coalescences.

    Args:
        particles: List of Particle objects
        weights: Optional weights

    Returns:
        List of (position, coalescence_time, weight) tuples
    """
    if weights is None:
        weights = np.array([p.weight for p in particles])

    weights = weights / np.sum(weights)

    events = []

    for particle, weight in zip(particles, weights):
        # For each tree in the sequence
        for tree in particle.ts.trees():
            # Get the root time (TMRCA for this tree)
            root_time = tree.time(tree.root)
            position = tree.interval.left

            events.append((position, root_time, weight))

    return events


def estimate_Ne_from_coalescence_times(
    coalescence_events: List[Tuple[float, float]],
    num_samples: int,
    method: str = "mean"
) -> Tuple[float, float]:
    """
    Estimate Ne from coalescence times.

    Under the coalescent, the expected time to the first coalescence
    in a sample of n haploids is 2*Ne / (n choose 2) = 4*Ne / (n*(n-1))

    Args:
        coalescence_events: List of (time, weight) tuples
        num_samples: Number of samples
        method: Estimation method ("mean", "median", "weighted_mean")

    Returns:
        (Ne_estimate, Ne_std) tuple
    """
    if len(coalescence_events) == 0:
        return 0.0, 0.0

    times, weights = zip(*coalescence_events)
    times = np.array(times)
    weights = np.array(weights)

    # Normalize weights
    weights = weights / np.sum(weights)

    if method == "mean":
        mean_time = np.mean(times)
    elif method == "median":
        mean_time = np.median(times)
    elif method == "weighted_mean":
        mean_time = np.sum(times * weights)
    else:
        raise ValueError(f"Unknown method: {method}")

    # Under coalescent with n samples:
    # E[T_2] = 4*Ne / (n*(n-1))  for diploids
    # E[T_2] = 2*Ne / (n*(n-1))  for haploids

    # So: Ne = E[T] * n * (n-1) / 2  for haploids
    n = num_samples
    Ne_estimate = mean_time * n * (n - 1) / 2.0

    # Standard deviation (rough approximation)
    std_time = np.sqrt(np.sum(weights * (times - mean_time)**2))
    Ne_std = std_time * n * (n - 1) / 2.0

    return Ne_estimate, Ne_std


def estimate_Ne_from_tmrca(
    particles: List[Particle],
    weights: Optional[np.ndarray] = None
) -> Tuple[float, float]:
    """
    Estimate Ne from TMRCA (time to most recent common ancestor).

    This uses the mean TMRCA across all trees and particles.

    Args:
        particles: List of particles
        weights: Optional weights

    Returns:
        (Ne_estimate, Ne_std) tuple
    """
    if weights is None:
        weights = np.array([p.weight for p in particles])

    weights = weights / np.sum(weights)

    tmrca_values = []
    tmrca_weights = []

    for particle, weight in zip(particles, weights):
        # Get TMRCA for each tree
        for tree in particle.ts.trees():
            tmrca = tree.time(tree.root)
            tmrca_values.append(tmrca)
            tmrca_weights.append(weight)

    tmrca_values = np.array(tmrca_values)
    tmrca_weights = np.array(tmrca_weights)
    tmrca_weights = tmrca_weights / np.sum(tmrca_weights)

    # Weighted mean TMRCA
    mean_tmrca = np.sum(tmrca_values * tmrca_weights)
    std_tmrca = np.sqrt(np.sum(tmrca_weights * (tmrca_values - mean_tmrca)**2))

    # For haploid samples in msprime: E[TMRCA] = 2*Ne*(1 - 1/n)
    # Therefore: Ne = TMRCA / (2*(1 - 1/n))

    n = particles[0].num_samples

    Ne_estimate = mean_tmrca / (2.0 * (1.0 - 1.0/n))
    Ne_std = std_tmrca / (2.0 * (1.0 - 1.0/n))

    return Ne_estimate, Ne_std


def estimate_Ne_from_pairwise_coalescence(
    particles: List[Particle],
    weights: Optional[np.ndarray] = None,
    num_pairs: int = 100
) -> Tuple[float, float]:
    """
    Estimate Ne from pairwise coalescence times.

    For each pair of samples, find their coalescence time.
    Expected pairwise coalescence time = 2*Ne for haploids.

    Args:
        particles: List of particles
        weights: Optional weights
        num_pairs: Number of sample pairs to use (for efficiency)

    Returns:
        (Ne_estimate, Ne_std) tuple
    """
    if weights is None:
        weights = np.array([p.weight for p in particles])

    weights = weights / np.sum(weights)

    pairwise_times = []
    pairwise_weights = []

    for particle, weight in zip(particles, weights):
        ts = particle.ts
        samples = list(ts.samples())

        # Sample random pairs
        n = len(samples)
        if n < 2:
            continue

        # Use the first tree for simplicity (or could average across trees)
        tree = ts.first()

        # Sample random pairs
        num_sample_pairs = min(num_pairs, n * (n - 1) // 2)

        for _ in range(num_sample_pairs):
            i, j = np.random.choice(n, size=2, replace=False)
            sample_i = samples[i]
            sample_j = samples[j]

            # Find MRCA of this pair
            mrca = tree.mrca(sample_i, sample_j)
            if mrca != tskit.NULL:
                tmrca_pair = tree.time(mrca)
                pairwise_times.append(tmrca_pair)
                pairwise_weights.append(weight)

    if len(pairwise_times) == 0:
        return 0.0, 0.0

    pairwise_times = np.array(pairwise_times)
    pairwise_weights = np.array(pairwise_weights)
    pairwise_weights = pairwise_weights / np.sum(pairwise_weights)

    # Weighted mean
    mean_pairwise_time = np.sum(pairwise_times * pairwise_weights)
    std_pairwise_time = np.sqrt(
        np.sum(pairwise_weights * (pairwise_times - mean_pairwise_time)**2)
    )

    # For haploids: E[pairwise coalescence time] = 2*Ne
    Ne_estimate = mean_pairwise_time / 2.0
    Ne_std = std_pairwise_time / 2.0

    return Ne_estimate, Ne_std


def infer_constant_Ne(
    particles: List[Particle],
    weights: Optional[np.ndarray] = None,
    method: str = "tmrca"
) -> InferenceResults:
    """
    Infer constant population size Ne from particle ensemble.

    Args:
        particles: List of particles
        weights: Optional weights (default: use particle.weight)
        method: Estimation method:
            - "tmrca": Use TMRCA distribution
            - "pairwise": Use pairwise coalescence times
            - "events": Use all coalescence events

    Returns:
        InferenceResults object with Ne estimate
    """
    if weights is None:
        weights = np.array([p.weight for p in particles])

    results = InferenceResults()

    if method == "tmrca":
        Ne_estimate, Ne_std = estimate_Ne_from_tmrca(particles, weights)
    elif method == "pairwise":
        Ne_estimate, Ne_std = estimate_Ne_from_pairwise_coalescence(
            particles, weights
        )
    elif method == "events":
        # Extract all coalescence events
        events = extract_coalescence_events(particles, weights)
        results.weighted_coalescence_times = events
        results.num_coalescence_events = len(events)

        if len(events) > 0:
            times, event_weights = zip(*events)
            results.coalescence_times = list(times)
            results.mean_coalescence_time = np.mean(times)

        Ne_estimate, Ne_std = estimate_Ne_from_coalescence_times(
            events, particles[0].num_samples, method="weighted_mean"
        )
    else:
        raise ValueError(f"Unknown method: {method}")

    results.Ne_estimate = Ne_estimate
    results.Ne_std = Ne_std

    # Also get TMRCA for reporting
    tmrca_events = extract_coalescence_times_per_tree(particles, weights)
    if len(tmrca_events) > 0:
        _, times, _ = zip(*tmrca_events)
        results.tmrca = np.mean(times)

    return results


def compare_to_true_parameters(
    inference_results: InferenceResults,
    true_Ne: float
) -> Dict[str, float]:
    """
    Compare inferred parameters to true values.

    Args:
        inference_results: Results from inference
        true_Ne: True effective population size

    Returns:
        Dictionary with comparison statistics
    """
    estimated_Ne = inference_results.Ne_estimate

    # Relative error
    rel_error = (estimated_Ne - true_Ne) / true_Ne
    abs_rel_error = abs(rel_error)

    # Standard errors
    if inference_results.Ne_std > 0:
        z_score = (estimated_Ne - true_Ne) / inference_results.Ne_std
    else:
        z_score = 0.0

    return {
        "true_Ne": true_Ne,
        "estimated_Ne": estimated_Ne,
        "Ne_std": inference_results.Ne_std,
        "absolute_error": abs(estimated_Ne - true_Ne),
        "relative_error": rel_error,
        "absolute_relative_error": abs_rel_error,
        "relative_error_percent": rel_error * 100,
        "z_score": z_score,
        "within_1_std": abs(z_score) <= 1.0,
        "within_2_std": abs(z_score) <= 2.0,
    }


def summarize_inference(
    inference_results: InferenceResults,
    true_Ne: Optional[float] = None
) -> str:
    """
    Create a summary string of inference results.

    Args:
        inference_results: Results from inference
        true_Ne: Optional true Ne for comparison

    Returns:
        Formatted summary string
    """
    lines = []
    lines.append("=" * 60)
    lines.append("PARAMETER INFERENCE RESULTS")
    lines.append("=" * 60)
    lines.append(f"Estimated Ne: {inference_results.Ne_estimate:,.0f}")
    lines.append(f"Std deviation: {inference_results.Ne_std:,.0f}")
    lines.append(f"95% CI: [{inference_results.Ne_estimate - 2*inference_results.Ne_std:,.0f}, "
                f"{inference_results.Ne_estimate + 2*inference_results.Ne_std:,.0f}]")
    lines.append("")

    if true_Ne is not None:
        comparison = compare_to_true_parameters(inference_results, true_Ne)
        lines.append(f"True Ne: {true_Ne:,.0f}")
        lines.append(f"Absolute error: {comparison['absolute_error']:,.0f}")
        lines.append(f"Relative error: {comparison['relative_error_percent']:.1f}%")
        lines.append(f"Z-score: {comparison['z_score']:.2f}")

        if comparison['within_1_std']:
            lines.append("✓ True value within 1 std deviation")
        elif comparison['within_2_std']:
            lines.append("✓ True value within 2 std deviations")
        else:
            lines.append("✗ True value outside 2 std deviations")
        lines.append("")

    if inference_results.num_coalescence_events > 0:
        lines.append(f"Coalescence events analyzed: {inference_results.num_coalescence_events}")
        lines.append(f"Mean coalescence time: {inference_results.mean_coalescence_time:.1f} generations")

    if inference_results.tmrca > 0:
        lines.append(f"Mean TMRCA: {inference_results.tmrca:.1f} generations")

    lines.append("=" * 60)

    return "\n".join(lines)
