"""
Resampling methods for particle filters.

This module implements systematic resampling and related utilities for
maintaining particle diversity in the SMC particle filter.
"""

import numpy as np
from typing import List, Optional, Tuple
from .particle import Particle


def effective_sample_size(weights: np.ndarray) -> float:
    """
    Calculate the effective sample size (ESS).

    ESS = 1 / sum(normalized_weights^2)

    A low ESS indicates that most particles have negligible weight,
    suggesting that resampling is needed.

    Args:
        weights: Array of particle weights (not necessarily normalized)

    Returns:
        Effective sample size (between 1 and N)
    """
    # Normalize weights
    weights = np.array(weights)
    if np.sum(weights) == 0:
        return 0.0

    normalized_weights = weights / np.sum(weights)

    # Calculate ESS
    ess = 1.0 / np.sum(normalized_weights ** 2)

    return ess


def normalize_weights(weights: np.ndarray) -> np.ndarray:
    """
    Normalize weights to sum to 1.

    Args:
        weights: Array of weights

    Returns:
        Normalized weights
    """
    weights = np.array(weights)
    total = np.sum(weights)

    if total == 0:
        # If all weights are zero, return uniform weights
        return np.ones_like(weights) / len(weights)

    return weights / total


def systematic_resampling(
    particles: List[Particle],
    weights: Optional[np.ndarray] = None
) -> List[Particle]:
    """
    Perform systematic resampling of particles.

    Systematic resampling is a low-variance resampling method that
    uses a single random number and deterministic spacing to select
    particles. This reduces sampling variability compared to multinomial
    resampling.

    Args:
        particles: List of Particle objects
        weights: Optional array of weights (default: use particle.weight)

    Returns:
        List of resampled Particle objects (same length as input)
    """
    N = len(particles)

    # Extract weights if not provided
    if weights is None:
        weights = np.array([p.weight for p in particles])

    # Normalize weights
    weights = normalize_weights(weights)

    # Generate systematic sample positions
    # Single random offset, then deterministic spacing
    positions = (np.arange(N) + np.random.uniform(0, 1)) / N

    # Calculate cumulative sum of weights
    cumsum = np.cumsum(weights)

    # Resample particles
    new_particles = []
    j = 0

    for i, pos in enumerate(positions):
        # Find the particle corresponding to this position
        while j < N and cumsum[j] < pos:
            j += 1

        if j >= N:
            j = N - 1

        # Copy the selected particle
        new_particle = particles[j].copy()

        # Reset weights to uniform (1/N)
        new_particle.weight = 1.0 / N
        new_particle.pilot_weight = 1.0
        new_particle.multiplicity = 1

        new_particles.append(new_particle)

    return new_particles


def multinomial_resampling(
    particles: List[Particle],
    weights: Optional[np.ndarray] = None
) -> List[Particle]:
    """
    Perform multinomial resampling of particles.

    This is a simpler but higher-variance alternative to systematic resampling.
    Each new particle is independently drawn from the weighted distribution.

    Args:
        particles: List of Particle objects
        weights: Optional array of weights (default: use particle.weight)

    Returns:
        List of resampled Particle objects
    """
    N = len(particles)

    # Extract weights if not provided
    if weights is None:
        weights = np.array([p.weight for p in particles])

    # Normalize weights
    weights = normalize_weights(weights)

    # Sample indices
    indices = np.random.choice(N, size=N, replace=True, p=weights)

    # Create new particle list
    new_particles = []
    for idx in indices:
        new_particle = particles[idx].copy()
        new_particle.weight = 1.0 / N
        new_particle.pilot_weight = 1.0
        new_particle.multiplicity = 1
        new_particles.append(new_particle)

    return new_particles


def stratified_resampling(
    particles: List[Particle],
    weights: Optional[np.ndarray] = None
) -> List[Particle]:
    """
    Perform stratified resampling of particles.

    Stratified resampling divides the [0,1] interval into N equal strata
    and samples one point uniformly from each stratum. This has similar
    low-variance properties to systematic resampling.

    Args:
        particles: List of Particle objects
        weights: Optional array of weights

    Returns:
        List of resampled Particle objects
    """
    N = len(particles)

    # Extract weights if not provided
    if weights is None:
        weights = np.array([p.weight for p in particles])

    # Normalize weights
    weights = normalize_weights(weights)

    # Generate stratified sample positions
    positions = (np.arange(N) + np.random.uniform(0, 1, size=N)) / N

    # Calculate cumulative sum
    cumsum = np.cumsum(weights)

    # Resample particles
    new_particles = []
    j = 0

    for pos in positions:
        while j < N and cumsum[j] < pos:
            j += 1

        if j >= N:
            j = N - 1

        new_particle = particles[j].copy()
        new_particle.weight = 1.0 / N
        new_particle.pilot_weight = 1.0
        new_particle.multiplicity = 1
        new_particles.append(new_particle)

    return new_particles


def residual_resampling(
    particles: List[Particle],
    weights: Optional[np.ndarray] = None
) -> List[Particle]:
    """
    Perform residual resampling of particles.

    Residual resampling first deterministically replicates particles
    according to the integer part of N*weight, then uses another
    resampling scheme for the fractional parts.

    Args:
        particles: List of Particle objects
        weights: Optional array of weights

    Returns:
        List of resampled Particle objects
    """
    N = len(particles)

    # Extract weights if not provided
    if weights is None:
        weights = np.array([p.weight for p in particles])

    # Normalize weights
    weights = normalize_weights(weights)

    # Calculate number of deterministic copies
    N_weights = N * weights
    N_deterministic = np.floor(N_weights).astype(int)

    # Deterministic resampling
    new_particles = []
    for i, count in enumerate(N_deterministic):
        for _ in range(count):
            new_particles.append(particles[i].copy())

    # Residual weights for stochastic resampling
    residual = N_weights - N_deterministic
    residual_sum = np.sum(residual)

    # Stochastic resampling for remaining slots
    remaining = N - len(new_particles)
    if remaining > 0 and residual_sum > 0:
        residual_weights = residual / residual_sum
        indices = np.random.choice(N, size=remaining, replace=True, p=residual_weights)
        for idx in indices:
            new_particles.append(particles[idx].copy())

    # Reset weights
    for p in new_particles:
        p.weight = 1.0 / N
        p.pilot_weight = 1.0
        p.multiplicity = 1

    return new_particles


def should_resample(
    weights: np.ndarray,
    threshold: float = 0.5,
    ess_threshold: Optional[float] = None
) -> bool:
    """
    Determine if resampling should be performed.

    Resampling is triggered when the effective sample size drops below
    a threshold, indicating particle degeneracy.

    Args:
        weights: Array of particle weights
        threshold: ESS threshold as fraction of N (default: 0.5)
        ess_threshold: Absolute ESS threshold (overrides threshold if provided)

    Returns:
        True if resampling should be performed
    """
    N = len(weights)
    ess = effective_sample_size(weights)

    if ess_threshold is not None:
        return ess < ess_threshold
    else:
        return ess < threshold * N


def compute_resampling_stats(
    particles: List[Particle],
    weights: Optional[np.ndarray] = None
) -> dict:
    """
    Compute statistics about the particle ensemble.

    Args:
        particles: List of Particle objects
        weights: Optional array of weights

    Returns:
        Dictionary with statistics: ess, max_weight, min_weight, etc.
    """
    if weights is None:
        weights = np.array([p.weight for p in particles])

    normalized_weights = normalize_weights(weights)
    ess = effective_sample_size(weights)

    return {
        'ess': ess,
        'ess_ratio': ess / len(particles),
        'max_weight': np.max(normalized_weights),
        'min_weight': np.min(normalized_weights),
        'mean_weight': np.mean(normalized_weights),
        'std_weight': np.std(normalized_weights),
        'cv_weight': np.std(normalized_weights) / np.mean(normalized_weights)
        if np.mean(normalized_weights) > 0 else 0.0,
        'num_effective_particles': int(ess),
        'num_total_particles': len(particles)
    }
