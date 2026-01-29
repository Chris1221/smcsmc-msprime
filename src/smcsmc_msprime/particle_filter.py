"""
Sequential Monte Carlo (SMC) particle filter for demographic inference.

This module implements the main particle filter algorithm that propagates
particles along the genome, updating weights at variant sites and resampling
when particle degeneracy occurs.
"""

import numpy as np
import msprime
import tskit
from typing import Optional, List, Dict, Tuple
from dataclasses import dataclass, field

from .particle import Particle, create_initial_particles
from .likelihood import calculate_site_likelihood
from .resampling import (
    systematic_resampling,
    effective_sample_size,
    should_resample,
    compute_resampling_stats
)


@dataclass
class ParticleFilterConfig:
    """Configuration for the particle filter."""

    num_particles: int = 100
    """Number of particles in the ensemble"""

    resample_threshold: float = 0.5
    """Resample when ESS < threshold * num_particles"""

    mutation_rate: float = 1e-8
    """Mutation rate per base per generation"""

    recombination_rate: float = 1e-8
    """Recombination rate per base per generation"""

    sequence_length: float = 1e6
    """Length of the genomic sequence"""

    random_seed: Optional[int] = None
    """Random seed for reproducibility"""

    verbose: bool = True
    """Print progress information"""


@dataclass
class ParticleFilterResults:
    """Results from running the particle filter."""

    log_likelihood: float = 0.0
    """Total log likelihood along the sequence"""

    normalization_constants: List[float] = field(default_factory=list)
    """Normalization constants at each variant site"""

    ess_history: List[float] = field(default_factory=list)
    """Effective sample size at each step"""

    resample_positions: List[float] = field(default_factory=list)
    """Positions where resampling occurred"""

    final_particles: List[Particle] = field(default_factory=list)
    """Particles at the end of the filter"""

    positions: List[float] = field(default_factory=list)
    """Positions of variant sites processed"""


class ParticleFilter:
    """
    Sequential Monte Carlo particle filter for genealogical inference.

    This implements the SMCSMC algorithm using msprime/tskit for efficient
    representation of genealogies.
    """

    def __init__(
        self,
        num_samples: int,
        genotypes: np.ndarray,
        positions: np.ndarray,
        demography: Optional[msprime.Demography] = None,
        config: Optional[ParticleFilterConfig] = None
    ):
        """
        Initialize the particle filter.

        Args:
            num_samples: Number of haploid samples
            genotypes: Genotype matrix (num_sites x num_samples)
            positions: Positions of variant sites
            demography: Demographic model (default: constant Ne=10000)
            config: Configuration object
        """
        self.num_samples = num_samples
        self.genotypes = genotypes
        self.positions = positions
        self.demography = demography
        self.config = config or ParticleFilterConfig()

        # Use default constant demography if none provided
        if self.demography is None:
            from .demography import constant_Ne
            self.demography = constant_Ne(Ne=10000)

        # Initialize particles
        self.particles: List[Particle] = []
        self.results = ParticleFilterResults()

        # Validate inputs
        self._validate_inputs()

    def _validate_inputs(self):
        """Validate input data."""
        num_sites, num_samples = self.genotypes.shape

        if num_samples != self.num_samples:
            raise ValueError(
                f"Genotypes have {num_samples} samples but expected {self.num_samples}"
            )

        if len(self.positions) != num_sites:
            raise ValueError(
                f"Have {len(self.positions)} positions but {num_sites} sites"
            )

        # Check positions are sorted
        if not np.all(self.positions[:-1] <= self.positions[1:]):
            raise ValueError("Positions must be sorted")

    def initialize_particles(self):
        """
        Initialize the particle ensemble.

        Each particle is created by simulating from the prior demographic model.
        """
        if self.config.verbose:
            print(f"Initializing {self.config.num_particles} particles...")

        self.particles = create_initial_particles(
            num_particles=self.config.num_particles,
            sample_size=self.num_samples,
            sequence_length=self.config.sequence_length,
            demography=self.demography,
            recombination_rate=self.config.recombination_rate,
            random_seed=self.config.random_seed
        )

        if self.config.verbose:
            print(f"Initialized {len(self.particles)} particles")

    def update_weights_at_site(self, site_idx: int):
        """
        Update particle weights at a variant site using likelihoods.

        Args:
            site_idx: Index of the variant site
        """
        position = self.positions[site_idx]
        site_genotypes = self.genotypes[site_idx, :]

        total_weight = 0.0

        for particle in self.particles:
            # Get the tree at this position
            tree = particle.ts.at(position)

            # Calculate likelihood
            likelihood = calculate_site_likelihood(
                tree=tree,
                genotypes=site_genotypes,
                mutation_rate=self.config.mutation_rate,
                use_log=False
            )

            # Update particle weight
            particle.update_weight(likelihood)
            total_weight += particle.weight

        # Store normalization constant
        if total_weight > 0:
            self.results.normalization_constants.append(total_weight)

            # Normalize weights
            for particle in self.particles:
                particle.normalize_weight(total_weight)

            # Update log likelihood
            self.results.log_likelihood += np.log(total_weight)
        else:
            # All particles have zero weight - this is bad!
            print(f"WARNING: All particles have zero weight at position {position}")
            self.results.normalization_constants.append(0.0)

    def resample_if_needed(self, position: float):
        """
        Check if resampling is needed and perform it if necessary.

        Args:
            position: Current genomic position (for logging)
        """
        weights = np.array([p.weight for p in self.particles])
        ess = effective_sample_size(weights)

        # Store ESS
        self.results.ess_history.append(ess)

        # Check if resampling is needed
        if should_resample(weights, threshold=self.config.resample_threshold):
            if self.config.verbose:
                print(f"  Resampling at position {position:.0f} (ESS={ess:.1f})")

            # Perform systematic resampling
            self.particles = systematic_resampling(self.particles)

            # Record resample position
            self.results.resample_positions.append(position)

    def run(self) -> ParticleFilterResults:
        """
        Run the particle filter along the sequence.

        Returns:
            ParticleFilterResults object with filtering results
        """
        if self.config.verbose:
            print("=" * 60)
            print("Running SMCSMC Particle Filter")
            print("=" * 60)
            print(f"Num particles: {self.config.num_particles}")
            print(f"Num samples: {self.num_samples}")
            print(f"Num sites: {len(self.positions)}")
            print(f"Sequence length: {self.config.sequence_length:.0f}")
            print(f"Mutation rate: {self.config.mutation_rate:.2e}")
            print(f"Recombination rate: {self.config.recombination_rate:.2e}")
            print("=" * 60)

        # Initialize particles
        self.initialize_particles()

        # Process each variant site
        for site_idx in range(len(self.positions)):
            position = self.positions[site_idx]

            if self.config.verbose and site_idx % 100 == 0:
                print(f"Processing site {site_idx}/{len(self.positions)} "
                      f"at position {position:.0f}")

            # Store position
            self.results.positions.append(position)

            # Update weights based on data likelihood
            self.update_weights_at_site(site_idx)

            # Check if resampling is needed
            self.resample_if_needed(position)

        # Store final particles
        self.results.final_particles = self.particles

        if self.config.verbose:
            print("=" * 60)
            print("Particle Filter Complete")
            print(f"Log Likelihood: {self.results.log_likelihood:.2f}")
            print(f"Num resamplings: {len(self.results.resample_positions)}")
            if len(self.results.ess_history) > 0:
                print(f"Mean ESS: {np.mean(self.results.ess_history):.1f}")
            print("=" * 60)

        return self.results

    def get_weighted_mean(self, values: np.ndarray) -> float:
        """
        Calculate weighted mean across particles.

        Args:
            values: Array of values (one per particle)

        Returns:
            Weighted mean
        """
        weights = np.array([p.weight for p in self.particles])
        weights = weights / np.sum(weights)
        return np.sum(values * weights)

    def get_statistics(self) -> Dict:
        """
        Get statistics about the particle ensemble.

        Returns:
            Dictionary of statistics
        """
        stats = compute_resampling_stats(self.particles)
        stats['log_likelihood'] = self.results.log_likelihood
        stats['num_resamplings'] = len(self.results.resample_positions)
        return stats


def run_particle_filter(
    genotypes: np.ndarray,
    positions: np.ndarray,
    num_particles: int = 100,
    demography: Optional[msprime.Demography] = None,
    mutation_rate: float = 1e-8,
    recombination_rate: float = 1e-8,
    sequence_length: Optional[float] = None,
    random_seed: Optional[int] = None,
    verbose: bool = True
) -> ParticleFilterResults:
    """
    Convenience function to run the particle filter.

    Args:
        genotypes: Genotype matrix (num_sites x num_samples)
        positions: Positions of variant sites
        num_particles: Number of particles
        demography: Demographic model
        mutation_rate: Mutation rate per base per generation
        recombination_rate: Recombination rate per base per generation
        sequence_length: Length of sequence (default: max position)
        random_seed: Random seed
        verbose: Print progress

    Returns:
        ParticleFilterResults
    """
    num_sites, num_samples = genotypes.shape

    if sequence_length is None:
        sequence_length = positions[-1] * 1.1  # Add 10% buffer

    # Create config
    config = ParticleFilterConfig(
        num_particles=num_particles,
        mutation_rate=mutation_rate,
        recombination_rate=recombination_rate,
        sequence_length=sequence_length,
        random_seed=random_seed,
        verbose=verbose
    )

    # Create and run filter
    pf = ParticleFilter(
        num_samples=num_samples,
        genotypes=genotypes,
        positions=positions,
        demography=demography,
        config=config
    )

    return pf.run()
