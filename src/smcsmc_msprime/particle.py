"""
Particle class for SMCSMC particle filter.

A Particle wraps a tskit.TreeSequence (representing an ARG) along with
importance sampling weights for the particle filter.
"""

import copy as copy_module
import tskit
import numpy as np
from typing import Optional


class Particle:
    """
    A particle in the SMC particle filter.

    Each particle represents a possible genealogical history (ARG) for the sample,
    along with importance sampling weights.

    Attributes:
        ts: The tree sequence representing the genealogy/ARG
        weight: Posterior importance weight
        pilot_weight: Pilot weight for importance sampling
        multiplicity: Number of copies of this particle (for resampling)
        log_likelihood: Cached log likelihood value
    """

    def __init__(
        self,
        ts: tskit.TreeSequence,
        weight: float = 1.0,
        pilot_weight: float = 1.0,
        multiplicity: int = 1
    ):
        """
        Initialize a particle.

        Args:
            ts: Tree sequence representing the genealogy
            weight: Initial posterior weight (default 1.0)
            pilot_weight: Initial pilot weight (default 1.0)
            multiplicity: Initial multiplicity (default 1)
        """
        self.ts = ts
        self.weight = weight
        self.pilot_weight = pilot_weight
        self.multiplicity = multiplicity
        self.log_likelihood = 0.0

    def copy(self) -> "Particle":
        """
        Create a copy of this particle.

        Tree sequences are immutable, so we can share the reference.
        We only copy the weights and metadata.

        Returns:
            A new Particle with the same tree sequence but independent weights
        """
        return Particle(
            ts=self.ts,  # TreeSequences are immutable, no need to copy
            weight=self.weight,
            pilot_weight=self.pilot_weight,
            multiplicity=self.multiplicity
        )

    def __repr__(self) -> str:
        """String representation of the particle."""
        return (
            f"Particle(num_trees={self.ts.num_trees}, "
            f"num_samples={self.ts.num_samples}, "
            f"weight={self.weight:.6f}, "
            f"pilot_weight={self.pilot_weight:.6f}, "
            f"multiplicity={self.multiplicity})"
        )

    @property
    def num_samples(self) -> int:
        """Number of samples in the tree sequence."""
        return self.ts.num_samples

    @property
    def sequence_length(self) -> float:
        """Length of the sequence."""
        return self.ts.sequence_length

    @property
    def num_trees(self) -> int:
        """Number of local trees in the tree sequence."""
        return self.ts.num_trees

    @property
    def num_mutations(self) -> int:
        """Number of mutations in the tree sequence."""
        return self.ts.num_mutations

    def normalize_weight(self, normalization_constant: float):
        """
        Normalize the particle weight.

        Args:
            normalization_constant: The sum of all particle weights
        """
        self.weight /= normalization_constant

    def update_weight(self, likelihood: float):
        """
        Update the particle weight by multiplying with a likelihood.

        Args:
            likelihood: The likelihood to multiply the weight by
        """
        self.weight *= likelihood

    def update_log_weight(self, log_likelihood: float):
        """
        Update the particle weight using log likelihood.

        This is more numerically stable than multiplying likelihoods directly.

        Args:
            log_likelihood: The log likelihood to add to the log weight
        """
        # Convert weight to log space, add log_likelihood, convert back
        if self.weight > 0:
            log_weight = np.log(self.weight) + log_likelihood
            self.weight = np.exp(log_weight)
        else:
            self.weight = 0.0
        self.log_likelihood += log_likelihood


def create_initial_particles(
    num_particles: int,
    sample_size: int,
    sequence_length: float,
    demography: Optional[object] = None,
    recombination_rate: float = 1e-8,
    random_seed: Optional[int] = None
) -> list[Particle]:
    """
    Create initial particle ensemble by simulating from the prior.

    Each particle is initialized with a tree sequence simulated under
    the demographic model.

    Args:
        num_particles: Number of particles to create
        sample_size: Number of samples (haplotypes) in each particle
        sequence_length: Length of the sequence
        demography: msprime Demography object (None = constant population)
        recombination_rate: Recombination rate per base per generation
        random_seed: Random seed for reproducibility

    Returns:
        List of initialized Particle objects
    """
    import msprime

    particles = []

    for i in range(num_particles):
        # Set seed for each particle (deterministic but different per particle)
        seed = None if random_seed is None else random_seed + i

        # Simulate tree sequence from prior
        ts = msprime.sim_ancestry(
            samples=sample_size,
            demography=demography,
            sequence_length=sequence_length,
            recombination_rate=recombination_rate,
            random_seed=seed,
            ploidy=1  # Haploid individuals
        )

        # Create particle with equal initial weights
        particle = Particle(ts=ts, weight=1.0 / num_particles)
        particles.append(particle)

    return particles
