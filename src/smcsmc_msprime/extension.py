"""
ARG extension methods for particle propagation.

This module handles extending genealogies (ARGs) as particles are propagated
along the genome. This is the most algorithmically challenging part of SMCSMC.
"""

import numpy as np
import msprime
import tskit
from typing import Optional, Tuple
from .particle import Particle


def extend_particle_simple(
    particle: Particle,
    start_pos: float,
    end_pos: float,
    demography: msprime.Demography,
    recombination_rate: float,
    mutation_rate: Optional[float] = None,
    random_seed: Optional[int] = None
) -> Particle:
    """
    Extend particle using the simple "simulate-then-observe" approach.

    This is the simplest extension strategy: if the particle's tree sequence
    doesn't cover the target region, we assume it was simulated with the full
    chromosome and just use it directly.

    For the initial prototype, we assume particles were initialized with
    tree sequences covering the full sequence length.

    Args:
        particle: Particle to extend
        start_pos: Start position
        end_pos: End position
        demography: Demographic model
        recombination_rate: Recombination rate
        mutation_rate: Mutation rate (optional, for adding mutations)
        random_seed: Random seed

    Returns:
        Extended particle (may be same particle if already covers region)
    """
    # Check if particle already covers this region
    if particle.sequence_length >= end_pos:
        # Already covers the region, no extension needed
        return particle

    # If not, we have a problem - this simple approach assumes full simulation
    raise ValueError(
        f"Particle sequence length ({particle.sequence_length}) "
        f"does not cover target position ({end_pos}). "
        f"For simple extension, initialize particles with full sequence length."
    )


def extend_particle_incremental(
    particle: Particle,
    start_pos: float,
    end_pos: float,
    demography: msprime.Demography,
    recombination_rate: float,
    mutation_rate: Optional[float] = None,
    random_seed: Optional[int] = None
) -> Particle:
    """
    Extend particle incrementally using msprime's initial_state parameter.

    This is a more realistic approach that extends the ARG incrementally.
    EXPERIMENTAL - may not work correctly yet.

    Args:
        particle: Particle to extend
        start_pos: Start position
        end_pos: End position
        demography: Demographic model
        recombination_rate: Recombination rate
        mutation_rate: Optional mutation rate
        random_seed: Random seed

    Returns:
        Extended particle with updated tree sequence and importance weight
    """
    # TODO: Implement proper incremental extension
    # This would use msprime.sim_ancestry with initial_state parameter
    # or a custom extension mechanism

    raise NotImplementedError(
        "Incremental ARG extension not yet implemented. "
        "Use simple extension (simulate full chromosome upfront) for now."
    )


def calculate_extension_weight(
    old_ts: tskit.TreeSequence,
    new_ts: tskit.TreeSequence,
    start_pos: float,
    end_pos: float
) -> float:
    """
    Calculate importance sampling weight for ARG extension.

    When using biased proposal distributions (e.g., importance sampling
    on recombination rate), we need to calculate the weight adjustment.

    Args:
        old_ts: Tree sequence before extension
        new_ts: Tree sequence after extension
        start_pos: Start position of extension
        end_pos: End position of extension

    Returns:
        Importance weight (ratio of target to proposal density)
    """
    # For simple extension with no biased sampling, weight is 1.0
    # In full SMCSMC, this would calculate the ratio of probabilities
    # under the target vs. proposal distributions

    # TODO: Implement proper importance weight calculation
    return 1.0


def extract_recombination_events(
    ts: tskit.TreeSequence,
    start_pos: float,
    end_pos: float
) -> list:
    """
    Extract recombination events from a tree sequence in a region.

    Recombination events are identified by edges that start/end within
    the specified region.

    Args:
        ts: Tree sequence
        start_pos: Start of region
        end_pos: End of region

    Returns:
        List of recombination events (position, time, nodes)
    """
    events = []

    for edge in ts.edges():
        # Check if edge starts or ends in our region
        if start_pos <= edge.left < end_pos or start_pos < edge.right <= end_pos:
            # This edge represents a recombination event
            events.append({
                'left': edge.left,
                'right': edge.right,
                'parent': edge.parent,
                'child': edge.child,
                'time': ts.node(edge.parent).time
            })

    return events


def extract_coalescence_events(
    ts: tskit.TreeSequence,
    start_pos: float,
    end_pos: float
) -> list:
    """
    Extract coalescence events from a tree sequence.

    Coalescence events are identified by internal nodes and their children.

    Args:
        ts: Tree sequence
        start_pos: Start of region
        end_pos: End of region

    Returns:
        List of coalescence events (time, nodes)
    """
    events = []

    # For each tree in the region
    for tree in ts.trees():
        if tree.interval.left >= end_pos:
            break
        if tree.interval.right <= start_pos:
            continue

        # For each internal node
        for node in tree.nodes():
            if not tree.is_leaf(node) and node != tree.root:
                children = list(tree.children(node))
                if len(children) >= 2:
                    events.append({
                        'position': tree.interval.left,
                        'time': tree.time(node),
                        'node': node,
                        'children': children,
                        'num_children': len(children)
                    })

    return events


def simplify_particle(
    particle: Particle,
    keep_samples: Optional[list] = None
) -> Particle:
    """
    Simplify a particle's tree sequence to reduce memory usage.

    This removes ancient parts of the genealogy that are no longer
    relevant for likelihood calculations.

    Args:
        particle: Particle to simplify
        keep_samples: Sample nodes to keep (default: all samples)

    Returns:
        New particle with simplified tree sequence
    """
    if keep_samples is None:
        keep_samples = particle.ts.samples()

    # Simplify the tree sequence
    simplified_ts = particle.ts.simplify(samples=keep_samples)

    # Create new particle with simplified tree sequence
    new_particle = Particle(
        ts=simplified_ts,
        weight=particle.weight,
        pilot_weight=particle.pilot_weight,
        multiplicity=particle.multiplicity
    )

    return new_particle


def trim_particle(
    particle: Particle,
    max_time: float
) -> Particle:
    """
    Trim a particle's tree sequence to remove very ancient nodes.

    This can reduce memory usage by removing genealogical information
    beyond a certain time depth.

    Args:
        particle: Particle to trim
        max_time: Maximum time to keep

    Returns:
        New particle with trimmed tree sequence
    """
    # TODO: Implement tree sequence trimming by time
    # This would involve filtering nodes and edges by time
    # Not critical for initial prototype

    raise NotImplementedError("Tree sequence trimming not yet implemented")


class ParticleExtender:
    """
    Helper class for extending particles with various strategies.
    """

    def __init__(
        self,
        demography: msprime.Demography,
        recombination_rate: float,
        mutation_rate: Optional[float] = None,
        strategy: str = "simple"
    ):
        """
        Initialize the extender.

        Args:
            demography: Demographic model
            recombination_rate: Recombination rate
            mutation_rate: Mutation rate (optional)
            strategy: Extension strategy ("simple" or "incremental")
        """
        self.demography = demography
        self.recombination_rate = recombination_rate
        self.mutation_rate = mutation_rate
        self.strategy = strategy

    def extend(
        self,
        particle: Particle,
        start_pos: float,
        end_pos: float,
        random_seed: Optional[int] = None
    ) -> Particle:
        """
        Extend a particle from start_pos to end_pos.

        Args:
            particle: Particle to extend
            start_pos: Start position
            end_pos: End position
            random_seed: Random seed

        Returns:
            Extended particle
        """
        if self.strategy == "simple":
            return extend_particle_simple(
                particle=particle,
                start_pos=start_pos,
                end_pos=end_pos,
                demography=self.demography,
                recombination_rate=self.recombination_rate,
                mutation_rate=self.mutation_rate,
                random_seed=random_seed
            )
        elif self.strategy == "incremental":
            return extend_particle_incremental(
                particle=particle,
                start_pos=start_pos,
                end_pos=end_pos,
                demography=self.demography,
                recombination_rate=self.recombination_rate,
                mutation_rate=self.mutation_rate,
                random_seed=random_seed
            )
        else:
            raise ValueError(f"Unknown extension strategy: {self.strategy}")
