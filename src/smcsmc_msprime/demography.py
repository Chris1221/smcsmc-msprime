"""
Demographic model specification for SMCSMC.

This module provides helper functions for creating and working with
demographic models using msprime's Demography API.
"""

import msprime
import numpy as np
from typing import List, Optional, Tuple, Dict


class DemographyHelper:
    """
    Helper class for creating and manipulating demographic models.
    """

    @staticmethod
    def constant_size(Ne: float, population_name: str = "pop0") -> msprime.Demography:
        """
        Create a constant population size demographic model.

        Args:
            Ne: Effective population size
            population_name: Name of the population

        Returns:
            msprime Demography object
        """
        demography = msprime.Demography()
        demography.add_population(name=population_name, initial_size=Ne)
        return demography

    @staticmethod
    def exponential_growth(
        Ne_present: float,
        growth_rate: float,
        population_name: str = "pop0"
    ) -> msprime.Demography:
        """
        Create an exponentially growing/declining population model.

        Args:
            Ne_present: Present-day effective population size
            growth_rate: Exponential growth rate (positive = growth, negative = decline)
            population_name: Name of the population

        Returns:
            msprime Demography object
        """
        demography = msprime.Demography()
        demography.add_population(
            name=population_name,
            initial_size=Ne_present,
            growth_rate=growth_rate
        )
        return demography

    @staticmethod
    def piecewise_constant(
        sizes: List[float],
        change_times: List[float],
        population_name: str = "pop0"
    ) -> msprime.Demography:
        """
        Create a piecewise constant population size model.

        The population has different constant sizes in different epochs.

        Args:
            sizes: List of population sizes (most recent first)
            change_times: List of times (in generations ago) when size changes
            population_name: Name of the population

        Returns:
            msprime Demography object

        Example:
            # Population of size 10000 from present to 1000 generations ago,
            # then size 5000 from 1000 to 5000 generations ago,
            # then size 20000 before that
            demography = piecewise_constant(
                sizes=[10000, 5000, 20000],
                change_times=[1000, 5000]
            )
        """
        if len(sizes) != len(change_times) + 1:
            raise ValueError(
                f"Need {len(change_times) + 1} sizes for {len(change_times)} "
                f"change times, got {len(sizes)}"
            )

        demography = msprime.Demography()
        demography.add_population(
            name=population_name,
            initial_size=sizes[0]
        )

        # Add population size changes
        for time, size in zip(change_times, sizes[1:]):
            demography.add_population_parameters_change(
                time=time,
                initial_size=size,
                population=population_name
            )

        return demography

    @staticmethod
    def bottleneck(
        Ne_present: float,
        Ne_bottleneck: float,
        bottleneck_start: float,
        bottleneck_end: float,
        Ne_ancestral: Optional[float] = None,
        population_name: str = "pop0"
    ) -> msprime.Demography:
        """
        Create a population bottleneck model.

        Args:
            Ne_present: Present-day population size
            Ne_bottleneck: Population size during bottleneck
            bottleneck_start: Time (generations ago) when bottleneck started
            bottleneck_end: Time (generations ago) when bottleneck ended
            Ne_ancestral: Ancestral population size (default: same as present)
            population_name: Name of the population

        Returns:
            msprime Demography object
        """
        if Ne_ancestral is None:
            Ne_ancestral = Ne_present

        return DemographyHelper.piecewise_constant(
            sizes=[Ne_present, Ne_bottleneck, Ne_ancestral],
            change_times=[bottleneck_start, bottleneck_end],
            population_name=population_name
        )

    @staticmethod
    def two_epoch(
        Ne_present: float,
        Ne_ancestral: float,
        change_time: float,
        population_name: str = "pop0"
    ) -> msprime.Demography:
        """
        Create a simple two-epoch model.

        Args:
            Ne_present: Recent population size
            Ne_ancestral: Ancestral population size
            change_time: Time of size change (generations ago)
            population_name: Name of the population

        Returns:
            msprime Demography object
        """
        return DemographyHelper.piecewise_constant(
            sizes=[Ne_present, Ne_ancestral],
            change_times=[change_time],
            population_name=population_name
        )

    @staticmethod
    def from_epochs(
        epochs: List[Dict[str, float]],
        population_name: str = "pop0"
    ) -> msprime.Demography:
        """
        Create a demographic model from a list of epochs.

        Args:
            epochs: List of dicts with 'start_time', 'end_time', and 'Ne'
            population_name: Name of the population

        Returns:
            msprime Demography object

        Example:
            epochs = [
                {'start_time': 0, 'end_time': 1000, 'Ne': 10000},
                {'start_time': 1000, 'end_time': 5000, 'Ne': 5000},
                {'start_time': 5000, 'end_time': float('inf'), 'Ne': 20000}
            ]
        """
        # Sort epochs by start time
        sorted_epochs = sorted(epochs, key=lambda x: x['start_time'])

        sizes = [epoch['Ne'] for epoch in sorted_epochs]
        change_times = [epoch['start_time'] for epoch in sorted_epochs[1:]]

        return DemographyHelper.piecewise_constant(
            sizes=sizes,
            change_times=change_times,
            population_name=population_name
        )


def get_epoch_boundaries(demography: msprime.Demography) -> List[float]:
    """
    Extract the epoch boundaries from a demographic model.

    Args:
        demography: msprime Demography object

    Returns:
        List of epoch boundary times (sorted)
    """
    # Get all demographic events
    events = demography.events
    times = [0.0] + sorted(set(event.time for event in events))
    return times


def get_population_size_at_time(
    demography: msprime.Demography,
    time: float,
    population_name: str = "pop0"
) -> float:
    """
    Get the population size at a specific time.

    Args:
        demography: msprime Demography object
        time: Time (generations ago)
        population_name: Name of the population

    Returns:
        Population size at the specified time
    """
    # This is a simplified implementation
    # In practice, would need to trace through all demographic events
    pop = demography[population_name]
    return pop.initial_size


def create_human_demographic_model() -> msprime.Demography:
    """
    Create a simple human demographic model (approximate).

    This is a simplified model with:
    - Recent exponential growth
    - Constant ancestral size

    Returns:
        msprime Demography object
    """
    # Approximate human demography
    # Recent size ~10^6, ancestral size ~10^4
    # Growth started ~400 generations ago
    return DemographyHelper.piecewise_constant(
        sizes=[1e6, 1e4],
        change_times=[400],
        population_name="human"
    )


def rescale_demography(
    demography: msprime.Demography,
    scaling_factor: float
) -> msprime.Demography:
    """
    Rescale all population sizes in a demographic model.

    Useful for sensitivity analysis or changing the reference Ne.

    Args:
        demography: Original demography
        scaling_factor: Multiply all population sizes by this factor

    Returns:
        New demography with rescaled sizes
    """
    # This is a placeholder - proper implementation would need to
    # deep copy the demography and rescale all sizes
    raise NotImplementedError(
        "Demography rescaling not yet implemented. "
        "Create a new demography with scaled sizes instead."
    )


# Convenience functions
def constant_Ne(Ne: float = 10000, **kwargs) -> msprime.Demography:
    """Shorthand for constant size demographic model."""
    return DemographyHelper.constant_size(Ne, **kwargs)


def exponential(Ne_present: float, growth_rate: float, **kwargs) -> msprime.Demography:
    """Shorthand for exponential growth model."""
    return DemographyHelper.exponential_growth(Ne_present, growth_rate, **kwargs)


def bottleneck(
    Ne_present: float,
    Ne_bottleneck: float,
    bottleneck_start: float,
    bottleneck_end: float,
    **kwargs
) -> msprime.Demography:
    """Shorthand for bottleneck model."""
    return DemographyHelper.bottleneck(
        Ne_present, Ne_bottleneck, bottleneck_start, bottleneck_end, **kwargs
    )


def two_epoch(
    Ne_present: float,
    Ne_ancestral: float,
    change_time: float,
    **kwargs
) -> msprime.Demography:
    """Shorthand for two-epoch model."""
    return DemographyHelper.two_epoch(Ne_present, Ne_ancestral, change_time, **kwargs)
