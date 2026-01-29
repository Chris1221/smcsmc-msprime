"""
SMCSMC-msprime: Efficient demographic inference using tree sequences.

A Python implementation of SMCSMC (Sequential Monte Carlo for population genetics)
using msprime and tskit for scalable genealogical inference.
"""

__version__ = "0.1.0"

# Core classes
from .particle import Particle, create_initial_particles
from .particle_filter import (
    ParticleFilter,
    ParticleFilterConfig,
    ParticleFilterResults,
    run_particle_filter
)

# Likelihood calculation
from .likelihood import (
    calculate_likelihood,
    calculate_site_likelihood,
    felsenstein_pruning,
    calculate_likelihood_at_variants
)

# Resampling
from .resampling import (
    systematic_resampling,
    multinomial_resampling,
    stratified_resampling,
    residual_resampling,
    effective_sample_size,
    should_resample,
    normalize_weights,
    compute_resampling_stats
)

# Demography helpers
from .demography import (
    DemographyHelper,
    constant_Ne,
    exponential,
    bottleneck,
    two_epoch,
    get_epoch_boundaries,
    create_human_demographic_model
)

# Extension methods
from .extension import (
    extend_particle_simple,
    extract_recombination_events,
    extract_coalescence_events,
    simplify_particle,
    ParticleExtender
)

# Inference
from .inference import (
    InferenceResults,
    extract_coalescence_events,
    extract_coalescence_times_per_tree,
    estimate_Ne_from_coalescence_times,
    estimate_Ne_from_tmrca,
    estimate_Ne_from_pairwise_coalescence,
    infer_constant_Ne,
    compare_to_true_parameters,
    summarize_inference
)

__all__ = [
    # Version
    "__version__",
    # Core classes
    "Particle",
    "create_initial_particles",
    "ParticleFilter",
    "ParticleFilterConfig",
    "ParticleFilterResults",
    "run_particle_filter",
    # Likelihood
    "calculate_likelihood",
    "calculate_site_likelihood",
    "felsenstein_pruning",
    "calculate_likelihood_at_variants",
    # Resampling
    "systematic_resampling",
    "multinomial_resampling",
    "stratified_resampling",
    "residual_resampling",
    "effective_sample_size",
    "should_resample",
    "normalize_weights",
    "compute_resampling_stats",
    # Demography
    "DemographyHelper",
    "constant_Ne",
    "exponential",
    "bottleneck",
    "two_epoch",
    "get_epoch_boundaries",
    "create_human_demographic_model",
    # Extension
    "extend_particle_simple",
    "extract_recombination_events",
    "simplify_particle",
    "ParticleExtender",
    # Inference
    "InferenceResults",
    "extract_coalescence_times_per_tree",
    "estimate_Ne_from_coalescence_times",
    "estimate_Ne_from_tmrca",
    "estimate_Ne_from_pairwise_coalescence",
    "infer_constant_Ne",
    "compare_to_true_parameters",
    "summarize_inference",
]
