"""
Likelihood calculation for tree sequences using Felsenstein's pruning algorithm.

This module implements efficient likelihood calculations for genotype data
given a genealogical tree using the classic Felsenstein pruning algorithm.
"""

import numpy as np
import tskit
from typing import Optional, Union


def calculate_likelihood(
    ts: tskit.TreeSequence,
    genotypes: np.ndarray,
    mutation_rate: float,
    positions: Optional[np.ndarray] = None,
    use_log: bool = True
) -> float:
    """
    Calculate P(genotypes | tree_sequence) using Felsenstein's pruning.

    This calculates the likelihood of observed genotype data given a tree sequence
    under an infinite sites mutation model.

    Args:
        ts: Tree sequence (genealogy)
        genotypes: Array of shape (num_sites, num_samples) with genotypes (0/1)
        mutation_rate: Mutation rate per base per generation
        positions: Positions of the variant sites (default: evenly spaced)
        use_log: If True, return log likelihood (more stable)

    Returns:
        (Log) likelihood value
    """
    num_sites, num_samples = genotypes.shape

    if num_samples != ts.num_samples:
        raise ValueError(
            f"Genotypes have {num_samples} samples but tree sequence has "
            f"{ts.num_samples} samples"
        )

    # Default positions: evenly spaced along sequence
    if positions is None:
        positions = np.linspace(
            0, ts.sequence_length, num_sites, endpoint=False
        )

    log_likelihood = 0.0

    # For each variant site, find the tree and calculate likelihood
    for site_idx, pos in enumerate(positions):
        tree = ts.at(pos)
        site_genotypes = genotypes[site_idx, :]

        # Calculate likelihood for this site
        site_log_lik = calculate_site_likelihood(
            tree=tree,
            genotypes=site_genotypes,
            mutation_rate=mutation_rate,
            use_log=True
        )

        log_likelihood += site_log_lik

    if use_log:
        return log_likelihood
    else:
        return np.exp(log_likelihood)


def calculate_site_likelihood(
    tree: tskit.Tree,
    genotypes: np.ndarray,
    mutation_rate: float,
    use_log: bool = True
) -> float:
    """
    Calculate likelihood for a single site using Felsenstein's pruning.

    This implements the classic pruning algorithm that recursively calculates
    conditional likelihoods from leaves to root.

    Args:
        tree: A tskit Tree object
        genotypes: Array of length num_samples with genotypes (0/1) at leaves
        mutation_rate: Mutation rate per base per generation
        use_log: If True, return log likelihood

    Returns:
        (Log) likelihood of the genotypes given the tree
    """
    # Conditional likelihoods: L[node, state] = P(data below | node has state)
    # state 0 = ancestral, state 1 = derived
    num_nodes = tree.tree_sequence.num_nodes
    L = np.zeros((num_nodes, 2))

    # Initialize leaves with observed genotypes
    for sample_id in tree.tree_sequence.samples():
        genotype = genotypes[sample_id]
        if genotype == 0:
            L[sample_id, 0] = 1.0
            L[sample_id, 1] = 0.0
        elif genotype == 1:
            L[sample_id, 0] = 0.0
            L[sample_id, 1] = 1.0
        else:
            # Missing data or uncertain genotype
            L[sample_id, 0] = 1.0
            L[sample_id, 1] = 1.0

    # Process nodes in postorder (children before parents)
    for node in tree.nodes(order="postorder"):
        # Skip leaves (already initialized)
        if tree.is_leaf(node):
            continue

        # For internal nodes, calculate conditional likelihood
        children = tree.children(node)

        # Start with probability 1 for both states
        L[node, 0] = 1.0
        L[node, 1] = 1.0

        for child in children:
            # Get branch length
            branch_length = tree.time(node) - tree.time(child)

            # Calculate transition probabilities
            # P(no mutation) = exp(-mutation_rate * branch_length)
            # Under infinite sites: P(mutation) = 1 - P(no mutation)
            p_no_mut = np.exp(-mutation_rate * branch_length)
            p_mut = 1.0 - p_no_mut

            # For each parent state, sum over child states
            # P(child data | parent=0) = p_no_mut * L[child,0] + p_mut * L[child,1]
            # P(child data | parent=1) = p_mut * L[child,0] + p_no_mut * L[child,1]
            L[node, 0] *= (p_no_mut * L[child, 0] + p_mut * L[child, 1])
            L[node, 1] *= (p_mut * L[child, 0] + p_no_mut * L[child, 1])

    # At root, sum over possible root states (equal prior)
    root = tree.root
    likelihood = 0.5 * (L[root, 0] + L[root, 1])

    if use_log:
        if likelihood > 0:
            return np.log(likelihood)
        else:
            # Very small likelihood, return large negative log likelihood
            return -1e100
    else:
        return likelihood


def felsenstein_pruning(
    tree: tskit.Tree,
    genotypes: np.ndarray,
    mutation_rate: float,
    span: float = 1.0
) -> float:
    """
    Alias for calculate_site_likelihood for backward compatibility.

    Args:
        tree: A tskit Tree object
        genotypes: Array of genotypes at leaves
        mutation_rate: Mutation rate per base per generation
        span: Genomic span of the tree (currently not used, for compatibility)

    Returns:
        Likelihood of the genotypes given the tree
    """
    return calculate_site_likelihood(
        tree=tree,
        genotypes=genotypes,
        mutation_rate=mutation_rate,
        use_log=False
    )


def log_sum_exp(log_values: np.ndarray) -> float:
    """
    Numerically stable log-sum-exp.

    Computes log(sum(exp(log_values))) in a numerically stable way.

    Args:
        log_values: Array of log values

    Returns:
        log(sum(exp(log_values)))
    """
    max_val = np.max(log_values)
    if np.isinf(max_val):
        return max_val
    return max_val + np.log(np.sum(np.exp(log_values - max_val)))


def calculate_likelihood_at_variants(
    ts: tskit.TreeSequence,
    mutation_rate: float,
    use_log: bool = True
) -> float:
    """
    Calculate likelihood using the actual variants stored in the tree sequence.

    This is useful when the tree sequence already has mutations from simulation.

    Args:
        ts: Tree sequence with mutations
        mutation_rate: Mutation rate per base per generation
        use_log: If True, return log likelihood

    Returns:
        (Log) likelihood of the variants given the trees
    """
    if ts.num_sites == 0:
        # No variants, likelihood is 1 (log likelihood is 0)
        return 0.0 if use_log else 1.0

    log_likelihood = 0.0

    for variant in ts.variants():
        tree = ts.at(variant.site.position)

        # Calculate likelihood for this site
        site_log_lik = calculate_site_likelihood(
            tree=tree,
            genotypes=variant.genotypes,
            mutation_rate=mutation_rate,
            use_log=True
        )

        log_likelihood += site_log_lik

    if use_log:
        return log_likelihood
    else:
        return np.exp(log_likelihood)


def calculate_likelihood_parallel(
    ts: tskit.TreeSequence,
    genotypes: np.ndarray,
    mutation_rate: float,
    positions: Optional[np.ndarray] = None,
    num_threads: int = 1
) -> float:
    """
    Calculate likelihood using parallel processing (future optimization).

    For now, this just calls the serial version. Can be optimized later
    using multiprocessing or numba.

    Args:
        ts: Tree sequence
        genotypes: Genotype matrix
        mutation_rate: Mutation rate
        positions: Variant positions
        num_threads: Number of threads to use

    Returns:
        Log likelihood
    """
    # TODO: Implement parallel version
    return calculate_likelihood(ts, genotypes, mutation_rate, positions, use_log=True)
