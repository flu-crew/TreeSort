# -*- coding: utf-8 -*-
import math
import random as rnd
from typing import List, Optional
from dendropy import Tree, Node, Edge
import numpy as np
from scipy.optimize import minimize, LinearConstraint
import warnings

from treesort.tree_indexer import TreeIndexer
from treesort.helpers import sibling_distance


def compute_rea_rate_simple(annotated_tree: Tree, evol_rate: float, ignore_top_edges=1) -> float:
    """
    A simpler way to compute the reassortment rate: the number of detected events divided by the total size of the tree (in years)
    :ignore_top_edges: the longest x percent of edges will not be counted for.
    """
    edge_cutoff = math.inf
    if ignore_top_edges > 0:
        edge_lengths = sorted([node.edge_length for node in annotated_tree.postorder_node_iter() if node.edge_length])
        top_percentile = min(len(edge_lengths) - 1, int(round(len(edge_lengths) * (1.0 - ignore_top_edges / 100))))
        edge_cutoff = edge_lengths[top_percentile]

    # Compute the number of reassortment events detected (a ?-only edge counts as 0.5).
    rea_events = 0
    for node in annotated_tree.postorder_node_iter():
        if node is annotated_tree.seed_node:
            continue  # Skip the root edge
        edge: Edge = node.edge
        if edge and node.edge_length >= edge_cutoff:
            continue  # Skip the edge if its in the top percentile.
        if edge.annotations.get_value('is_reassorted', '0') == '1':
            rea_annotation = edge.annotations.get_value('rea').strip('"')
            is_uncertain = all([g_str.startswith('?') for g_str in rea_annotation.split(',')])  # Is this 100% uncertain reassortment?
            if not is_uncertain:
                rea_events += 1
            else:
                rea_events += 0.5

    # Compute the total tree length (phylogenetic diversity)
    tree_length = 0
    for node in annotated_tree.postorder_node_iter():
        if node is not annotated_tree.seed_node:
            if node.edge_length and node.edge_length >= edge_cutoff:
                continue  # Skip the edge if its in the top percentile.
            tree_length += node.edge_length

    # print(f'{rea_events}, {tree_length}, {evol_rate}')
    rea_rate_per_lineage_per_year = (rea_events / tree_length * evol_rate) if tree_length > 0 else 0.0
    return rea_rate_per_lineage_per_year


def likelihood_binary(x, rea_events, edge_lengths):
    func = 0
    if x < 1e-10:
        return np.inf
    for i in range(len(rea_events)):
        if edge_lengths[i] > 0:
            if rea_events[i] > 0:
                with warnings.catch_warnings():
                    warnings.filterwarnings('error')
                    try:
                        func -= np.log(1 - np.exp(-1 * x * edge_lengths[i]))
                    except Warning:
                        # print(x)
                        func += np.inf
            else:
                func -= (-1 * x * edge_lengths[i])
        elif rea_events[i] > 0:
            # print('+1')
            pass
    return func


def compute_rea_rate_binary_mle(annotated_tree: Tree, evol_rate: float, ref_seg_len=1700) -> Optional[float]:
    rea_events = []  # reassortment events per branches (1 - at least one event, 0 - no events).
    edge_lengths = []  # Corresponding branch lengths (the two arrays are coupled).
    processed_uncertain = set()
    node: Node
    for node in annotated_tree.postorder_node_iter():
        if node is annotated_tree.seed_node:
            continue  # Skip the root edge
        is_uncertain = False
        edge: Edge = node.edge
        if edge.annotations.get_value('is_reassorted', '0') == '1':
            rea_annotation = edge.annotations.get_value('rea').strip('"')
            is_uncertain = all(
                [g_str.startswith('?') for g_str in rea_annotation.split(',')])  # Is this 100% uncertain reassortment?
            if not is_uncertain:  # Uncertain branches are handled below.
                rea_events.append(1)
        else:
            rea_events.append(0)

        edge_length = node.edge_length
        if is_uncertain:
            # check if the sister edge was already processed.
            siblings = node.parent_node.child_nodes()
            sibling = siblings[0] if siblings[0] is not Node else siblings[1]
            if sibling not in processed_uncertain:
                # log the event over the two sister branches
                rea_events.append(1)
                edge_length = sibling_distance(node.parent_node)
                processed_uncertain.add(node)
            else:
                continue  # Skip if already processed.

        if edge_length > 1e-7:
            edge_lengths.append(edge_length / evol_rate)
        elif rea_events[-1] > 0:
            # If reassortment happened on too short of an edge, this can mess up the likelihood function.
            # Replace branch length with (1 / ref_seg_len), e.g., 1 / 1700 for HA (1 substitution).
            edge_lengths.append((1 / ref_seg_len) / evol_rate)
        else:
            edge_lengths.append(0)

    # print(len(rea_events), len(edge_lengths))
    est = compute_rea_rate_simple(annotated_tree, evol_rate, ignore_top_edges=1)
    np_est = np.array([est])
    linear_constraint = LinearConstraint([[1]], [0])
    num_est = minimize(likelihood_binary, np_est, args=(rea_events, edge_lengths), tol=1e-9,
                       constraints=[linear_constraint])
    if num_est.success:
        return num_est.x[0]
    else:
        return None
