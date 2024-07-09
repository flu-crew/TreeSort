# -*- coding: utf-8 -*-
import math
import random as rnd
from dendropy import Tree, Node, Edge
from treesort.tree_indexer import TreeIndexer


def compute_rea_rate_fast(annotated_tree: Tree, evol_rate: float, ignore_top_edges=1) -> float:
    """
    Computes the reassortment rate per lineage per year: O(n) implementation.
    This method counts the number of reassortment events in the parental lineage of each strain/taxon.
    Assumes strict molecular clock with the given evolutionary rate.
    :ignore_top_edges: the longest x percent of edges will be ignored.
    """
    edge_cutoff = math.inf
    if ignore_top_edges > 0:
        edge_lengths = sorted([node.edge_length for node in annotated_tree.postorder_node_iter() if node.edge_length])
        top_percentile = int(round(len(edge_lengths) * (1.0 - ignore_top_edges / 100)))
        edge_cutoff = edge_lengths[top_percentile]

    if not TreeIndexer.is_indexed(annotated_tree):
        indexer = TreeIndexer(annotated_tree.taxon_namespace)
        indexer.index_tree(annotated_tree)

    # rea_per_strain = [0 for leaf in annotated_tree.leaf_nodes()]
    # time_per_strain = [0 for leaf in annotated_tree.postorder_node_iter()]
    rea_to_root = [0 for node in annotated_tree.preorder_node_iter()]
    time_to_root = [0 for node in annotated_tree.preorder_node_iter()]
    node: Node
    for node in annotated_tree.preorder_node_iter():
        if node is annotated_tree.seed_node:
            continue
        if node.edge and node.edge_length >= edge_cutoff:
            continue  # Skip the edge if its in the top percentile (the counts will be reset to 0 here).
        parent_ind = node.parent_node.index
        time_to_root[node.index] = time_to_root[parent_ind] + node.edge_length
        rea_to_root[node.index] = rea_to_root[parent_ind]

        edge: Edge = node.edge
        if edge.annotations.get_value('is_reassorted', '0') == '1':
            rea_annotation = edge.annotations.get_value('rea').strip('"')
            is_uncertain = all([g_str.startswith('?') for g_str in rea_annotation.split(',')])  # Is this 100% uncertain reassortment?
            if not is_uncertain:
                rea_to_root[node.index] += 1
            elif rnd.uniform(0, 1) < 0.5:  # Includes probabilistically with .5 probability. TODO: reconsider?
                rea_to_root[node.index] += 1

    rea_rates = [rea_to_root[leaf.index] / time_to_root[leaf.index] * evol_rate for leaf in annotated_tree.leaf_nodes()
                 if time_to_root[leaf.index] > 0]
    rea_rate_per_lineage_per_year = sum(rea_rates) / len(rea_rates)
    return rea_rate_per_lineage_per_year


def compute_rea_rate_simple(annotated_tree: Tree, evol_rate: float, ignore_top_edges=1) -> float:
    """
    A simpler way to compute the reassortment rate: the number of detected events divided by the total size of the tree (in years)
    :ignore_top_edges: the longest x percent of edges will not be counted for.
    """
    edge_cutoff = math.inf
    if ignore_top_edges > 0:
        edge_lengths = sorted([node.edge_length for node in annotated_tree.postorder_node_iter() if node.edge_length])
        top_percentile = int(round(len(edge_lengths) * (1.0 - ignore_top_edges / 100)))
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
