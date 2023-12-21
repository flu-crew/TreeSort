# -*- coding: utf-8 -*-
import random as rnd
from dendropy import Tree, Node, Edge
from treesort.tree_indexer import TreeIndexer


def compute_rea_rate(annotated_tree: Tree, evol_rate: float) -> float:
    """
    Computes the reassortment rate per lineage per year.
    Assumes strict molecular clock with the given evolutionary rate.
    """
    if not TreeIndexer.is_indexed(annotated_tree):
        indexer = TreeIndexer(annotated_tree.taxon_namespace)
        indexer.index_tree(annotated_tree)

    rea_per_strain = [0 for leaf in annotated_tree.leaf_nodes()]
    time_per_strain = [0 for leaf in annotated_tree.postorder_node_iter()]
    node: Node
    for node in annotated_tree.postorder_node_iter():
        if node is annotated_tree.seed_node:
            continue  # Skip the root edge
        edge: Edge = node.edge
        if edge.annotations.get_value('is_reassorted', '0') == '1':
            rea_annotation = edge.annotations.get_value('rea').strip('"')
            is_uncertain = all([g_str.startswith('?') for g_str in rea_annotation.split(',')])  # Is this 100% uncertain reassortment?
            if not is_uncertain:  # TODO: include probabilistically?
                for leaf in node.leaf_nodes():
                    rea_per_strain[leaf.taxon.index] += 1
        for leaf in node.leaf_nodes():
            time_per_strain[leaf.taxon.index] += edge.length
    rea_rates = [rea_per_strain[leaf.taxon.index] / time_per_strain[leaf.taxon.index] * evol_rate for
                 leaf in annotated_tree.leaf_nodes()]
    rea_rate_per_lineage_per_year = sum(rea_rates) / len(rea_rates)
    return rea_rate_per_lineage_per_year
