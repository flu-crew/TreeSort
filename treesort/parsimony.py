# -*- coding: utf-8 -*-
from typing import Tuple

import numpy as np
from dendropy import Tree, Node, DnaCharacterMatrix
from dendropy.model import parsimony
import random as rnd

# character_sets_annotation = 'character_sets'


def compute_parsimony_edge_lengths(tree: Tree, aln_path: str) -> np.ndarray:
    """
    Compute the parsimony score of the tree given an alignment and find associated edge-lengths.
    :param tree: Tree topology to be scored by parsimony. Must be BINARY.
    :param aln_path: DNA alignment for the tips of the tree
    :return: A dictionary that specifies # of parsimony substitutions per node (except the root).
    """
    tree_copy: Tree = tree.clone()
    taxon_characters: DnaCharacterMatrix = DnaCharacterMatrix.get_from_path(aln_path, schema='fasta',
                                                                            taxon_namespace=tree_copy.taxon_namespace)
    taxon_states = taxon_characters.taxon_state_sets_map(gaps_as_missing=True)
    p_score = parsimony.fitch_down_pass(tree_copy.postorder_node_iter(), taxon_state_sets_map=taxon_states)
    print(p_score)

    edge_lengths = np.zeros(len(tree.nodes()), dtype=int)
    p_score_2 = 0
    node: Node
    for node in tree_copy.preorder_node_iter():
        edge_len = 0
        parent: Node = node.parent_node
        for site in range(taxon_characters.sequence_size):
            if parent:
                parent_state = parent.state_sets[site]
                if parent_state in node.state_sets[site]:
                    node.state_sets[site] = parent_state
                    continue
                else:
                    edge_len += 1
                    p_score_2 += 1
            # choose a random state and assign
            state_sets = list(node.state_sets[site])
            rnd_state = rnd.choice(state_sets)
            node.state_sets[site] = rnd_state
        if parent:
            cluster = {leaf.taxon.label for leaf in node.leaf_nodes()}
            original_node = tree.find_node(filter_fn=lambda tree_node: {leaf.taxon.label for leaf in tree_node.leaf_nodes()} == cluster)
            edge_lengths[original_node.index] = edge_len
    print(p_score, p_score_2)
    return edge_lengths


def compute_parsimony_sibling_dist(tree: Tree, aln_path: str, schema='fasta') -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    # tree must be binary!
    tree_copy: Tree = tree.clone()
    taxon_characters: DnaCharacterMatrix = DnaCharacterMatrix.get_from_path(aln_path, schema=schema,
                                                                            taxon_namespace=tree_copy.taxon_namespace)
    taxon_states = taxon_characters.taxon_state_sets_map(gaps_as_missing=True)
    p_score = parsimony.fitch_down_pass(tree_copy.postorder_node_iter(), taxon_state_sets_map=taxon_states)
    # print(p_score)

    children_dists = np.zeros(len(tree.internal_nodes()), dtype=int)
    child1_to_sibling_dists, child2_to_sibling_dists = np.zeros(len(tree.internal_nodes()), dtype=int),\
                                                       np.zeros(len(tree.internal_nodes()), dtype=int)
    p_score_2 = 0
    node: Node
    for node in tree_copy.preorder_internal_node_iter():
        children_dist = 0
        child1_to_sibling, child2_to_sibling = 0, 0
        sibling = node.sibling_nodes()[0] if node.parent_node else None
        if not sibling:
            child1_to_sibling, child2_to_sibling = -1, -1
        child1, child2 = node.child_nodes()
        for site in range(taxon_characters.sequence_size):
            if len(child1.state_sets[site].intersection(child2.state_sets[site])) == 0:
                children_dist += 1
                p_score_2 += 1
            if sibling and len(child1.state_sets[site].intersection(sibling.state_sets[site])) == 0:
                child1_to_sibling += 1
            if sibling and len(child2.state_sets[site].intersection(sibling.state_sets[site])) == 0:
                child2_to_sibling += 1
        # cluster = {leaf.taxon.label for leaf in node.leaf_nodes()}
        # original_node = tree.find_node(filter_fn=lambda tree_node: {leaf.taxon.label for leaf in tree_node.leaf_nodes()} == cluster)
        children_dists[node.index] = children_dist
        child1_to_sibling_dists[node.index] = child1_to_sibling
        child2_to_sibling_dists[node.index] = child2_to_sibling
    # print(p_score, p_score_2)
    # TODO: if p_score != p_score_2: log a warning (debug only)
    return children_dists, child1_to_sibling_dists, child2_to_sibling_dists


def get_cluster_str(node: Node) -> str:
    return ';'.join(sorted([leaf.taxon.label for leaf in node.leaf_nodes()]))

#
# if __name__ == '__main__':
#     # seg1_tree_path = '../simulations/segs2/l1500/sim_250_10/sim_1.trueSeg1.tre'
#     # schema = 'nexus'
#     # seg1_path = '../simulations/segs2/l1500/sim_250_10/sim_1.seg1.alignment.fasta'
#     # seg2_path = '../simulations/segs2/l1500/sim_250_10/sim_1.seg2.alignment.fasta'
#     # simulated = True
#     seg1_tree_path = '../../gammas/HAs.fast.rooted.tre'
#     schema = 'newick'
#     seg1_path = '../../gammas/HAs_unique.aln'
#     seg2_path = '../../gammas/NAs_unique.aln'
#     simulated = False
#     na_ha_ratio = 1.057
#
#     tree: Tree = Tree.get(path=seg1_tree_path, schema=schema, preserve_underscores=True)
#     binarize_tree(tree)  # Randomly binarize.
#     tree_indexer = TreeIndexer(tree.taxon_namespace)
#     tree_indexer.index_tree(tree)
#     if simulated:
#         node: Node
#         for node in tree.nodes():
#             if node.edge_length:
#                 node.edge_length *= 0.00474
#
#     # lengths_by_node_s1 = compute_parsimony_edge_lengths(tree, 'testdata/l1000_50_5/sim_1.seg4.alignment.fasta')
#     # lengths_by_node_s2 = compute_parsimony_edge_lengths(tree, 'testdata/l1000_50_5/sim_1.seg1.alignment.fasta')
#     child_dists_s1, child1_dists_s1, child2_dists_s1 = compute_parsimony_sibling_dist(tree, seg1_path)
#     child_dists_s2, child1_dists_s2, child2_dists_s2 = compute_parsimony_sibling_dist(tree, seg2_path)
#     seg2_aln = list(SeqIO.parse(seg2_path, format='fasta'))
#     seg2_len = len(seg2_aln[0])
#
#     node_by_index = {}
#     cluster_by_index = {}
#     for node in tree.postorder_node_iter():
#         node_by_index[node.index] = node
#         cluster_by_index[node.index] = get_cluster_str(node)
#
#     # s1_lengths = np.zeros(len(lengths_by_node_s1))
#     # node: Node
#     # for node in tree.postorder_internal_node_iter():
#     #     cluster = [leaf.taxon.label for leaf in node.leaf_nodes()]
#     #     # print(cluster, lengths_by_node_s1[node.index], lengths_by_node_s2[node.index])
#     #     print(node.index, sorted(cluster),
#     #           f'c1 {node.child_nodes()[0].index}({child1_dists_s1[node.index]}, {child1_dists_s2[node.index]})',
#     #           f'c2 {node.child_nodes()[1].index}({child2_dists_s1[node.index]}, {child2_dists_s2[node.index]})')
#
#     outlier_detector = LMOutlierDetector(child_dists_s1, child_dists_s2)
#     outliers = [(ind, outlier_detector.get_residual(child_dists_s1[ind], child_dists_s2[ind]))
#                 for ind in range(child_dists_s1.size) if
#                 outlier_detector.is_outlier(child_dists_s1[ind], child_dists_s2[ind], iqd_mult=3)]
#     jc_outliers = [(node.index, -1) for node in tree.internal_nodes() if
#                    is_jc_outlier(child_dists_s2[node.index], seg2_len, helpers.sibling_distance(node),
#                                  rate_ratio=na_ha_ratio, pvalue_threshold=0.001)]
#     jc_outlier_indices = [x[0] for x in jc_outliers]
#     outliers = sorted(outliers, key=lambda x: x[1], reverse=True)
#     print(len(outliers))
#     print(len(jc_outliers))
#     for outlier_ind, residual in jc_outliers:
#         outlier_node = node_by_index[outlier_ind]
#         is_c1_rea = outlier_detector.is_outlier(child1_dists_s1[outlier_ind], child1_dists_s2[outlier_ind], iqd_mult=1.7)
#         is_c2_rea = outlier_detector.is_outlier(child2_dists_s1[outlier_ind], child2_dists_s2[outlier_ind], iqd_mult=1.7)
#         print(outlier_ind, residual, is_c1_rea, cluster_by_index[outlier_node.child_nodes()[0].index],
#               is_c2_rea, cluster_by_index[outlier_node.child_nodes()[1].index])
#
#     jc_colors = ['red' if i in jc_outlier_indices else 'blue' for i in range(len(tree.internal_nodes()))]
#     plt.scatter(child_dists_s1, child_dists_s2, c=jc_colors)
#     for ind in range(len(child_dists_s1)):
#         plt.annotate(str(ind), (child_dists_s1[ind], child_dists_s2[ind] + 0.2))
#     plt.show()
