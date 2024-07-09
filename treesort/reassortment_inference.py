# -*- coding: utf-8 -*-
import math
from dendropy import Tree, Node, Edge, DnaCharacterMatrix
from dendropy.model import parsimony
import random as rnd
from typing import List, Set, Optional, Union, Tuple

from treesort import helpers
from treesort.jc_reassortment_test import JCReassortmentTester, jc_pvalue
from treesort.helpers import node_distance_w_lca, sibling_distance_n2, binarize_tree
from treesort.parsimony import compute_parsimony_sibling_dist
from treesort.tree_indexer import TreeIndexer


REA_FIELD = 'rea_events'  # This field of Edge will store the list of inferred reassortment events.
# PROPAGATED_FROM = 'propagated_from'  # The field of Node that indicates whether the state sets were copied from a lower node.
MINCUT_STATES = 'mincut_states'  # This field of Node will store the list of 'MinCutState'


class MinCutState(object):
    node: Node
    state_sets: List[Set]
    child_states: Optional[List['MinCutState']]
    pvalue: float

    def __init__(self, node: Node, state_sets: List[Set], child_states=None, pvalue=-1):
        self.node = node
        self.state_sets = state_sets
        self.child_states = child_states
        self.pvalue = pvalue

    @classmethod
    def merge_states(cls, node: Node, child_state1: 'MinCutState', child_state2: 'MinCutState', pvalue: float):
        state_sets = []  # Fitch parsimony step up from the two children state sets.
        for site in range(len(child_state1.state_sets)):  # Perform Fitch-up.
            intersection = child_state1.state_sets[site].intersection(child_state2.state_sets[site])
            if len(intersection) == 0:
                union = child_state1.state_sets[site].union(child_state2.state_sets[site])
                state_sets.append(union)
            else:
                state_sets.append(intersection)
        return cls(node, state_sets, [child_state1, child_state2], pvalue)


class ReassortmentDetector(object):

    def __init__(self, tree: Tree, aln_path: str, segment: str, rea_tester: JCReassortmentTester, schema='fasta'):
        # tree_copy: Tree = tree.clone()  # Process the copy of the original tree.
        for node in tree.postorder_node_iter():  # Remove previous dendropy Fitch-parsimony annotations
            if hasattr(node, 'state_sets'):
                delattr(node, 'state_sets')
        self.taxon_characters: DnaCharacterMatrix = DnaCharacterMatrix.get_from_path(aln_path, schema=schema,
                                                                                taxon_namespace=tree.taxon_namespace)
        taxon_states = self.taxon_characters.taxon_state_sets_map(gaps_as_missing=True)
        parsimony.fitch_down_pass(tree.postorder_node_iter(),
                                  taxon_state_sets_map=taxon_states)  # Do the first Fitch parsimony pass (annotates the tree).
        self.tree = tree
        self.segment = segment
        self.aln_path = aln_path
        self.rea_tester = rea_tester

    def parsimony_distance(self, node1: Union[Node, MinCutState], node2: Union[Node, MinCutState]) -> int:
        parsimony_dist = 0
        for site in range(self.taxon_characters.sequence_size):
            if len(node1.state_sets[site].intersection(node2.state_sets[site])) == 0:
                parsimony_dist += 1
        return parsimony_dist

    def are_reassorted(self, node1: Union[Node, MinCutState], node2: Union[Node, MinCutState], lca: Node) -> Tuple[bool, float]:
        """
        Returns whether the results of a statistical test for reassortment between the two nodes + the pvalue.
        """
        # Compute the parsimony score between the two nodes.
        parsimony_dist = self.parsimony_distance(node1, node2)
        # Compute the ML distance on the reference tree and test for reassortment.
        if isinstance(node1, MinCutState) and isinstance(node2, MinCutState):
            ml_distance = node_distance_w_lca(node1.node, node2.node, lca)
        else:
            ml_distance = node_distance_w_lca(node1, node2, lca)
        return self.rea_tester.is_reassorted(parsimony_dist, ml_distance)

    def merge_siblings(self, sib1: Node, sib2: Node) -> Node:
        parent: Node = sib1.parent_node
        parent.remove_child(sib1)
        parent.remove_child(sib2)

        new_node = Node(edge_length=0)
        new_node.set_child_nodes([sib1, sib2])
        parent.add_child(new_node)
        return new_node

    def propagate_parsimony(self, node1: Node, node2: Node, lca: Node):
        lca_state_sets = []
        for site in range(self.taxon_characters.sequence_size):
            intersection = node1.state_sets[site].intersection(node2.state_sets[site])
            if len(intersection) == 0:
                union = node1.state_sets[site].union(node2.state_sets[site])
                lca_state_sets.append(union)
            else:
                lca_state_sets.append(intersection)
        lca.state_sets = lca_state_sets

    # def propagate_to_parent(self, node: Node):
    #     parent_state_sets = [s.copy() for s in node.state_sets]  # Copy state sets to the parent node.
    #     node.parent_node.state_sets = parent_state_sets
    #     setattr(node.parent_node, PROPAGATED_FROM, node)  # Set the PROPAGATED_FROM field.

    def add_rea_annotation(self, edge: Edge, parsimony_dist: int, is_uncertain: bool):
        # node = edge.head_node
        # while getattr(node, PROPAGATED_FROM, None):  # Go down to the lowest non-propagated node (if needed).
        #     node = getattr(node, PROPAGATED_FROM)
        # edge = node.edge
        annotation = f'{"?" if is_uncertain else ""}{self.segment}({parsimony_dist})'
        rea_events = getattr(edge, REA_FIELD, [])
        rea_events.append(annotation)
        setattr(edge, REA_FIELD, rea_events)

    def binarize_tree_greedy(self):
        """
        Greedily resolve multifurcations by grouping siblings into non-reassortant groups.
        """
        node: Node
        for node in self.tree.postorder_node_iter():
            if node.child_nodes() and len(node.child_nodes()) > 2:
                # Found a multifurcation.
                siblings = node.child_nodes()
                rnd.shuffle(siblings)  # Randomly shuffle all the siblings.

                new_siblings = [siblings[0]]  # We will group all the siblings into reassortment-free blocks greedily.
                for sibling in siblings[1:]:
                    placed = False
                    for i, new_sibling in enumerate(new_siblings):
                        reassorted, pvalue = self.are_reassorted(new_sibling, sibling, node)
                        if not reassorted:
                            merged_node = self.merge_siblings(new_sibling, sibling)
                            self.propagate_parsimony(new_sibling, sibling, merged_node)
                            new_siblings[i] = merged_node
                            placed = True
                            break
                    if not placed:
                        new_siblings.append(sibling)

                # for sibling in node.child_nodes():  # Remove all the children from "node".
                #     node.remove_child(sibling)

                if len(new_siblings) == 1:
                    # No reassortment: split out the top and add as children to the node.
                    merged_node = new_siblings[0]
                    node.state_sets = merged_node.state_sets
                    node.set_child_nodes(merged_node.child_nodes())
                else:
                    # Merge reassortant blocks in a caterpillar structure in the reverse order of their size
                    new_siblings.sort(key=lambda new_sib: len(new_sib.leaf_nodes()), reverse=True)
                    while len(new_siblings) > 2:
                        new1, new2 = new_siblings.pop(), new_siblings.pop()
                        merged_node = self.merge_siblings(new1, new2)
                        self.propagate_parsimony(new1, new2, merged_node)
                        new_siblings.append(merged_node)
                    node.set_child_nodes(new_siblings)
                    self.propagate_parsimony(new_siblings[0], new_siblings[1], node)

    def infer_reassortment_mincut(self) -> int:
        """
        MinCut algorithm that annotates the branches of the tree with the inferred reassortment events.
        This algorithm cuts the tree into the smallest number of reassortment-free subtrees possible.
        :return: The number of inferred reassortment events.
        """
        tree_indexer = TreeIndexer(self.tree.taxon_namespace)
        tree_indexer.index_tree(self.tree)
        node: Node  # every node will get a 'mincut_states' list.
        for node in self.tree.leaf_nodes():
            setattr(node, MINCUT_STATES, [MinCutState(node, node.state_sets)])  # Initialize a new mincutstate for the leaf.
        for node in self.tree.postorder_internal_node_iter():
            child1, child2 = node.child_nodes()
            compatible_pairs: List[Tuple[MinCutState, MinCutState, float]] = []  # (left_state, right_state, pvalue)
            for left_state in getattr(child1, MINCUT_STATES):
                best_match: Tuple[MinCutState, float] = (None, math.inf)  # right_state and pvalue.
                for right_state in getattr(child2, MINCUT_STATES):
                    reassorted, pvalue = self.are_reassorted(left_state, right_state, node)
                    if not reassorted and abs(0.5 - pvalue) < abs(0.5 - best_match[1]):  # check if pvalue is closer to median.
                        best_match = (right_state, pvalue)
                if best_match[0]:
                    compatible_pairs.append((left_state, best_match[0], best_match[1]))
            if compatible_pairs:
                mincut_states = []
                for left_state, right_state, pvalue in compatible_pairs:
                    # Merge the pairs.
                    mincut_states.append(MinCutState.merge_states(node, left_state, right_state, pvalue))
            else:
                # Get the union of mincut_states of children.
                mincut_states: List = getattr(child1, MINCUT_STATES).copy()
                mincut_states.extend(getattr(child2, MINCUT_STATES))
            setattr(node, MINCUT_STATES, mincut_states)
        rea_events = self._backtrack_mincut()
        print(f'\tInferred reassortment events with {self.segment}: {rea_events}.')
        return rea_events

    def _backtrack_mincut(self) -> int:
        """
        A backtracking subroutine for 'infer_reassortment_mincut'.
        Works top-down and identifies branches with reassortment.
        :return: The number of identified reassortment events.
        """
        rea_events = 0
        root: Node = self.tree.seed_node
        node: Node
        for node in self.tree.preorder_node_iter():
            if node == root:
                # Choose the option with the largest number of leaf nodes (and highest pvalue, if equal) as the most likely ancestral state.
                root_states: List[MinCutState] = getattr(node, MINCUT_STATES)
                sorted_states = sorted(root_states, key=lambda state: (len(state.node.leaf_nodes()), -abs(0.5 - state.pvalue)), reverse=True)
                major_state = sorted_states[0]
                setattr(node, MINCUT_STATES, major_state)
            else:
                parent_state: MinCutState = getattr(node.parent_node, MINCUT_STATES)
                node_states: List[MinCutState] = getattr(node, MINCUT_STATES)
                if node_states[0].node is node:
                    node_states.sort(key=lambda state: abs(0.5 - state.pvalue))  # Sort the states by pvalue closeness to 0.5
                    # Find the node's states that agree with the parent assignment:
                    relevant_states = [state for state in node_states if state in parent_state.child_states]
                    if relevant_states:
                        # Assign the "relevant" state with the highest pvalue.
                        best_state = relevant_states[0]
                    else:
                        # Cut off the edge (add rea annotation) and assign the state with the highest pvalue.
                        rea_events += 1
                        best_state = node_states[0]
                        pars_dist = self.parsimony_distance(parent_state, best_state)
                        self.add_rea_annotation(node.edge, pars_dist, is_uncertain=False)
                    setattr(node, MINCUT_STATES, best_state)
                else:
                    # Just propagate the parent state below.
                    setattr(node, MINCUT_STATES, parent_state)
        return rea_events

    def infer_reassortment_local(self, pval_threshold: float, add_uncertain=True):
        """
        The first (local) implementation, where the reassortment placement is determined by the aunt node.
        If reassortment placement in unclear, both branches get marked as potential reassortment events.
        """
        tree_indexer = TreeIndexer(self.tree.taxon_namespace)
        tree_indexer.index_tree(self.tree)
        child_dists_s2, child1_dists_s2, child2_dists_s2 = compute_parsimony_sibling_dist(self.tree, self.aln_path)
        node_by_index = {}
        node: Node
        for node in self.tree.postorder_node_iter():
            node_by_index[node.index] = node

        pvalues = [(node.index, self.rea_tester.is_reassorted(child_dists_s2[node.index], helpers.sibling_distance(node))[1])
                   for node in self.tree.internal_nodes()]
        jc_outliers = [(index, pvalue) for index, pvalue in pvalues if pvalue < pval_threshold]
        jc_outlier_indices = [x[0] for x in jc_outliers]
        # print(len(jc_outliers))
        total_rea, certain_rea = 0, 0
        for outlier_ind, pvalue in jc_outliers:
            outlier_node: Node = node_by_index[outlier_ind]
            specific_edge: Node = None
            c1, c2 = outlier_node.child_nodes()
            annotation = f'{self.segment}({child_dists_s2[outlier_ind]})'
            if outlier_node != self.tree.seed_node:
                c1_pvalue = self.rea_tester.is_reassorted(child1_dists_s2[outlier_ind], helpers.aunt_distance(c1))[1]
                c2_pvalue = self.rea_tester.is_reassorted(child2_dists_s2[outlier_ind], helpers.aunt_distance(c2))[1]
                c1_outlier = c1_pvalue < pval_threshold
                c2_outlier = c2_pvalue < pval_threshold
                if (not c1_outlier) and (not c2_outlier):
                    # print('Neither', child_dists_s2[outlier_ind], helpers.sibling_distance(outlier_node))
                    continue
                if c1_outlier and c2_outlier:
                    # print('Both', child_dists_s2[outlier_ind], helpers.sibling_distance(outlier_node))
                    pass
                if c1_outlier ^ c2_outlier:
                    specific_edge = c1 if c1_outlier else c2
            total_rea += 1
            if specific_edge:
                certain_rea += 1
                edge: Edge = specific_edge.edge
                self.add_rea_annotation(edge, child_dists_s2[outlier_ind], is_uncertain=False)
            elif add_uncertain:
                for edge in [c1.edge, c2.edge]:
                    self.add_rea_annotation(edge, child_dists_s2[outlier_ind], is_uncertain=True)
        print(f'\tInferred reassortment events with {self.segment}: {total_rea}.\n'
              f'\tIdentified exact branches for {certain_rea}/{total_rea} of them')
