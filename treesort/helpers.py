# -*- coding: utf-8 -*-
from datetime import datetime
import re

from Bio import SeqIO
from dendropy import Tree, Node


def parse_dates(aln_path: str):
    records = SeqIO.parse(aln_path, 'fasta')
    dates = {}
    for record in records:
        name = record.name
        # date_str = name.split('|')[-1]
        date = None
        for token in name.split('|'):
            if re.fullmatch(r'[\d\-/]{4,}', token) and not re.fullmatch(r'\d{5,}', token):
                if token.count('/') == 2:
                    date = datetime.strptime(token, '%m/%d/%Y')
                elif token.count('/') == 1:
                    date = datetime.strptime(token, '%m/%Y')
                elif token.count('-') == 2:
                    date = datetime.strptime(token, '%Y-%m-%d')
                elif token.count('-') == 1:
                    date = datetime.strptime(token, '%Y-%m')
                else:
                    date = datetime.strptime(token, '%Y')
        dec_date = date.year + ((date.month - 1) * 30 + date.day) / 365.0
        dates[name] = dec_date
    return dates


def sibling_distance(parent_node: Node) -> float:
    return parent_node.child_nodes()[0].edge_length + parent_node.child_nodes()[1].edge_length


def aunt_distance(node: Node) -> float:
    # We assume that the tree is binary
    assert node.parent_node and node.parent_node.parent_node
    parent: Node = node.parent_node
    aunt: Node = parent.sibling_nodes()[0]
    return node.edge_length + parent.edge_length + aunt.edge_length


def node_distance(node1: Node, node2: Node) -> float:
    """
    Linear-time algorithm to find a distance between two nodes on the same tree.
    Note: with constant-time LCA computation, one can compute distance in constant time.
    """
    node1_depth = get_node_depth(node1)
    node2_depth = get_node_depth(node2)
    distance = 0
    p1, p2 = node1, node2
    if node1_depth > node2_depth:
        for step in range(node1_depth - node2_depth):
            distance += p1.edge_length
            p1 = p1.parent_node
    elif node2_depth > node1_depth:
        for step in range(node2_depth - node1_depth):
            distance += p2.edge_length
            p2 = p1.parent_node

    while p1 != p2:
        distance += p1.edge_length
        distance += p2.edge_length
        p1 = p1.parent_node
        p2 = p2.parent_node
    return distance


def get_node_depth(node: Node) -> int:
    depth = 0
    p: Node = node.parent_node
    while p:
        depth += 1
        p = p.parent_node
    return depth


def binarize_tree(tree: Tree, edge_length=0):
    """
    Adds/removes nodes from the tree to make it fully binary (added edges will have length 'edge_length')
    :param tree: Dendropy tree to be made bifurcating.
    """

    # First suppress unifurcations.
    tree.suppress_unifurcations()

    # Now binarize multifurcations.
    node: Node
    for node in tree.postorder_node_iter():
        if node.child_nodes() and len(node.child_nodes()) > 2:
            num_children = len(node.child_nodes())
            children = node.child_nodes()
            interim_node = node
            # Creates a caterpillar structure with children on the left of the trunk:
            for child_ind in range(len(children) - 2):
                new_node = Node(edge_length=edge_length)
                interim_node.set_child_nodes([children[child_ind], new_node])
                interim_node = new_node
            interim_node.set_child_nodes(children[num_children - 2:])
