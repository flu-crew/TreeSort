# -*- coding: utf-8 -*-
import sys

from Bio import SeqIO
from dendropy import Tree, Node, Edge

from treesort import options, helpers
from treesort.helpers import binarize_tree
from treesort.jc_outlier_detector import is_jc_outlier, jc_pvalue
from treesort.parsimony import compute_parsimony_sibling_dist
from treesort.tree_indexer import TreeIndexer

PVALUE_THRESHOLD = 0.0001
ADD_UNCERTAIN = True

# TODO: taxa should be identified by strain name only (then substituted before writing the output tree).
#  a regex can be used to specify a capture pattern


def run_treesort_cli():
    # Each segment has format (name, aln_path, tree_path, rate)
    sys.setrecursionlimit(100000)
    segments, ref_segment_i, output_path, clades_out_path = options.parse_args()
    ref_tree_path = segments[ref_segment_i][2]
    tree: Tree = Tree.get(path=ref_tree_path, schema='newick', preserve_underscores=True)
    binarize_tree(tree)  # Randomly binarize.
    tree_indexer = TreeIndexer(tree.taxon_namespace)
    tree_indexer.index_tree(tree)

    ref_seg = segments[ref_segment_i]
    # child_dists_s1, child1_dists_s1, child2_dists_s1 = compute_parsimony_sibling_dist(tree, ref_seg[1])

    node_by_index = {}
    node: Node
    for node in tree.postorder_node_iter():
        node_by_index[node.index] = node

    edge_annotations = {}
    for i, seg in enumerate(segments):
        if i == ref_segment_i:
            continue

        print(f'Inferring reassortment with the {seg[0]} segment...')
        child_dists_s2, child1_dists_s2, child2_dists_s2 = compute_parsimony_sibling_dist(tree, seg[1])
        seg2_aln = list(SeqIO.parse(seg[1], format='fasta'))
        seg2_len = len(seg2_aln[0])
        rate_ratio = seg[3] / ref_seg[3]

        # print('Segment rate ratio: %.5f' % rate_ratio)

        pvalues = [(node.index, jc_pvalue(child_dists_s2[node.index], seg2_len, helpers.sibling_distance(node),
                                          rate_ratio=rate_ratio)) for node in tree.internal_nodes()]
        jc_outliers = [(index, pvalue) for index, pvalue in pvalues if pvalue < PVALUE_THRESHOLD]
        jc_outlier_indices = [x[0] for x in jc_outliers]
        # print(len(jc_outliers))
        total_rea, certain_rea = 0, 0
        for outlier_ind, pvalue in jc_outliers:
            outlier_node: Node = node_by_index[outlier_ind]
            specific_edge: Node = None
            c1, c2 = outlier_node.child_nodes()
            annotation = f'{seg[0]}({child_dists_s2[outlier_ind]})'
            if outlier_node != tree.seed_node:
                c1_pvalue = jc_pvalue(child1_dists_s2[outlier_ind], seg2_len, helpers.aunt_distance(c1),
                                      rate_ratio=rate_ratio)
                c2_pvalue = jc_pvalue(child2_dists_s2[outlier_ind], seg2_len, helpers.aunt_distance(c2),
                                      rate_ratio=rate_ratio)
                c1_outlier = c1_pvalue < PVALUE_THRESHOLD
                c2_outlier = c2_pvalue < PVALUE_THRESHOLD
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
                node_id = edge.head_node.index
                rea_annotation = edge_annotations.get(node_id, '')
                rea_annotation = annotation if not rea_annotation else rea_annotation + ',' + annotation
                edge_annotations[node_id] = rea_annotation
            elif ADD_UNCERTAIN:
                for edge in [c1.edge, c2.edge]:
                    node_id = edge.head_node.index
                    rea_annotation = edge_annotations.get(node_id, '')
                    rea_annotation = '?' + annotation if not rea_annotation else rea_annotation + ',?' + annotation
                    edge_annotations[node_id] = rea_annotation
        print(f'\tInferred {ref_seg[0]}-{seg[0]} reassortment events: {total_rea}.\n'
              f'\tIdentified exact branches for {certain_rea}/{total_rea} of them')

    clades_out = None
    reported_rea = set()
    if clades_out_path:
        clades_out = open(clades_out_path, 'w')

    for node in tree.postorder_node_iter():
        annotation = edge_annotations.get(node.index, None)
        if annotation:
            if len(node.leaf_nodes()) >= 20:
                leaf = node.leaf_nodes()[0]
                # print(annotation, len(node.leaf_nodes()), leaf.taxon.label)  # TODO: comment out
            else:
                # print(annotation, len(node.leaf_nodes()), ';'.join([leaf.taxon.label for leaf in node.leaf_nodes()]))
                pass
            edge: Edge = node.edge
            node.edge.annotations.add_new('rea', f'"{annotation}"')
            node.edge.annotations.add_new('is_reassorted', '1')

            if clades_out:
                # Report reassortment associated with the clade.
                clade = ';'.join(sorted([leaf.taxon.label for leaf in node.leaf_nodes()]))
                sister_clade = ';'.join(sorted([leaf.taxon.label for leaf in node.sister_nodes()[0].leaf_nodes()]))
                reassorted_genes = [g_str[:g_str.find('(')] for g_str in annotation.split(',')]
                # Drop out already reported ?-genes: we want to report ?-genes only once.
                # TODO: ideally report ?-genes with the more likely clade
                report_genes = [gene for gene in reassorted_genes if (sister_clade, gene) not in reported_rea]
                for rea_gene in report_genes:
                    if rea_gene.startswith('?'):
                        reported_rea.add((clade, rea_gene))  # mark the ?-genes that we report here
                report_genes_str = ';'.join(report_genes)
                clades_out.write(f'{clade},{report_genes_str}'
                                 # if there are ?-genes - add alternative clade (sister clade)
                                 f'{"," + sister_clade if report_genes_str.count("?") > 0 else ","}\n')
        else:
            node.edge.annotations.add_new('is_reassorted', '0')

    if clades_out:
        clades_out.close()

    tree.write_to_path(output_path, schema='nexus')
    # tree.write_to_path(output_path + 'phylo.xml', schema='phyloxml')
    print(f'Saved the annotated tree file to {output_path}')
