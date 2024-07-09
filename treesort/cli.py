# -*- coding: utf-8 -*-
import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from dendropy import Tree, Node, Edge

from treesort import options, helpers
from treesort.helpers import binarize_tree, collapse_zero_branches
from treesort.tree_indexer import TreeIndexer
from treesort.reassortment_utils import compute_rea_rate_simple
from treesort.reassortment_inference import REA_FIELD, ReassortmentDetector
from treesort.jc_reassortment_test import JCReassortmentTester

ADD_UNCERTAIN = True  # For local method only.
# METHOD = 'LOCAL'  # MINCUT or LOCAL.
RESOLVE_GREEDY = True  # Whether to use the greedy multifurcation reslution algorithm
# COLLAPSE_ZERO_BRANCHES = True


# TODO: taxa should be identified by strain name only (then substituted before writing the output tree).
#  a regex can be used to specify a capture pattern


def run_treesort_cli():
    # Each segment has format (name, aln_path, tree_path, rate)
    sys.setrecursionlimit(100000)
    descriptor_name, segments, ref_segment_i, output_path, clades_out_path, pval_threshold, allowed_deviation, \
        method, collapse_branches = options.parse_args()
    ref_tree_path = segments[ref_segment_i][2]
    tree: Tree = Tree.get(path=ref_tree_path, schema='newick', preserve_underscores=True)
    ref_seg = segments[ref_segment_i]
    if collapse_branches:
        collapse_zero_branches(tree, 1e-7)

    if RESOLVE_GREEDY:
        print('Optimally resolving the multifurcations to minimize reassortment...')
        tree.suppress_unifurcations()  # remove the unifurcations.
        # Read all alignments and compute the averaged sub rate.
        aln_by_seg = []
        total_rate = 0
        total_sites = 0
        for i, seg in enumerate(segments):
            if i == ref_segment_i:
                aln_by_seg.append(None)
                continue
            seg_list = list(SeqIO.parse(seg[1], format='fasta'))
            seg_aln = {seq.id: seq.seq for seq in seg_list}
            aln_len = len(seg_list[0].seq)
            total_rate += seg[3] * aln_len
            total_sites += aln_len
            aln_by_seg.append(seg_aln)
        overall_rate = total_rate / total_sites
        # print('Concatenated rate:', overall_rate)

        # Concatenate non-ref alignments.
        concatenated_seqs = []
        taxa = [leaf.taxon.label for leaf in tree.leaf_nodes()]
        for taxon in taxa:
            # Find this taxon across all (non-reference) segments and concatenate aligned sequences.
            concat_seq = ''
            for i, seg in enumerate(segments):
                if i == ref_segment_i:
                    continue
                concat_seq += aln_by_seg[i][taxon]
            concatenated_seqs.append(SeqRecord(concat_seq, id=taxon, name=taxon, description=''))
        concat_path = descriptor_name + '.concatenated.fasta'
        SeqIO.write(concatenated_seqs, concat_path, 'fasta')

        # Binarize the tree
        reassortment_tester = JCReassortmentTester(total_sites, overall_rate / ref_seg[3], pval_threshold, allowed_deviation)
        rea_detector = ReassortmentDetector(tree, concat_path, 'concatenated', reassortment_tester)
        rea_detector.binarize_tree_greedy()
        # tree = rea_detector.tree  # use the binarized tree
    else:
        binarize_tree(tree)  # simple binarization, where polytomies are resolved as caterpillars.

    tree_indexer = TreeIndexer(tree.taxon_namespace)
    tree_indexer.index_tree(tree)

    for i, seg in enumerate(segments):
        if i == ref_segment_i:
            continue

        print(f'Inferring reassortment with the {seg[0]} segment...')
        seg2_aln = list(SeqIO.parse(seg[1], format='fasta'))
        seg2_len = len(seg2_aln[0])
        rate_ratio = seg[3] / ref_seg[3]
        reassortment_tester = JCReassortmentTester(seg2_len, rate_ratio, pval_threshold, allowed_deviation)
        rea_detector = ReassortmentDetector(tree, seg[1], seg[0], reassortment_tester)

        if method == 'MINCUT':
            rea_detector.infer_reassortment_mincut()
        else:
            rea_detector.infer_reassortment_local(pval_threshold, add_uncertain=ADD_UNCERTAIN)

    clades_out = None
    reported_rea = set()
    if clades_out_path:
        clades_out = open(clades_out_path, 'w')

    for node in tree.postorder_node_iter():
        annotation = ','.join(getattr(node.edge, REA_FIELD, []))
        if annotation:
            # if len(node.leaf_nodes()) >= 20:
            #     leaf = node.leaf_nodes()[0]
            #     # print(annotation, len(node.leaf_nodes()), leaf.taxon.label)  # TODO: comment out
            # else:
            #     # print(annotation, len(node.leaf_nodes()), ';'.join([leaf.taxon.label for leaf in node.leaf_nodes()]))
            #     pass
            # edge: Edge = node.edge
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
                if report_genes_str:
                    clades_out.write(f'{clade},{report_genes_str}'
                                     # if there are ?-genes - add alternative clade (sister clade)
                                     f'{"," + sister_clade if report_genes_str.count("?") > 0 else ","}\n')
        else:
            node.edge.annotations.add_new('is_reassorted', '0')

    if clades_out:
        clades_out.close()

    rea_rate = compute_rea_rate_simple(tree, ref_seg[3], ignore_top_edges=1)
    print(f'Estimated reassortment rate per lineage per year: {round(rea_rate, 6)}')

    tree.write_to_path(output_path, schema='nexus')
    # tree.write_to_path(output_path + 'phylo.xml', schema='phyloxml')
    print(f'Saved the annotated tree file to {output_path}')
