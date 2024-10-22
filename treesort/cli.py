# -*- coding: utf-8 -*-
import sys
from typing import List, Optional, Dict, Tuple, Set
import re
import os

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from dendropy import Tree, Node, Edge, Taxon

from treesort import options, helpers
from treesort.helpers import binarize_tree, collapse_zero_branches
from treesort.tree_indexer import TreeIndexer
from treesort.reassortment_utils import compute_rea_rate_simple
from treesort.reassortment_inference import REA_FIELD, ReassortmentDetector
from treesort.jc_reassortment_test import JCReassortmentTester

ADD_UNCERTAIN = True  # For local method only.
RESOLVE_GREEDY = True  # Whether to use the greedy multifurcation resolution algorithm
MIN_TAXA_THRESHOLD = 10


def extract_join_regex(label: str, join_on_regex: str, print_error=True) -> Optional[str]:
    """
    Extracts the portion of the strain label captured by 'join_on_regex'.
    """
    re_search = re.search(join_on_regex, label)
    if re_search and re_search.group(0):
        return re_search.group(0)
    else:
        if print_error:
            print(f'Cannot match pattern {join_on_regex} to {label}. Skipping this strain.')
        return None


def is_taxon_in_tree(label: str, tree: Tree, join_on_regex=None) -> bool:
    if join_on_regex:
        label = extract_join_regex(label, join_on_regex)
        tree_labels = {extract_join_regex(leaf.taxon.label, join_on_regex) for leaf in tree.leaf_nodes()}
    else:
        tree_labels = {leaf.taxon.label for leaf in tree.leaf_nodes()}
    if label:
        return label in tree_labels
    else:
        return False


def get_aln_labels(aln: Dict[str, SeqRecord], join_on_regex=None) -> Set[str]:
    """
    Get a set of strain labels from the alignment
    (returns the label portion matched by 'join_on_regex' only, if specified).
    """
    if join_on_regex:
        return {extract_join_regex(strain, join_on_regex) for strain in aln.keys()}
    else:
        return set(aln.keys())


def is_taxon_in_aln(label: str, aln_labels: Set[str], join_on_regex=None) -> bool:
    if join_on_regex:
        label = extract_join_regex(label, join_on_regex)
    if label:
        return label in aln_labels
    else:
        return False


def find_common_taxa(aln_by_seg: List[Dict[str, SeqRecord]], ref_segment_i: int, join_on_regex=None) -> List[str]:
    # Find taxa in common.
    aln_labels_by_seg = [get_aln_labels(aln, join_on_regex) for aln in aln_by_seg]
    common_taxa = [extract_join_regex(strain, join_on_regex) if join_on_regex else strain
                   for strain in aln_by_seg[ref_segment_i] if
                   all([is_taxon_in_aln(strain, aln_labels, join_on_regex) for aln_labels in aln_labels_by_seg])]
    return common_taxa


def prune_tree_to_taxa(tree: Tree, common_taxa: List[str], join_on_regex=None) -> Optional[Dict[str, str]]:
    """
    Prune the tree and rename the taxa according to 'join_on_regex' regex, if provided.
    Returns a dictionary that maps the new names to the old names (if subs were made).
    """
    name_map: Optional[Dict[str, str]] = None
    if join_on_regex:
        # Need to rename all taxa first.
        name_map = {}
        new_taxa = set()
        taxon: Taxon
        for taxon in tree.taxon_namespace:
            new_label = extract_join_regex(taxon.label, join_on_regex, print_error=False)
            if new_label:
                if new_label not in new_taxa:
                    name_map[new_label] = taxon.label
                    taxon.label = new_label
                    new_taxa.add(new_label)
                else:
                    print(f'REPEATED strain {new_label} - discarding the copy')

    # Prune the tree.
    tree.retain_taxa_with_labels(common_taxa)
    return name_map


def prune_and_update_alignments(aln_by_seg: List[Dict[str, SeqRecord]], segments: List[Tuple[str, str, str, float]],
                                common_taxa: List[str], outdir: str, join_on_regex=None) -> List[Dict[str, SeqRecord]]:
    if not join_on_regex:
        # Don't need to do anything.
        return aln_by_seg
    else:
        upd_aln_by_seg: List[Dict[str, SeqRecord]] = []
        for i, seg in enumerate(segments):
            aln_map = aln_by_seg[i]
            new_aln: List[SeqRecord] = []
            upd_aln_map: Dict[str, SeqRecord] = {}
            added_labels_upper = set()
            upd_aln_by_seg.append(upd_aln_map)
            new_aln_path = os.path.join(outdir, f'{seg[0]}_unified.aln')
            for label in aln_map:
                new_label = extract_join_regex(label, join_on_regex, print_error=False)
                if new_label and new_label in common_taxa:
                    if new_label.upper() in added_labels_upper:
                        # print(f'REPEATED strain {new_label} in segment {seg[0]} - discarding')
                        continue  # do not add to the alignment
                    else:
                        record = aln_map[label]
                        record.id = record.name = new_label
                        record.description = ''
                        new_aln.append(record)
                        upd_aln_map[new_label] = record
                        added_labels_upper.add(new_label.upper())
            SeqIO.write(new_aln, new_aln_path, 'fasta')
            segments[i] = (seg[0], new_aln_path, seg[2], seg[3])
        return upd_aln_by_seg


def run_treesort_cli():
    # Each segment has format (name, aln_path, tree_path, rate)
    sys.setrecursionlimit(100000)
    segments: List[Tuple[str, str, str, float]]  # name, aln path, tree path, rate.
    descriptor_name, outdir, segments, ref_segment_i, output_path, clades_out_path, pval_threshold, allowed_deviation, \
        method, collapse_branches, join_on_regex = options.parse_args()
    ref_tree_path = segments[ref_segment_i][2]
    tree: Tree = Tree.get(path=ref_tree_path, schema='newick', preserve_underscores=True)
    ref_seg = segments[ref_segment_i]
    if collapse_branches:
        collapse_zero_branches(tree, 1e-7)

    # Parse the alignments into a list of dictionaries.
    aln_by_seg: List[Dict[str, SeqRecord]] = []
    for i, seg in enumerate(segments):
        seg_list = list(SeqIO.parse(seg[1], format='fasta'))
        seg_aln = {seq.id: seq for seq in seg_list}
        aln_by_seg.append(seg_aln)

    # Find taxa in common and prune trees/alignments if needed.
    common_taxa = find_common_taxa(aln_by_seg, ref_segment_i, join_on_regex)
    name_map: Optional[Dict[str, str]] = None # New to old label names map (if subs were made)
    if len(common_taxa) >= MIN_TAXA_THRESHOLD:
        print(f'Found {len(common_taxa)} strains in common across the alignments.')
        name_map = prune_tree_to_taxa(tree, common_taxa, join_on_regex)
        aln_by_seg = prune_and_update_alignments(aln_by_seg, segments, common_taxa, outdir, join_on_regex)
    else:
        # Print an error and exit.
        print(f'Found {len(common_taxa)} strains in common across the segment alignments - insufficient for a reassortment analysis.')
        exit(-1)

    if RESOLVE_GREEDY:
        print('Optimally resolving the multifurcations to minimize reassortment...')
        tree.suppress_unifurcations()  # remove the unifurcations.
        # Compute the averaged sub rate.
        total_rate = 0
        total_sites = 0
        for i, seg in enumerate(segments):
            if i == ref_segment_i:
                continue
            aln_len = len(next(iter(aln_by_seg[i].values())).seq)
            total_rate += seg[3] * aln_len
            total_sites += aln_len
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
                concat_seq += aln_by_seg[i][taxon].seq
            concatenated_seqs.append(SeqRecord(concat_seq, id=taxon, name=taxon, description=''))
        concat_path = os.path.join(outdir, descriptor_name + '.concatenated.fasta')
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
        # seg2_aln = list(SeqIO.parse(seg[1], format='fasta'))
        seg2_len = len(next(iter(aln_by_seg[i].values())).seq)
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
            #     # print(annotation, len(node.leaf_nodes()), leaf.taxon.label)
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

    if ref_seg[3] < 1:  # Do not estimate the reassortment rate if --equal-rates was given.
        rea_rate = compute_rea_rate_simple(tree, ref_seg[3], ignore_top_edges=1)
        print(f'Estimated reassortment rate per lineage per year: {round(rea_rate, 6)}')

    if name_map:
        # Substitute the labels back.
        for leaf in tree.leaf_node_iter():
            if leaf.taxon:
                leaf.taxon.label = name_map[leaf.taxon.label]

    tree.write_to_path(output_path, schema='nexus')
    # tree.write_to_path(output_path + 'phylo.xml', schema='phyloxml')
    print(f'Saved the annotated tree file to {output_path}')
