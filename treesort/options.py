# -*- coding: utf-8 -*-
import argparse
from argparse import RawDescriptionHelpFormatter
import random
import os

import matplotlib

matplotlib.use('Agg')
from matplotlib import pyplot as plt

from treetime import TreeTime, TreeTimeError, utils
from dendropy import Tree

from treesort.helpers import parse_dates

DEFAULT_PVALUE_THRESHOLD = 0.001
DEFAULT_DEVIATION = 2
DEFAULT_METHOD = 'LOCAL'
STRAIN_NAME_REGEX = r'([ABCD](/[^/\|]+){3,5})'
EPI_ISL_REGEX = r'(EPI_ISL_\d+)'

# Program interface:
parser = argparse.ArgumentParser(description='TreeSort: fast and effective reassortment detection in '
                                             'segmented RNA viruses (primarily influenza)',
                                 formatter_class=RawDescriptionHelpFormatter)
parser._optionals.title = "Arguments"
parser.add_argument('-i', type=str, action='store', dest='descriptor',
                    help='Path to the descriptor file. The descriptor file provides paths to the alignments and '
                         'phylogenetic trees for different virus segments (see examples/)', required=True)
parser.add_argument('-o', type=str, action='store', dest='output',
                    help='Path to the output file (tree will be save in nexus format)', required=True)
parser.add_argument('-m', type=str, action='store', dest='method',
                    help='Methods are "local" (default) or "mincut". The mincut method always determines the most '
                         'parsimonious reassortment placement even in ambiguous circumstances.', required=False)
parser.add_argument('--dev', type=float, action='store', dest='deviation',
                    help='Maximum deviation from the estimated substitution rate within each segment. '
                         'The default is 2 - the substitution rate on a particular tree branch is allowed to be '
                         'twice as high or twice as low as the estimated rate. '
                         'The default value was estimated from the empirical influenza A data', required=False)
parser.add_argument('--pvalue', type=float, action='store', dest='pvalue',
                    help='The cutoff p-value for the reassortment tests: the default is 0.001 (0.1 percent). '
                         'You may want to decrease or increase this parameter depending on how stringent you want '
                         'the analysis to be', required=False)
parser.add_argument('--match-on-strain', action='store_true', dest='match_strain',
                    help='Match the names (deflines in fastas) across the segments based on the strain name. '
                         'E.g., "A/Texas/32/2021" or "A/swine/A0229832/Iowa/2021". Works for flu A, B, C, and D.'
                         'This way no preprocessing is needed to standardize the names before the analysis.')
parser.add_argument('--match-on-epi', action='store_true', dest='match_epi',
                    help='Similar to "--match-on-strain", but here segments are matched based on the "EPI_ISL_XXX" '
                         'field (if present in deflines)')
parser.add_argument('--match-on-regex', action='store', dest='match_regex',
                    help='Provide your own custom regex to match the segments across the alignments.',
                    required=False)
parser.add_argument('--no-collapse', action='store_true', dest='no_collapse',
                    help='Do not collapse near-zero length branches into multifurcations '
                         '(by default, TreeSort collapses all branches shorter than 1e-7 and then optimizes '
                         'the multifurcations).')
parser.add_argument('--equal-rates', action='store_true', dest='equal_rates',
                    help='Do not estimate molecular clock rates for different segments: assume equal rates',
                    required=False)
parser.add_argument('--clades', type=str, action='store', dest='clades_path',
                    help='Path to an output file, where clades with evidence of reassrotment will be saved',
                    required=False)


def make_outdir(descriptor_path: str) -> str:
    descriptor_path = descriptor_path.split(os.path.sep)[-1]
    if descriptor_path.count('.') > 0:
        descriptor_name = '.'.join(descriptor_path.split('.')[:-1])
    else:
        descriptor_name = descriptor_path
    i = 1
    outdir = f'treesort-{descriptor_name}-{i}'
    while os.path.exists(outdir) and i <= 50:
        i += 1
        outdir = f'treesort-{descriptor_name}-{i}'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    return outdir


def estimate_clock_rate(segment: str, tree_path: str, aln_path: str, plot=False, outdir=None) -> (float, float):
    # This code was adapted from the  'estimate_clock_model' method in the treetime/wrappers.py.
    # print(f"\tExecuting TreeTime on segment {segment}...")
    tree: Tree = Tree.get(path=tree_path, schema='newick', preserve_underscores=True)
    if len(tree.leaf_nodes()) > 1000:
        # TODO: implement Bio.Phylo tree subsampling to avoid creating temporary files
        # Downsample the tree for rate estimation
        tree_path = os.path.join(outdir, os.path.split(tree_path)[-1] + '.sample1k.tre')
        taxa_labels = [t.label for t in tree.taxon_namespace]
        random.shuffle(taxa_labels)
        subtree: Tree = tree.extract_tree_with_taxa_labels(taxa_labels[:1000])
        subtree.write(path=tree_path, schema='newick')

    dates = parse_dates(aln_path)
    try:
        timetree = TreeTime(dates=dates, tree=tree_path, aln=aln_path, gtr='JC69', verbose=-1)  # TODO: JC->GTR?
    except TreeTimeError as e:
        parser.error(f"TreeTime exception on the input files {tree_path} and {aln_path}: {e}\n "
                     f"Please make sure that the specified alignments and trees are correct.")
    timetree.clock_filter(n_iqd=3, reroot='least-squares')
    timetree.use_covariation = True
    timetree.reroot()
    timetree.get_clock_model(covariation=True)
    r_val = timetree.clock_model['r_val']
    if plot:
        timetree.plot_root_to_tip()
        plt.savefig(f'{outdir}/{segment}-treetime-clock.pdf')
    d2d = utils.DateConversion.from_regression(timetree.clock_model)
    clock_rate = round(d2d.clock_rate, 7)
    print(f"\t{segment} estimated molecular clock rate: {clock_rate} (R^2 = {round(r_val, 3)})")
    return d2d.clock_rate, r_val


# Currently requiring a tree for all segments
def parse_descriptor(path: str, outdir: str, estimate_rates=True):
    segments = []
    ref_segment = -1
    with open(path) as descriptor:
        for line in descriptor:
            line = line.strip('\n').strip()
            if line:
                tokens = [token.strip() for token in line.split(',')]
                if len(tokens) != 3:
                    parser.error(f'The descriptor file should have 3 columns: {line}')
                else:
                    seg_name, aln_path, tree_path = tokens
                    if seg_name.startswith('*'):
                        ref_segment = len(segments)
                        seg_name = seg_name[1:]
                    # seg_rate = estimate_clock_rate(seg_name, tree_path, aln_path)
                    segments.append((seg_name, aln_path, tree_path, 1))
    if len(segments) <= 1:
        parser.error('The descriptor should specify at least two segments.')
    if ref_segment < 0:
        parser.error('The descriptor should specify one of the segments as a reference segment, e.g., like "*HA".')

    print(f'Read {len(segments)} segments: {", ".join([seg[0] for seg in segments])}')
    if estimate_rates:
        print('Estimating molecular clock rates for each segment (TreeTime)...')
        for i, seg in enumerate(segments):
            seg_name, aln_path, tree_path, _ = seg
            seg_rate, r_val = estimate_clock_rate(seg_name, tree_path, aln_path, plot=True, outdir=outdir)
            segments[i] = (seg_name, aln_path, tree_path, seg_rate)
    return segments, ref_segment


def parse_args():
    args = parser.parse_args()

    # Validate PVALUE
    pval = DEFAULT_PVALUE_THRESHOLD
    if args.pvalue is not None:
        pval = args.pvalue
        if pval >= 0.5 or pval <= 1e-9:
            parser.error('the PVALUE cutoff has to be positive and cannot exceed 0.5 (50%).')
    print(f'P-value threshold for significance was set to {pval}')

    # Validate DEVIATION
    deviation = DEFAULT_DEVIATION
    if args.deviation:
        deviation = args.deviation
        if deviation > 10 or deviation < 1:
            parser.error('DEVIATION has to be in the [1, 10] interval.')

    # Validate METHOD
    method = DEFAULT_METHOD
    if args.method:
        if args.method in ('local', 'mincut'):
            method = args.method.upper()
        else:
            parser.error(f'Unknown method "{args.method}". Known methods are "local" or "mincut".')

    # Check for a MATCH regex
    match_regex = None
    if args.match_strain:
        match_regex = STRAIN_NAME_REGEX
    if args.match_epi:
        match_regex = EPI_ISL_REGEX
    if args.match_regex:
        match_regex = args.match_regex
    if match_regex:
        print(f'Using the "{match_regex}" REGEX to match segments across trees and alignments.')

    collapse_branches = False if args.no_collapse else True

    outdir = make_outdir(args.descriptor)
    segments, ref_segment = parse_descriptor(args.descriptor, outdir, not args.equal_rates)

    return args.descriptor, outdir, segments, ref_segment, args.output, args.clades_path, pval, deviation, method, \
           collapse_branches, match_regex
