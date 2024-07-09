# -*- coding: utf-8 -*-
from scipy.stats import binomtest
import math
from typing import Tuple


def jc_pvalue(subs: int, sites: int, ml_distance: float, rate_ratio=1.0, allowed_deviation=1.5):
    """
    We assume the Jukes-Cantor substitution model and test whether the observed number of substitutions was likely
    to come from the observed time interval (ml_distance). The method assumes strict molecular clock (but with deviation)
    :param subs: Number of observed substitutions in the second gene segment
    :param sites: Number of sites in the second gene segment
    :param ml_distance: Expected number of substitutions per site in the first gene segment
    :param rate_ratio: Ratio in global substitution rates between the second and first segments
    :param pvalue_threshold: p-values below this threshold will be inferred as reassortments
    :param allowed_deviation: Should be >=1: allowed deviation from the strict molecular clock in each segment
    :return: the pvalue of observing the number of subs over the ml_distance edge.
    """
    if ml_distance < 1 / sites:
        ml_distance = 1 / sites
    max_deviation = allowed_deviation * allowed_deviation
    sub_probability = 0.75 - 0.75 * (math.exp(-(4 * ml_distance * rate_ratio * max_deviation) / 3))
    pvalue = binomtest(subs, sites, p=sub_probability, alternative='greater').pvalue
    # if pvalue < 0.001:
    #     print(subs, sites, ml_distance, sub_probability, pvalue)
    return pvalue


class JCReassortmentTester(object):

    def __init__(self, sites: int, rate_ratio: float, pvalue_threshold: float, allowed_deviation: float):
        self.sites = sites
        self.rate_ratio = rate_ratio
        self.pvalue_threshold = pvalue_threshold
        self.allowed_deviation = allowed_deviation

    def is_reassorted(self, subs: int, ml_distance: float) -> Tuple[bool, float]:
        pvalue = jc_pvalue(subs, self.sites, ml_distance, self.rate_ratio, self.allowed_deviation)
        return (pvalue < self.pvalue_threshold), pvalue
