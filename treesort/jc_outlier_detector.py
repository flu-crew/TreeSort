# -*- coding: utf-8 -*-
from scipy.stats import binomtest
import math


def is_jc_outlier(subs: int, sites: int, ml_distance: float, rate_ratio=1, pvalue_threshold=0.01,
                  allowed_deviation=1.7):
    """
    We assume the Jukes-Cantor substitution model and test whether the observed number of substitutions was likely
    to come from the observed time interval (ml_distance). The method assumes strict molecular clock.
    :param subs: Number of observed substitutions in the second gene segment
    :param sites: Number of sites in the second gene segment
    :param ml_distance: Expected number of substitutions per site in the first gene segment
    :param rate_ratio: Ratio in global substitution rates between the second and first segments
    :param pvalue_threshold: p-values below this threshold will be inferred as reassortments
    :param allowed_deviation: Should be >1: allowed deviation from the strict molecular clock in the second segment.
    :return: True if reassortment is likely
    """
    sub_probability = 0.75 - 0.75 * (math.exp(-(4 * ml_distance * rate_ratio * allowed_deviation)/3))
    pvalue = binomtest(subs, sites, p=sub_probability, alternative='greater').pvalue
    # print(subs, sites, ml_distance, sub_probability, pvalue)
    return pvalue < pvalue_threshold


def jc_pvalue(subs: int, sites: int, ml_distance: float, rate_ratio=1, allowed_deviation=1.7):
    if ml_distance < 1 / sites:
        ml_distance = 1 / sites
    sub_probability = 0.75 - 0.75 * (math.exp(-(4 * ml_distance * rate_ratio * allowed_deviation) / 3))
    pvalue = binomtest(subs, sites, p=sub_probability, alternative='greater').pvalue
    # if pvalue < 0.001:
    #     print(subs, sites, ml_distance, sub_probability, pvalue)
    return pvalue
