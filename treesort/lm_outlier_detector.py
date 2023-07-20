# -*- coding: utf-8 -*-
import numpy as np
from sklearn.linear_model import LinearRegression


# This is an earlier idea of using linear regression outliers for reassortment detection.
class LMOutlierDetector(object):
    trained_reg: LinearRegression
    iqd: float
    q3: float

    def __init__(self, sibling_dists_s1: np.ndarray, sibling_dists_s2: np.ndarray):
        assert len(sibling_dists_s1) >= 10
        self.sibling_dists_s1 = sibling_dists_s1
        self.sibling_dists_s2 = sibling_dists_s2

        reg: LinearRegression = LinearRegression(fit_intercept=True).fit(
            sibling_dists_s1.reshape(-1, 1), sibling_dists_s2.reshape(-1, 1))
        residuals: np.ndarray = sibling_dists_s2 - reg.predict(sibling_dists_s1.reshape(-1, 1)).reshape(-1)
        residuals.sort()
        print(residuals)
        q1 = residuals[round(len(residuals) / 4) - 1]
        q3 = residuals[round(len(residuals) * 3 / 4) - 1]
        self.iqd = q3 - q1
        self.trained_reg = reg
        self.q3 = q3
        print(f'IQD {self.iqd}, Q1 {q1}, Q3 {self.q3}')

    def get_residual(self, x: float, y: float) -> float:
        residual = y - self.trained_reg.predict(np.array([[x]], dtype=float))
        return residual[0, 0]

    def is_outlier(self, x: float, y: float, iqd_mult=2) -> bool:
        residual = self.get_residual(x, y)
        return residual >= self.q3 + self.iqd * iqd_mult
