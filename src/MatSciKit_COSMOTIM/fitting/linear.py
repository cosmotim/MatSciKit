"""
Weighted linear regression with error propagation.

Translated from linear_fit_with_errors.m
"""

import numpy as np
from typing import Tuple


def fit_with_errors(x: np.ndarray, y: np.ndarray, y_err: np.ndarray
                    ) -> Tuple[float, float, float]:
    """
    Weighted linear fit to data with measurement errors.

    Performs a weighted least-squares linear fit y = slope * (x - x0) + intercept,
    where weights are 1/y_err^2.

    Parameters
    ----------
    x : np.ndarray
        Independent variable values.
    y : np.ndarray
        Dependent variable values.
    y_err : np.ndarray
        Standard errors on y values.

    Returns
    -------
    slope : float
        Fitted slope.
    se_slope : float
        Standard error of the slope.
    intercept : float
        Fitted intercept (at x = mean(x)).

    Notes
    -----
    Translated from MATLAB ``linear_fit_with_errors.m``. The fit is performed
    with x shifted by its mean (x0 = mean(x)) for numerical stability.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    y_err = np.asarray(y_err, dtype=float)

    n = len(x)
    weights = 1.0 / y_err**2

    x0 = np.mean(x)
    dx = x - x0

    sum_w = np.sum(weights)
    sum_wx = np.sum(weights * dx)
    sum_wv = np.sum(weights * y)
    sum_wxx = np.sum(weights * dx**2)
    sum_wvx = np.sum(weights * y * dx)

    delta = sum_w * sum_wxx - sum_wx**2
    slope = (sum_wvx * sum_w - sum_wx * sum_wv) / delta
    intercept = (sum_wxx * sum_wv - sum_wx * sum_wvx) / delta

    # Standard error of slope
    residuals = y - slope * dx - intercept
    sum_wrr = np.sum(weights * residuals**2)
    se_slope = np.sqrt(sum_wrr * sum_w / delta / (n - 2))

    return slope, se_slope, intercept
