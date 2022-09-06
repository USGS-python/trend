"""
Tests for trends in time series data

Nonparametric Trend Tests
-------------------------
Nonparametric tests should be used where one or more of the following
conditions exist in the data.

    1. The data are nonnormal.

    2. There are missing values in the data.

    3. The data are censored. Censored data are those observations report
    at as being less than or greater than some threshold value.

In some cases, it may be desirable to remove the affect of covariates from a
trend test. To do so, the variable under consideration is adjusted
by one or several covariates by building a regression model and applying the
trend test to the residuals. A common example of this approach is
flow adjustment in hydrology, which attempts to remove the effect of natural variablity
in flow on constituents in a stream.
Nonparametric procedures cannot be applied to flow-adjusted records
containing censored data since regression residuals cannot be computed for
censored values (Hirsch et al., 1991).

References
----------
.. [1] Hirsch, R.M., R.B. Alexander, R.A. Smith. 1991. Selection of Methods
       for the Detection and Estimation of Trends in Water Quality Data.
       Water Resources Research.

Resources
----------
.. [1] https://up-rs-esp.github.io/mkt/#
.. [2] https://cran.r-project.org/web/packages/trend/vignettes/trend.pdf
"""

import pandas as pd
import numpy as np

from scipy.stats import mannwhitneyu, norm, rankdata


def sen_diff(x):
    """Sen's difference operator.

    Paramaters
    ----------
    x : array_like
        Observations taken at a fixed frequency.

    Returns
    -------
    Sen difference
    """
    #x = x[~np.isnan(x)]
    n = len(x)
    N = int(n*(n-1)/2)  # number of slope estimates
    s = np.zeros(N)
    i = 0
    for j in np.arange(1, n):
        #s[i:j+i] = (x[j] - x[0:j])/np.arange(1, j+1)
        s[i:j+i] = (x[j] - x[0:j])/np.arange(j, 0, -1)
        i += j

    return s


def sen_slope(x, alpha=None):
    """A nonparametric estimate of trend.

        Parameters
    ----------
    x : array_like
        Observations taken at a fixed frequency.

    Notes
    -----
    This method works with missing or censored data, as long as less <20% of
    observations are censored.

    References
    ----------
    .. [1] Helsel and Hirsch, R.M. 2002. Statistical Methods in Water Resources.
    .. [2] https://vsp.pnnl.gov/help/vsample/nonparametric_estimate_of_trend.htm
    """
    s = sen_diff(x)
    s.sort()

    if alpha:
        N = len(s)
        # calculate confidence limits
        C_alpha = norm.ppf(1-alpha/2)*np.sqrt(np.nanvar(x))
        U = int(np.round(1 + (N + C_alpha)/2))
        L = int(np.round((N - C_alpha)/2))
        return np.nanmedian(s), s[L], s[U]

    else:
        return np.nanmedian(s)


def seasonal_sen_slope(x, period=12, alpha=None):
    """A nonparametric estimate of trend for seasonal time series.

    Paramters
    ---------
    x : array_like
        Observations taken at a fixed frequency.

    period : int
        Number of observations in a cycle. The number of seasons. 
    """
    s = 0

    for season in np.arange(0, period):
        x_season = x[season::period]
        s = np.append(s, sen_diff(x_season))

    s.sort()

    if alpha:
        #  XXX This code needs to be verified
        N = len(s)
        # calculate confidence limits
        C_alpha = norm.ppf(1-alpha/2)*np.sqrt(np.nanvar(x))
        U = int(np.round(1 + (N + C_alpha)/2))
        L = int(np.round((N - C_alpha)/2))
        return np.nanmedian(s), s[L], s[U]

    else:
        return np.nanmedian(s)


def pettitt(x, alpha=0.05):
    """Pettitt's change-point test

    A nonparameteric test for detecting change points in a time series.

    Parameters
    ----------
    x : array_like
        Observations taken at a fixed frequency.

    alpha : float
        Significance level

    Return
    ------
    The index of the change point of the series, provided that it is
    statistically significant.
    """
    U_t = np.zeros_like(x)
    n = len(x)

    r = rankdata(x)
    for i in np.arange(n):
        U_t[i] = 2 * np.sum(r[:i+1]) - (i+1)*(n-1)

    t = np.argmax(np.abs(U_t))
    K_t = U_t[t]

    p = 2.0 * np.exp((-6.0 * K_t**2)/(n**3 + n**2))

    if p < alpha:
        return t
    else:
        return np.nan


def mk_z(s, var_s):
    """Computes the MK test statistic, Z.

    Parameters
    ----------
    s : float
        The MK trend statistic, S.

    var_s : float
        Variance of S.

    Returns
    -------
    MK test statistic, Z.
    """
    # calculate the MK test statistic
    if s > 0:
        z = (s - 1)/np.sqrt(var_s)
    elif s < 0:
        z = (s + 1)/np.sqrt(var_s)
    else:
        z = 0

    return z


def mk_score(x):
    """Computes S statistic used in Mann-Kendall tests.

    Parameters
    ----------
    x : array_like
        Chronologically ordered array of observations.

    Returns
    -------
    MK trend statistic (S).
    """
    x = x[~np.isnan(x)]
    n = len(x)
    s = 0

    for j in np.arange(1, n):
        s += np.sum(np.sign(x[j] - x[0:j]))

    return s


def mk_score_variance(x):
    """Computes corrected variance of S statistic used in Mann-Kendall tests.

    Equation 8.4 from Helsel and Hirsch (2002).

    Parameters
    ----------
    x : array_like

    Returns
    -------
    Variance of S statistic


    Note that this might be equivalent to:

        See page 728 of Hirsch and Slack

    References
    ----------
    .. [1] Helsel and Hirsch, R.M. 2002. Statistical Methods in Water Resources.
    """
    x = x[~np.isnan(x)]
    n = len(x)
    # calculate the unique data
    unique_x = np.unique(x)
    # calculate the number of tied groups
    g = len(unique_x)

    # calculate the var(s)
    if n == g:  # there is no tie
        var_s = (n*(n-1)*(2*n+5))/18

    else:  # there are some ties in data
        tp = np.zeros_like(unique_x)
        for i in range(len(unique_x)):
            tp[i] = sum(x == unique_x[i])
        var_s = (n*(n-1)*(2*n+5) - np.sum(tp*(tp-1)*(2*tp+5)))/18

    return var_s


def mann_kendall(x, alpha=0.05):
    """Mann-Kendall (MK) is a nonparametric test for monotonic trend.

    Parameters
    ----------
    x : array
        Observations taken at a fixed frequency.

    Returns
    -------
    z : float
        Normalized MK test statistic.

    Examples
    --------
    >>> x = np.random.rand(100) + np.linspace(0,.5,100)
    >>> z,p = kendall(x)


    Attribution
    -----------
    Modified from code by Michael Schramn available at
    https://github.com/mps9506/Mann-Kendall-Trend/blob/master/mk_test.py
    """
    n = len(x)

    s = mk_score(x)
    var_s = mk_score_variance(x)

    z = mk_z(s, var_s)
    # calculate the p_value
    p_value = 2*(1-norm.cdf(abs(z)))  # two tail test

    return p_value


def seasonal_mann_kendall(x, period=12):
    """ Seasonal nonparametric test for detecting a monotonic trend.

    Parameters
    ----------
    x : array
        A sequence of chronologically ordered observations with fixed
        frequency.

    period : int
        The number of observations that define period. This is the number of seasons.
    """
    # Compute the SK statistic, S, for each season
    s = 0
    var_s = 0

    for season in np.arange(period):
        x_season = x[season::period]
        s += mk_score(x_season)
        var_s += mk_score_variance(x_season)

    # Compute the SK test statistic, Z, for each season.
    z = mk_z(s, var_s)

    # calculate the p_value
    p_value = 2*(1-norm.cdf(abs(z)))  # two tail test

    return p_value


def mann_whitney(x, y, use_continuity=True, alternative=None):
    """
    Compute the Mann-Whitney rank test on samples x and y.

    Parameters
    ----------
    x, y : array_like
        Array of samples, should be one-dimensional.
    use_continuity : bool, optional
            Whether a continuity correction (1/2.) should be taken into
            account. Default is True.
    alternative : None (deprecated), 'less', 'two-sided', or 'greater'
            Whether to get the p-value for the one-sided hypothesis ('less'
            or 'greater') or for the two-sided hypothesis ('two-sided').
            Defaults to None, which results in a p-value half the size of
            the 'two-sided' p-value and a different U statistic. The
            default behavior is not the same as using 'less' or 'greater':
            it only exists for backward compatibility and is deprecated.
    Returns
    -------
    statistic : float
        The Mann-Whitney U statistic, equal to min(U for x, U for y) if
        `alternative` is equal to None (deprecated; exists for backward
        compatibility), and U for y otherwise.
    pvalue : float
        p-value assuming an asymptotic normal distribution. One-sided or
        two-sided, depending on the choice of `alternative`.
    Notes
    -----
    Use only when the number of observation in each sample is > 20 and
    you have 2 independent samples of ranks. Mann-Whitney U is
    significant if the u-obtained is LESS THAN or equal to the critical
    value of U.
    This test corrects for ties and by default uses a continuity correction.
    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Mann-Whitney_U_test
    .. [2] H.B. Mann and D.R. Whitney, "On a Test of Whether one of Two Random
           Variables is Stochastically Larger than the Other," The Annals of
           Mathematical Statistics, vol. 18, no. 1, pp. 50-60, 1947.
    """
    return mannwhitneyu(x, y, use_continuity, alternative)
