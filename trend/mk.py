import pandas as pd
import numpy as np

from scipy.stats import norm
from trend import mk_score, mk_score_variance

class MKT:
    @staticmethod
    def _z(s, var_s):
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
        return (s-np.sign(s)) / np.sqrt(var_s)
    
    @staticmethod
    def _s(x):
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
    
    @staticmethod
    def _s_variance(x):
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

class RKT:
    @classmethod
    def test(cls, df, block=None, cv=None, correct=False, rep="e"):
        """Regional Mann-Kendal Test

        Parameters
        ----------
        df : DataFrame
            Dateframe with datetime index and columns for each parameter to test

        block : str
            A column in df containing positive integers representing blocks, i.e. sites,
            seasons, or months, or a combination.

        cv : str
            Name of covariable column in df

        correct : boolean
            If False, no correction for correlation between blocks is performed.

        rep : str
            TODO

        Attribution
        -----------
        Aldo Marchetto's rkt package
        """
        df = df.copy()
        # Check that df is correctly formated

        # if correcting for a covariable, eliminate records where missing
        df = df.dropna(subset=[cv])

        # if no block given create a dummy block, could do the same for cv

        # replace values sharing the same date with median or average

        # don't elimate records with missing values for y, rather use a threshold
        # below which to eliminate/ignore records with too few values

        # make blocks consecutive integers from 1 to maxblock
        max_block = df[block].max()
        for i in range(max_block):
            if(df[df[block]==i].empty):
                df[block>=i] = df[block>=i] - 1

        # initialize valriables
        tau = np.nan
        max_block = df[block].max()
        minyear = df.index.min()
        ny = df.index.max() - minyear + 1 #XXX check the +1

        S = 0 # MK S score
        S_cv = 0 # MK S score of covariate
        var_S = 0 # variance of MK S Score
        var_S_cv = 0

        pc = np.nan # not sure
        var_Sc = np.nan # not sure

        partial_S = np.nan
        partial_S_var = np.nan
        partial_p = np.nan

        # Perform SKT/RKT on y
        for b in range(max_block):
            block_length = len(df[df[block] == b])
            if block_length > 3:
                block_S = df[block]
                S = S + block_S.apply(cls._s, axis=1)

                var_S = var_S + block_S.apply(cls._s_variance, axis=1)

        # this probably doesn't work XXX
        z = cls._z(S, var_S)
        p_value = 2*(1-norm.cdf(abs(z)))  # two tail test

        S_cv = S[cv]
        var_S_cv = var_S[cv]



        # for each block calculate sen and S

        # get a series with number of nans
        # get a series with number of samples

        # calculate sumties

        # calculate varS

        # If all blocks too short, throw an error





    def _skt_with_cv(self, df, block, cv):
        self.S = 0
        # eliminate records where covariable is missing
        df = df.dropna(subset=[cv])
        
        for b in range(self.maxblock):
            #calculate Scv for block
            self.S = df[df[block]==b]

        # for each block

        pass

    def _correct_for_intrablock_correlation():
        pass

    def _get_maxblock():
        pass



def rkt(y, block=None, cv=None, correct=False, rep="e"):
    """Regional Mann-Kendal Test

    Parameters
    ----------
    df : DataFrame
        Dateframe with datetime index and columns for each parameter to test

    block : str
        A column in df containing positive integers representing blocks, i.e. sites,
        seasons, or months, or a combination.

    cv : str
        Name of covariable column in df

    correct : boolean
        If False, no correction for correlation between blocks is performed.

    rep : str
        TODO

    Attribution
    -----------
    Aldo Marchetto's rkt package
    """
    # Check that df is correctly formated

    # if correcting for a covariable, eliminate records where missing
    df = df.dropna(subset=[cv])

    # if no block given create a dummy block, could do the same for cv

    # replace values sharing the same date with median or average

    # don't elimate records with missing values for y, rather use a threshold
    # below which to eliminate/ignore records with too few values

    # make blocks consecutive integers from 1 to maxblock

    # for each block calculate sen and S

    # get a series with number of nans
    # get a series with number of samples

    # calculate sumties

    # calculate varS

    # If all blocks too short, throw an error