import pandas as pd

from trend import mk_score, mk_score_variance

class RKT:
    def test(self, df, block=None, cv=None, correct):
        pass

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



def rkt(df, block=None, cv=None, correct=False, rep="e"):
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