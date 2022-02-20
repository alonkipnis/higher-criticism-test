import numpy as np


class HCtest(object):
    """
    Higher Criticism test

    References:
    [1] Donoho, D. L. and Jin, J.,
     "Higher criticism for detecting sparse hetrogenous mixtures", 
     Annals of Stat. 2004
    [2] Donoho, D. L. and Jin, J. "Higher critcism thresholding: Optimal 
    feature selection when useful features are rare and weak", proceedings
    of the national academy of sciences, 2008.
    ========================================================================

    Args:
    -----
        pvals    list of p-values. P-values that are np.nan are exluded.
        stbl     normalize by expected P-values (stbl=True) or observed 
                 P-values (stbl=False). stbl=True was suggested in [2].
                 stbl=False in [1]. 
        gamma    lower fruction of p-values to use.
        
    Methods :
    -------
        HC       HC and P-value attaining it
        HCstar   sample adjustet HC (HCdagger in [1])
        HCjin    a version of HC from 
                [2] Jiashun Jin and Wanjie Wang, "Influential features PCA for
                 high dimensional clustering"
        
    """

    def __init__(self, pvals, stbl=True):
        self._N = len(pvals)
        assert (self._N > 0)

        self._EPS = 1 / (1e2 + self._N ** 2)
        self._istar = 1

        self._pvals = np.sort(np.asarray(pvals.copy()))
        self._uu = np.linspace(1 / self._N, 1 - self._EPS, self._N)

        if stbl:
            denom = np.sqrt(self._uu * (1 - self._uu))
        else:
            denom = np.sqrt(self._pvals * (1 - self._pvals))

        denom = np.maximum(denom, self._EPS)
        self._zz = np.sqrt(self._N) * (self._uu - self._pvals) / denom

        self._imin_star = np.argmax(self._pvals > (1 - self._EPS) / self._N)
        self._imin_jin = np.argmax(self._pvals > np.log(self._N) / self._N)

    def __call__(self, gamma=0.2):
        return self.HC(gamma=gamma)

    def _compute_hc(self, imin, imax):
        if imin > imax:
            return np.nan
        if imin == imax:
            self._istar = imin
        else:
            self._istar = np.argmax(self._zz[imin:imax]) + imin
        zMaxStar = self._zz[self._istar]
        return zMaxStar, self._pvals[self._istar]

    def HC(self, gamma=0.2):
        """
        Higher Criticism test statistic

        Args:
        -----
        'gamma' : lower fraction of P-values to consider

        Return:
        -------
        HC test score, P-value attaining it

        """
        imin = 0
        imax = np.maximum(imin, int(gamma * self._N + 0.5))
        return self._compute_hc(imin, imax)

    def HCjin(self, gamma=0.2):
        """sample-adjusted higher criticism score from [2]

        Args:
        -----
        'gamma' : lower fraction of P-values to consider

        Return:
        -------
        HC score, P-value attaining it

        """

        imin = self._imin_jin
        imax = np.maximum(imin + 1, int(gamma * self._N + 0.5))
        return self._compute_hc(imin, imax)

    def HCstar(self, gamma=0.2):
        """sample-adjusted higher criticism score

        Args:
        -----
        'gamma' : lower fraction of P-values to consider

        Returns:
        -------
        HC score, P-value attaining it

        """

        imin = self._imin_star
        imax = np.maximum(imin + 1, int(gamma * self._N + 0.5))
        return self._compute_hc(imin, imax)

    def get_state(self):
        return {'pvals': self._pvals,
                'u': self._uu,
                'z': self._zz,
                'imin_star': self._imin_star,
                'imin_jin': self._imin_jin,
                }
