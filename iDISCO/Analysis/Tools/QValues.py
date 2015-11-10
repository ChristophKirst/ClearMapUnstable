# -*- coding: utf-8 -*-
"""
Q-Value estimation from p values

Created on Mon Nov  9 11:28:09 2015

Source:
    - https://github.com/nfusi/qvalue

Reference:
    - Storey and Tibshirani, 2003

Note:
    - The q-value of a particular feature can be described as the expected proportion of
      false  positives  among  all  features  as  or  more  extreme  than  the observed one
    
    - The estimated q-values are increasing in the same order as the p-values


modified by: ckirst
"""

import scipy

def estimateQValues(pvalues, m = None, pi0 = None, verbose = False, lowMemory = False):
    """Estimates q-values from p-values

    m: number of tests. if None, m = pvalues.size
    pi0: is estimate of m_0 / m which is the (true null / total tests) ratio, if None estimation via cubic spline as in Storey and Tibshirani, 2003.
    """

    if not (pvalues.min() >= 0 and pvalues.max() <= 1):
        raise RuntimeError("estimateQValues: p-values should be between 0 and 1");

    original_shape = pvalues.shape
    pvalues = pvalues.ravel() # flattens the array in place, more efficient than flatten() 

    if m == None:
        m = float(len(pvalues))
    else:
        # the user has supplied an m
        m *= 1.0

    # if the number of hypotheses is small, just set pi0 to 1
    if len(pvalues) < 100 and pi0 == None:
        pi0 = 1.0
    elif pi0 != None:
        pi0 = pi0
    else:
        # evaluate pi0 for different lambdas
        pi0 = []
        lam = scipy.arange(0, 0.90, 0.01)
        counts = scipy.array([(pvalues > i).sum() for i in lam])
        
        for l in range(len(lam)):
            pi0.append(counts[l]/(m*(1-lam[l])))

        pi0 = scipy.array(pi0)

        # fit natural cubic scipyline
        tck = scipy.interpolate.scipylrep(lam, pi0, k = 3)
        pi0 = scipy.interpolate.scipylev(lam[-1], tck)
        
        if pi0 > 1:
            if verbose:
                raise Warning("estimateQValues: got pi0 > 1 (%.3f) while estimating qvalues, setting it to 1" % pi0);
            pi0 = 1.0
        
    if not (pi0 >= 0 and pi0 <= 1):
        raise RuntimeError("estimateQValues: pi0 is not between 0 and 1: %f" % pi0);


    if lowMemory:
        # low memory version, only uses 1 pvalues and 1 qv matrices
        qv = scipy.zeros((len(pvalues),))
        last_pvalues = pvalues.argmax()
        qv[last_pvalues] = (pi0*pvalues[last_pvalues]*m)/float(m)
        pvalues[last_pvalues] = -scipy.inf
        prev_qv = last_pvalues
        for i in xrange(int(len(pvalues))-2, -1, -1):
            cur_max = pvalues.argmax()
            qv_i = (pi0*m*pvalues[cur_max]/float(i+1))
            pvalues[cur_max] = -scipy.inf
            qv_i1 = prev_qv
            qv[cur_max] = min(qv_i, qv_i1)
            prev_qv = qv[cur_max]

    else:
        p_ordered = scipy.argsort(pvalues)    
        pvalues = pvalues[p_ordered]
        qv = pi0 * m/len(pvalues) * pvalues
        qv[-1] = min(qv[-1],1.0)

        for i in xrange(len(pvalues)-2, -1, -1):
            qv[i] = min(pi0*m*pvalues[i]/(i+1.0), qv[i+1])
        
        # reorder qvalues
        qv_temp = qv.copy()
        qv = scipy.zeros_like(qv)
        qv[p_ordered] = qv_temp

        # reshape qvalues
        qv = qv.reshape(original_shape)
        
    return qv