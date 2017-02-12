# from AlleleSeq_v1.2a; no modification

import math
import scipy.stats

def binomtest(x, n, p):
    #return (scipy.stats.binom_test(x, n, p), normal_approx(x, n, p))
    if n*p > 50:
        return normal_approx(x, n, p)
    else:
        return scipy.stats.binom_test(x, n, p)


def normal_approx(x, n, p):
    if abs(x-n*p) < 1e-5:
        return 1.0
    u=p*n
    s=math.sqrt(n*p*(1-p))
    norm=scipy.stats.distributions.norm(u,s)
    if x<n*p:
        pval=2*norm.cdf(x+.5) # add 0.5 for continuity correction
    else:
        pval=2*(1-norm.cdf(x-.5)) 
    return pval


    
    
