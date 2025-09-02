''' 

Calulate 95% uncertainty intervals. Modified to calculate one-sided intervals if k=n or 0, 
to include the possibility that p=1 or 0.

'''
import numpy as np
import scipy.special
import matplotlib.pyplot as plt


def binomial_uncertainty(n, k):

    # Loop through probabilities from 0 to 1 with interval dp,
    dp = 0.001

    #Get start values:
    npoints = int(np.round(1 / dp) + 1)
    area_p = np.empty(npoints)
    true_p = np.empty(npoints)
    binomial_p = np.empty(npoints)

    area_p[0] = 0.0 #Cumulative area up to value p:
    true_p[0] = 0.0 #start value
    area_p_total = 0.0

    binomial_p[0] = scipy.special.binom(n, k)*(true_p[0]**k)*(1.0 - true_p[0])**(n-k)

    for i in range(1, npoints):
        true_p[i] = true_p[i-1] + dp
        # Clamp true_p[i] to avoid 0 or 1
        p = np.clip(true_p[i], 1e-10, 1 - 1e-10)

        # Compute log of binomial probability
        log_binom_coeff = np.log(scipy.special.binom(n, k))
        log_p_term = k * np.log(p)
        log_q_term = (n - k) * np.log(1.0 - p)
        log_prob = log_binom_coeff + log_p_term + log_q_term

        # Convert back to probability
        binomial_p[i] = np.exp(log_prob)

        this_area =  0.5 * dp * (binomial_p[i - 1] + binomial_p[i])
        area_p_total = area_p_total + this_area
        area_p[i] = area_p_total

    cdf = area_p / area_p_total

    if ((k!=0) and k!=n):
        p_low = np.interp(0.025,cdf,true_p)
        p_high = np.interp(0.975,cdf,true_p)
    else:
        if k == 0:
            p_low = 0
            p_high = np.interp(0.975, cdf, true_p)
        elif k==n: 
            p_low = np.interp(0.05, cdf, true_p)
            p_high = 1.0

    return p_low, p_high, true_p, binomial_p, cdf
