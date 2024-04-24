import numpy as np

def geometric_dist(height, nb_phy, q=1):
    """ returns distances between individual leaves along a geometric model """

    if q==1:
        u0=float(height)/nb_phy
    else:
        u0=height*(1.-q)/(1.-q**(nb_phy+1))

    return [u0*q**i for i in range(nb_phy)]


def bell_shaped_dist(max_leaf_length, nb_phy, rmax=.8, skew=0.15):
    """ returns leaf area of individual leaves along bell shaped model """

    k = -np.log(skew) * rmax
    r = np.linspace(1. / nb_phy, 1, nb_phy)
    relative_length = np.exp(-k / rmax * (2 * (r - rmax) ** 2 + (r - rmax) ** 3))
    # leaf_length = relative_length / relative_length.sum() * max_leaf_length
    leaf_length = relative_length * max_leaf_length
    return leaf_length.tolist()