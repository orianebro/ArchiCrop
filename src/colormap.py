# -*- python -*-
#
#       Copyright 2015 INRIA - CIRAD - INRA
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       WebSite : https://github.com/openalea-incubator/caribu
#
# ==============================================================================
""" A colormap class
    copied from  openalea.color.colormap because
    import openalea.color.colormap makes matplotlib
    with TkAgg frontend crash for some mysterious reason !!
"""

from math import isnan

class ColorMap(object):
    """A RGB color map, between 2 colors defined in HSV code

    :Examples:

    >>> minh,maxh = minandmax([height(i) for i in s2])
    >>> colormap = ColorMap(minh,maxh)
    >>> s3 = [ Shape(i.geometry, Material
    >>>    (Color3(colormap(height(i))), 1), i.id)
    >>>    for i in s2]

    """

    def __init__(self, minval=0., maxval=1.):
        self.minval = float(minval)
        self.maxval = float(maxval)

    def color(self, normedU):
        """

        :param normedU: todo

        """
        inter = 1 / 5.
        winter = int(normedU / inter)
        a = (normedU % inter) / inter
        b = 1 - a

        if winter < 0:
            col = (self.coul2, self.coul2, self.coul1)
        elif winter == 0:
            col = (self.coul2, self.coul2 * b + self.coul1 * a, self.coul1)
        elif winter == 1:
            col = (self.coul2, self.coul1, self.coul1 * b + self.coul2 * a)
        elif winter == 2:
            col = (self.coul2 * b + self.coul1 * a, self.coul1, self.coul2)
        elif winter == 3:
            col = (self.coul1, self.coul1 * b + self.coul2 * a, self.coul2)
        elif winter > 3:
            col = (self.coul1, self.coul2, self.coul2)
        return (int(col[0]), int(col[1]), int(col[2]))

    def greycolor(self, normedU):
        """

        :param normedU: todo
        :returns: todo
        """
        return (int(255 * normedU), int(255 * normedU), int(255 * normedU))

    def grey(self, u):
        """
        :param u:
        :returns: todo
        """
        return self.greycolor(self.normU(u))

    def normU(self, u):
        """
        :param u:
        :returns: todo
        """
        if self.minval == self.maxval:
            return 0.5
        return (u - self.minval) / (self.maxval - self.minval)

    def __call__(self, u, minval=0, maxval=1, coul1=80, coul2=20):
        self.coul1 = coul1
        self.coul2 = coul2
        self.minval = float(minval)
        self.maxval = float(maxval)
        return self.color(self.normU(u))


def nan_to_zero(values):
    return [0 if isnan(x) else x for x in values]


def jet_colors(values, minval=None, maxval=None):
    """return jet colors associated to values after gamme normalisation

    Args:
        values: (list of float) input values
        minval: (float) minimal value at lower bound of color range
        maxval: (float) maximal value at upper bound of color range

    Returns:
        a list of (r, g, b) tuples
    """

    values = nan_to_zero(values)
    if minval is None:
        minval = min(values)
    if maxval is None:
        maxval = max(values)
    cmap = ColorMap()
    return map(lambda x: cmap(x, minval, maxval, 250., 20.), values)
