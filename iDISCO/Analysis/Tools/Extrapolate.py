# -*- coding: utf-8 -*-
"""
Extrapolate

- extends interp1d to constantly / linearly extrapolate.

@author: ckirst
"""


import scipy.interpolate

from numpy import array


def extrap1d(x, y, interpolation = 'linear', exterpolation = 'constant'):
    interpolator = scipy.interpolate.interp1d(x, y, kind = interpolation);
    return extrap1dFromInterp1d(interpolator, exterpolation);   


def  extrap1dFromInterp1d(interpolator, exterpolation = 'constant'):
    xs = interpolator.x
    ys = interpolator.y
    cs = (exterpolation == 'constant');

    def pointwise(x):
        if cs:   #constant extrapolation
            if x < xs[0]:
                return ys[0];
            elif x > xs[-1]:
                return ys[-1];
            else:
                return interpolator(x);
        else:  # linear extrapolation
            if x < xs[0]:
                return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
            elif x > xs[-1]:
                return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
            else:
                return interpolator(x)

    def extrapfunc(xs):
        return array(map(pointwise, array(xs)))

    return extrapfunc