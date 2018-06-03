from math import pi, sqrt

import ROOT
from constants import k_fm, k_norm
from scipy import integrate
from scipy.special import j0


def decorate_pad(pad):
    pad.SetLogy()
    pad.SetTickx()
    pad.SetTicky()
    pad.SetGridy()
    pad.SetGridx()
    return pad


def canvas(name='name', title='canvas', x=5, y=5, scale=1.0):
    canvas = ROOT.TCanvas(name, title, int(
        128 * x * scale), int(96 * y * scale))
    return decorate_pad(canvas)
    # return adjust_canvas(canvas)


def hankel_transform(func):
    def impact_version(b, p, limits=(0, float("inf"))):
        def f(q):
            return q * j0(b * q / k_fm) * func(q * q, p) / 8 / pi
        # integral from zero to lower bound
        result = integrate.quad(f, *limits)[0]
        return result

    return impact_version


def impact_range(npoints=30, step=10.0, zero=1e-5):
    return (zero * (i == 0) + i / step for i in range(npoints))
