import ROOT
from scipy import integrate
from scipy.special import j0, j1
from constants import k_fm, k_norm
from math import sqrt, pi


def decorate_pad(pad):
    pad.SetLogy()
    pad.SetTickx()
    pad.SetTicky()
    pad.SetGridy()
    pad.SetGridx()

def canvas(name='name', x=6, y=6, scale=0.5):
    canvas = ROOT.TCanvas(name, 'Canvas', int(128 * x * scale) , int(96 * y * scale))
    return canvas
    # return adjust_canvas(canvas)

def hankel_transform(func):
    def impact_version(b, p, limits = (0, float("inf"))):
        f = lambda q : q * j0(b * q / k_fm) *  func(q * q, p) / sqrt(pi * k_norm)
        result = integrate.quad(f, *limits)[0]  # integral from zero to lower bound
        return result

    return impact_version 

def impact_range(npoints = 30, step = 10.0, zero = 1e-5):
    return (zero * (i == 0) + i / step for i in range(npoints))