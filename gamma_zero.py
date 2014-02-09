#!/usr/bin/python2.7
"""Runs script for error estimation"""

from ROOT import *

def read_data(filename):
    current_data = []
    with open(filename, 'r') as file:
        for line in file:
            current_data.append( [ float (l) for l in line.split() ] )
    return current_data


def getGraph(lst):
    graph = TGraphErrors()
    graph.SetName('graph')
    graph.SetTitle('#Gamma(0,s)')

    [graph.SetPoint(i, p[2], p[0]) for i, p in enumerate(lst)]
    [graph.SetPointError(i, 0, p[1]) for i, p in enumerate(lst)]
    graph.SetMarkerStyle(20)
    graph.SetMarkerColor(38)

    graph.GetXaxis().SetTitle('#sqrt{s}, GeV')
    graph.GetYaxis().SetTitle('#Gamma(0)')

    return graph


canvas = TCanvas('canvas', 'Impact analysis', 800, 900)
canvas.cd()

gamma = read_data('gamma_at_zero_errors.txt')
graph = getGraph(gamma)
graph.Draw('AP')
canvas.SetLogx()
canvas.Update()

raw_input('pease enter any key ...')
