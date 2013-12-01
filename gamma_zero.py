#!/usr/bin/python2.7
"""Runs script for error estimation"""

from ROOT import *

def read_data(filename1, filename2 ):
    current_data = []
    with open(filename1, 'r') as file, open(filename2, 'r') as file2:
        for line in file:
            current_data.append( [ float (l) for l in line.split() ] )
        for i, line in enumerate(file2):
            current_data[i].extend([ float(line.split()[0]) ])
    return current_data


def getGraph(lst):
    graph = TGraphErrors()
    graph.SetName('graph')
    graph.SetTitle('#Gamma(0,s)')

    [graph.SetPoint(i, p[0], p[1]) for i, p in enumerate(lst)]
    [graph.SetPointError(i, 0, p[2]) for i, p in enumerate(lst)]
    graph.SetMarkerStyle(20)
    graph.SetMarkerColor(38)

    graph.GetXaxis().SetTitle('#sqrt{s}, GeV')
    graph.GetYaxis().SetTitle('#Gamma(0)')

    return graph


canvas = TCanvas('canvas', 'Impact analysis', 800, 900)
canvas.cd()

gamma = read_data('gamma_at_zero.txt', 'gamma_at_zero_errors.txt')
graph = getGraph(gamma)
graph.Draw('AP')
canvas.Update()

raw_input('pease enter any key ...')
