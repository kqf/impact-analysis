#!/usr/bin/python2.7
"""Runs script for error estimation"""

from ROOT import *
from ComputeGamma import *

# TODO: change name for function
def processPoint(filename, my_list):
    current_data = []
    with open(filename, 'r') as file:
        for line in file:
            current_data.append( float (line.split()[1]) )
    my_list.append( current_data )

def getHist(lst):
    pivot = lst[0]  # chosing some random pivot
    hist = TH1F('hist', 'gamma distribution', 100, pivot - 0.2, pivot - 0.2)
    [ hist.Fill(p) for  p in lst ]
    return hist

def getGraph(lst):
    graph = TGraphErrors()
    graph.SetName('graph')
    graph.SetTitle('graph')
    [graph.SetPoint(i, i*3./100, p[0]) for i, p in enumerate(lst)]
    [graph.SetPointError(i, 0, p[1]) for i, p in enumerate(lst)]
    return graph

ENERGY = 7000
RHO    = 0.141
SIGMA  = 98.58
PROCESS = 'pp'

MC_AMOUNT = 100




# TODO: Rewrite  performComputationsMC return list
# Uncomment lines below to perform generation of files
for i in range(MC_AMOUNT):
    c = ComputeGamma(PROCESS, ENERGY, SIGMA, RHO , True)
    c.performComputationsMC(100, i)
    del c

mc = []
for i in range(MC_AMOUNT):
    processPoint(PROCESS + str(ENERGY) + str(i) + '.dat', mc)

gaussFunction = TF1('gaussFunction','gaus', 0, 1.4)
gaussFunction.SetParameter(0, 1)
gaussFunction.SetParameter(1, 0.1)
gaussFunction.SetParameter(2, 1)
    

file = open('errors.dat','w')
gamma_points = []


canvas_new = TCanvas('canvas_new', 'Canvas', 800, 600)

for idx in range(len(mc[0])):
    """iterate over each point"""
    ds_point = [mc[j][idx] for j in range(len(mc))]
    hist = getHist(ds_point)
    hist.Fit(gaussFunction,'r')
    ndf = gaussFunction.GetNDF()
    if idx == 89:
        canvas_new.cd()
        hist.DrawClone()
        gaussFunction.DrawClone('same')
    print 'Chi^2 = ', gaussFunction.GetChisquare(), idx

    mu = gaussFunction.GetParameter(1)
    sigma = gaussFunction.GetParameter(2)

    file.write(str(mu) + '\t' + str(sigma) + '\n')
    gamma_points.append([mu, sigma])
    del hist

my_gamma = ComputeGamma(PROCESS, ENERGY, SIGMA, RHO)
my_gamma.performComputations()

canvas_point = my_gamma.getCanvas()
graph = getGraph(gamma_points)

canvas_point.cd(2)
graph.Draw('same')
canvas_point.Update()

raw_input('pease enter any key ...')
