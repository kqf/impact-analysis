#!/usr/bin/python2.7
"""Runs script for error estimation"""

from ROOT import *
from ComputeGamma import *
from Formulas import getRealGammaError, getRealGamma
import os




ENERGY = 7000
RHO    = 0.14
SIGMA  = 98.3
PROCESS = 'pp'

DSIGMA = 2.23
DRHO = 0.007


MC_AMOUNT = 100
N_POINTS = 30
MAX_B = 3



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
    [graph.SetPoint(i, i*3./N_POINTS, p[0]) for i, p in enumerate(lst)]
    [graph.SetPointError(i, 0, p[1]) for i, p in enumerate(lst)]
    return graph

def imGammaPoints(comp):
    im_gamma = []

    p = comp.parameters
    pe = comp.parametersErrors
    covariance = comp.covariance
    t_max = comp.t_max()

    for i in range(N_POINTS + 1):
        if i == 0:
            b = 1e-05
        else:
            b = 3.*i/N_POINTS
        val = getRealGamma([b], p)
        im_gamma.append([val , getRealGammaError([b], p,  covariance, DSIGMA, DRHO)])

    return im_gamma

# TODO: Rewrite  performComputationsMC return list
# Uncomment lines below to perform generation of files

mc = []
for i in range(MC_AMOUNT):
    c = ComputeGamma(PROCESS, ENERGY, SIGMA, RHO)
    mc.append( c.performComputationsMC(N_POINTS, i) )

    if i == range(MC_AMOUNT)[-1]:
        print 'go'
        os.remove(c.parametersFile)

    del c


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
    if idx == 1:
        canvas_new.cd()
        hist.GetXaxis().SetTitle('#Gamma(7000 GeV, 0.1)')
        hist.DrawClone()
        gaussFunction.DrawClone('same')
        canvas_new.SaveAs('/home/sha/study/12_semestr/science_seminar/1/pictures/8.eps')
    print 'Chi^2 = ', gaussFunction.GetChisquare(), idx

    mu = gaussFunction.GetParameter(1)
    sigma = gaussFunction.GetParameter(2)

    file.write(str(mu) + '\t' + str(sigma) + '\n')
    gamma_points.append([mu, sigma])
    del hist

my_gamma = ComputeGamma(PROCESS, ENERGY, SIGMA, RHO)
my_gamma.performComputations()
canvas_point = my_gamma.getCanvas()

im_gamma_points = imGammaPoints(my_gamma)

graph = getGraph(gamma_points)

canvas_point.cd(2)
graph.Draw('same')
canvas_point.Update()
canvas_point.SaveAs(str(ENERGY) + PROCESS + '.eps')

with open('main_result' + '.txt' , 'a') as file:
    for i in range(N_POINTS +  1):
        if len(gamma_points) != (N_POINTS + 1):
            print 'Your number of errors is ', len(gamma_points), 'instead of ', (N_POINTS + 1)

        b = MAX_B*i*1./N_POINTS if i != 0 else 1e-05
        re_gamma = my_gamma.getReal_Gamma([b])
        delta_re_gamma = gamma_points[i][1]
        im_gamma = im_gamma_points[i][0]
        delta_im_gamma = im_gamma_points[i][1]
        
        file.write('%f\t%f\t%f\t%f\t%f\t%f\n' % (
            b,
            re_gamma,
            delta_re_gamma,
            im_gamma,
            delta_im_gamma,
            ENERGY
            ))
    file.write('\n')

raw_input('pease enter any key ...')
