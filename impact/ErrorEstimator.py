#!/usr/bin/python2.7
"""Runs script for error estimation"""

from ROOT import *
from ComputeGamma import *
from Formulas import getRealGammaError, getRealGamma
import os
import progressbar

# TODO: change the name of the function
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
    [graph.SetPoint(i, i * 3./100, p[0]) for i, p in enumerate(lst)]
    [graph.SetPointError(i, 0, p[1]) for i, p in enumerate(lst)]
    return graph



class ErrorEstimator(object):
    def __init__(self, ptype, energy, sigma, rho, dsigma, drho, nmc):
        super(ErrorEstimator, self).__init__()
        self.process = ptype
        self.energy = energy
        self.sigma = sigma
        self.rho = rho  
        self.dsigma = dsigma
        self.drho = drho
        self.nmc = nmc

        self.gausf = TF1('gaussFunction','gaus', 0, 1.4)
        self.gausf.SetParameter(0, 1)
        self.gausf.SetParameter(1, 0.1)
        self.gausf.SetParameter(2, 1)
 
    def generate_mc_data(self):
        mc = []
        bar = progressbar.ProgressBar()
        for i in bar(range(self.nmc)):
            c = ComputeGamma(self.process, self.energy, self.sigma, self.rho)
            mc.append( c.performComputationsMC(100, i, self.dsigma) )
            filename = c.gamma_fitter.par_file_name
            del c 
        os.remove(filename)
        return mc

    def main(self):
        # TODO: Rewrite  performComputationsMC return list
        # Uncomment lines below to perform generation of files

        mc = self.generate_mc_data()
        file = open('errors.dat','w')

        # Move average values to the separate function
        gamma_points = []
        ROOT.gROOT.SetBatch(True)
        canvas_new = TCanvas('canvas_new', 'Canvas', 800, 600)
        for idx in range(len(mc[0])):
            """iterate over each point"""
            ds_point = [mc[j][idx] for j in range(len(mc))]
            hist = getHist(ds_point)
            hist.Fit(self.gausf,'0qr')
            ndf = self.gausf.GetNDF()
            # if idx == 89:
                # canvas_new.cd()
                # hist.DrawClone()
                # self.gausf.DrawClone('same')
            print 'Chi^2 = ', self.gausf.GetChisquare(), idx

            mu = self.gausf.GetParameter(1)
            sigma = self.gausf.GetParameter(2)

            file.write(str(mu) + '\t' + str(sigma) + '\n')
            gamma_points.append([mu, sigma])
            del hist
        ROOT.gROOT.SetBatch(False)

        my_gamma = ComputeGamma(self.process, self.energy, self.sigma, self.rho )
        my_gamma.performComputations()

        canvas_point = my_gamma.gamma_fitter.canvas
        graph = getGraph(gamma_points)

        canvas_point.cd(2)
        graph.Draw('same')
        canvas_point.Update()
        canvas_point.SaveAs(str(self.energy) + self.process + '.eps')
        # raw_input('pease enter any key ...')
        return zip(*gamma_points)

        # TODO: Fix this
        with open('gamma_at_zero_errors.txt', 'a') as file:
            file.write('%f\t%f\t%f\t%f\t%f\n' % 
                                                (
                                                   my_gamma.gammaAtZero ,
                                                   gamma_points[0][1],
                                                   getRealGamma([0e-5], my_gamma.parameters),
                                                   getRealGammaError([0e-5], my_gamma.parameters, my_gamma.covariance, self.dsigma, self.drho),
                                                  self.energy 
                                                )
                      )

        raw_input('pease enter any key ...')

if __name__ == '__main__':
    main()