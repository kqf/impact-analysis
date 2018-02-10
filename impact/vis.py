#!/usr/bin/python2
import ROOT
import json
import utils as ut

from impact.model import Approx
from impact.datafit import DataFit

class Plots(object):

    def fit(self, model, dataset, conffile='config/datafit.json' ):
    	fitter = DataFit(model, conffile)
    	fitter.fit(dataset)

    	data = dataset.differential_cs()
    	approximation = fitter.fitfunction(dataset.parameters)

    	canvas = ut.canvas("testfit")
    	data.Draw()
    	approximation.Draw("same")

    	canvas.Update()
    	print 'Fitted parameters:'
    	print dataset.parameters[0:-2]
    	raw_input()
