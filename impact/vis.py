#!/usr/bin/python2
import ROOT
import json
import utils as ut

from impact.impactanalysis import ImpactAnalysis
from impact.model import Approx
from impact.datafit import DataFit

class Plots(object):

    def fit(self, model, dataset, conffile='config/datafit.json'):
        fitter = DataFit(model, conffile)
        fitter.fit(dataset)

        data = dataset.differential_cs()
        crossection = fitter.fitfunction(dataset.parameters)

        canvas = ut.canvas("testfit")
        data.Draw("AP")
        crossection.Draw("same")
        crossection.SetLineColor(3)

        ratio = fitter.fitfunction(dataset.parameters, 'ratio')
        ratio.SetLineColor(2)
        ratio.Draw("same")

        legend = ROOT.TLegend(0.7, 0.6, 0.9, 0.9);
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.AddEntry(crossection,"d#sigma/dt param","l");
        legend.AddEntry(ratio, "#rho(t) with the same parameters","l");
        legend.Draw("same")

        print 'Fitted parameters:'
        print dataset.parameters[0:-2]

        canvas.Update()
        raw_input()

    def draw_results(self, model, dataset, conffile='config/datafit.json' ):
        analysis = ImpactAnalysis(model, conffile)
        output = analysis.run(dataset)

        canvas = ut.canvas('The results')
        # canvas.SetLogy(False)
        real_gamma = self.graph(output.index, -output['real_gamma'], output['real_gamma_error'], '-Re #Gamma(b)', 36)
        real_gamma.Draw('AP')

        imag_gamma = self.graph(output.index, output['image_gamma'], output['image_error'], 'Im #Gamma(b)', 46)
        imag_gamma.Draw('same ap')

        canvas.Update()
        raw_input()



    @classmethod
    def graph(cls, b, data, errors, title, color):
        graph = ROOT.TGraphErrors(len(data))
        graph.SetTitle(title)

        for i, (b, d, e) in enumerate(zip(b, data, errors)):
            graph.SetPoint(i, b, d)
            graph.SetPointError(i, 0, e)

        graph.GetXaxis().SetTitle('b, fm')
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(color)
        return graph