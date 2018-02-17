#!/usr/bin/python2
import ROOT
import json
import utils as ut

from impact.impactanalysis import ImpactAnalysis
from impact.datafit import DataFit

ROOT.TH1.AddDirectory(False)

class Plots(object):

    def __init__(self):
        super(Plots, self).__init__()
        self._cache = []

    def draw_results(self, model, dataset, conffile='config/datafit.json'):
        canvas = ut.canvas('The results')
        self._cache = []

        canvas.Divide(2, 1)
        self.fit(model, dataset, conffile, canvas.cd(1))
        self.draw_gamma(model, dataset, conffile, canvas.cd(2))
        
        canvas.Update()
        canvas.SaveAs('impact-analysis-{0}.eps'.format(dataset.energy))
        raw_input('Press enter to continue ...')



    def draw_gamma(self,
            model,
            dataset,
            conffile='config/datafit.json',
            canvas=ut.canvas('The results')
        ):

        ut.decorate_pad(canvas)
        analysis = ImpactAnalysis(model, conffile)
        output = analysis.run(dataset)

        # canvas.SetLogy(False)
        real_gamma = self.graph(output.index, -output['real_gamma'], output['real_gamma_error'], '-Re #Gamma(b)', 36)
        real_gamma.Draw('AP')

        print output['image_gamma']
        image_gamma = self.graph(output.index, output['image_gamma'], output['image_error'], 'Im #Gamma(b)', 46)
        image_gamma.Draw('same ap')

        self._cache.append(real_gamma)
        self._cache.append(image_gamma)
        canvas.Update()
        output.to_csv('impact-analysis-{0}.csv'.format(dataset.energy))



    def fit(self,
            model,
            dataset,
            conffile='config/datafit.json',
            canvas=ut.canvas("testfit")
        ):

        ut.decorate_pad(canvas)
        canvas.SetLogy(True)
        fitter = DataFit(model, conffile)
        fitter.fit(dataset)

        data = dataset.differential_cs()
        crossection = fitter.fitfunction(dataset.parameters)

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

        self._cache.append(data)
        self._cache.append(ratio)
        self._cache.append(crossection)
        self._cache.append(legend)

        print 'Fitted parameters:'
        print dataset.parameters[0:-2]
        canvas.Update()

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