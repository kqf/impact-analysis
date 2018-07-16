#!/usr/bin/python2
import ROOT
import utils as ut

from impact.impactanalysis import ImpactAnalysis
from impact.datafit import DataFit

ROOT.TH1.AddDirectory(False)


class Plots(object):

    def __init__(self):
        super(Plots, self).__init__()
        self._cache = []

    def draw_results(self, model, dataset, conffile="config/standard.json"):
        canvas = ut.canvas("The results {0}".format(model.name), "Results")
        canvas.SetMargin(0, 0, 0, 0)
        self._cache = []

        canvas.Divide(2, 1, 0.01, 0)
        self.fit(model, dataset, canvas.cd(1), conffile)
        # self.draw_gamma(model, dataset, canvas.cd(2), conffile)

        canvas.Update()
        canvas.SaveAs(
            "impact-analysis-{0}-{1}.eps".format(dataset.energy, model.name))
        raw_input("Press enter to continue ...")

    def draw_gamma(self,
                   model,
                   dataset,
                   canvas,
                   conffile="config/standard.json",
                   ):

        ut.decorate_pad(canvas)
        canvas.SetMargin(0.15, 0.01, 0.09, 0.01)
        analysis = ImpactAnalysis(model, conffile)
        output = analysis.run(dataset)

        canvas.SetLogy(False)
        real_gamma = self.graph(
            output.index, -output["real_gamma"],
            output["real_gamma_error"], "-Re #Gamma(b)", 36)
        real_gamma.Draw("AP")

        # print output["imag_gamma"]
        imag_gamma = self.graph(
            output.index, output["imag_gamma"],
            output["imag_gamma_error"], "Im #Gamma(b)", 46)
        imag_gamma.Draw("same")

        g_inel = self.graph(
            output.index, output["g_inel"],
            output["g_inel_error"], "Im #Gamma(b)", 4)
        g_inel.Draw("same")

        multigraph = ROOT.TMultiGraph()
        multigraph.SetTitle(
            "; b (fm); H(s,b), G_{inel}(s, b)")  # H(s, b), G_{inel}(s, b)
        multigraph.Add(real_gamma, "cp")
        multigraph.Add(imag_gamma, "cp")
        multigraph.Add(g_inel, "cp")
        multigraph.Draw("ap")
        legend = ROOT.TLegend(0.6, 0.75, 0.9, 0.95)
        legend.SetBorderSize(0)
        legend.SetFillStyle(1001)
        legend.SetFillColor(0)
        legend.AddEntry(imag_gamma, "Im #Gamma(s, b)", "lep")
        legend.AddEntry(real_gamma, "Re #Gamma(s, b)", "lep")
        legend.AddEntry(g_inel, "G_{inel}(s, b)", "lep")
        legend.Draw("same")
        self._cache.append(legend)
        self._cache.append(multigraph)
        self._cache.append(real_gamma)
        self._cache.append(imag_gamma)
        self._cache.append(g_inel)

        canvas.Update()
        output = output[["im_h", "re_h", "im_h_error",
                         "re_h_error", "g_inel", "g_inel_error"]]
        ofile = "{1}-{0}GeV.csv".format(
            dataset.energy, model.name)

        print "saving the file", ofile
        output.to_csv(ofile, sep="\t")
        return model

    def fit(self,
            model,
            dataset,
            canvas,
            conffile="config/standard.json"
            ):

        canvas.Divide(1, 2, 0, 0)
        ut.decorate_pad(canvas.cd(1))
        canvas.cd(1).SetMargin(0.15, 0.01, 0, 0.01)
        canvas.cd(1).SetLogy(True)
        fitter = DataFit(model, conffile)
        fitter.fit(dataset)

        data = dataset.differential_cs()
        data.SetTitle()
        data.GetYaxis().SetTitleSize(1.5 * data.GetYaxis().GetTitleSize())
        data.GetYaxis().SetTitle("d#sigma/dt (mb/GeV^{2})")
        crossection = fitter.fitfunction(dataset.parameters)

        data.Draw("AP")
        crossection.Draw("same")
        crossection.SetLineColor(3)

        legend = ROOT.TLegend(0.3, 0.82, 0.5, 0.9)
        legend.SetBorderSize(0)
        legend.SetFillStyle(1001)
        legend.SetFillColor(0)
        legend.AddEntry(crossection, "d#sigma/dt fit", "l")
        legend.Draw("same")

        ut.decorate_pad(canvas.cd(2))
        canvas.cd(2).SetMargin(0.15, 0.01, 0.18, 0)
        canvas.cd(2).SetLogy(False)

        ratio = fitter.fitfunction(dataset.parameters, "ratio")
        ratio.SetLineColor(2)

        hax = ROOT.TH1D("hax-{0}".format(model.name), "", 10,
                        data.GetXaxis().GetXmin(), data.GetXaxis().GetXmax())
        hax.GetYaxis().SetTitleSize(1.5 * hax.GetYaxis().GetTitleSize())
        hax.GetXaxis().SetTitleSize(1.5 * hax.GetXaxis().GetTitleSize())
        hax.GetXaxis().SetLabelSize(1.5 * hax.GetXaxis().GetLabelSize())
        hax.GetXaxis().SetTitle("-t (GeV^{2}/c^{2})")
        hax.GetYaxis().SetTitle("#rho")
        # ratio.GetMinimum(data.GetXaxis().GetXmin() + 1,data.GetXaxis().GetXmax()) - 0.52,ratio.GetMaximum(data.GetXaxis().GetXmin() + 1,data.GetXaxis().GetXmax()) + 0.52
        hax.SetAxisRange(-5.23, 5.23, "Y")
        hax.SetStats(False)
        hax.Draw()

        ratio.Draw("same")

        legendr = ROOT.TLegend(0.3, 0.82, 0.5, 0.9)
        legendr.SetBorderSize(0)
        legendr.SetFillStyle(1001)
        legendr.SetFillColor(0)
        legendr.AddEntry(ratio, "#rho(t)", "l")
        legendr.Draw("same")

        self._cache.append(data)
        self._cache.append(ratio)
        self._cache.append(crossection)
        self._cache.append(legend)
        self._cache.append(hax)
        self._cache.append(legendr)

        print "Fitted parameters:"
        print dataset.parameters[0:-2]
        canvas.Update()

    @classmethod
    def graph(cls, b, data, errors, title, color):
        graph = ROOT.TGraphErrors(len(data))
        graph.SetTitle(title)

        for i, (b, d, e) in enumerate(zip(b, data, errors)):
            graph.SetPoint(i, b, d)
            graph.SetPointError(i, 0, e)

        graph.GetXaxis().SetTitle("b, fm")
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(color)
        return graph
