import ROOT
import pandas as pd

import impact.utils as ut


def prepare_dataset(name_color):
    name, color = name_color
    hist = ROOT.TH1F(
        name.replace("-", "_"),
        name.replace("-", " "),
        100, 0.49, 0.53)

    df = pd.read_csv("13000GeV/mc-gamma-{}.csv".format(name))
    data = df.values.T[0] / 2.
    print len(data)
    map(hist.Fill, data)
    hist.SetLineColor(color)
    return hist


def test_compare_mc_productions():
    models = [
        ("standard-parametrization", ROOT.kRed + 2),
        ("full-three-plus-one", ROOT.kGreen + 2),
        ("full-three-plus-two", ROOT.kBlue + 2),
    ]

    hists = map(prepare_dataset, models)
    legend = ROOT.TLegend(0.80, 0.7, 0.995, 0.8)
    stack = ROOT.THStack("stack", "Generated Im H(b=0); ImH(0); counts")
    canvas = ut.canvas()
    canvas.cd()
    for hist in hists:
        stack.Add(hist)
        legend.AddEntry(hist, hist.GetTitle())
    stack.Draw("nostack")
    legend.Draw("same")
    canvas.Update()
    raw_input("Press enter to continue")
