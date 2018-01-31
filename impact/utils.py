import ROOT

def decorate_pad(pad):
    pad.SetLogy()
    pad.SetTickx()
    pad.SetTicky()
    pad.SetGridy()
    pad.SetGridx()

def canvas(name='name', x=6, y=6, scale=0.5):
    canvas = ROOT.TCanvas(name, 'Canvas', int(128 * x * scale) , int(96 * y * scale))
    return canvas
    # return adjust_canvas(canvas)