
import ROOT
import progressbar

class ImageError(object):
    def __init__(self, gamma_estimator, nmc, dsigma):
        super(ImageError, self).__init__()
        self.nmc = nmc 
        self.dsigma = dsigma
        self.gamma_estimator = gamma_estimator 
        self.gausf = ROOT.TF1('gaussFunction','gaus', 0, 1.4)
        self.gausf.SetParameter(0, 1)
        self.gausf.SetParameter(1, 0.1)
        self.gausf.SetParameter(2, 1)       

    def average_and_deviation(self, data):
        pivot = data[0]  # chosing some random pivot
        hist = ROOT.TH1F('hist', 'gamma distribution', 100, pivot - 0.2, pivot - 0.2)
        map(hist.Fill, data)
        hist.Fit(self.gausf,'0qr')
        mu, sigma = self.gausf.GetParameter(1), self.gausf.GetParameter(2)
        return mu, sigma

    def generate_mc_data(self):
        bar = progressbar.ProgressBar()
        print 'Generating the sample of gamma points'
        mc = [self.gamma_estimator.generate_mc_gamma(100, i, self.dsigma) for i in bar(range(self.nmc))]
        return mc

    def estimate_deviations(self, mc):
        print 'Calculating mean, and deviation for all gamma values'
        bar = progressbar.ProgressBar()
        return [self.average_and_deviation(p) for p in bar(zip(*mc))]

    def evaluate(self):
        # Generate monte-carlo data
        mc = self.generate_mc_data()
        mc_av_and_deviation = self.estimate_deviations(mc)
        return mc_av_and_deviation