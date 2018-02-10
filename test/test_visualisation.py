import unittest

from test.configurable import Configurable
from impact.parametrization.numeric import Numeric
from impact.vis import Plots

class TestVisualRepresentation(Configurable):

    def testValues(self):
        visualisator = Plots()
        visualisator.fit(Numeric(), self.dataset)   
