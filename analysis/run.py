#!/usr/bin/python
import json
import unittest

from impact.datapoint import DataSet
from impact.impactanalysis import ImpactAnalysis
from impact.parametrization.symbolic import Symbolic
from impact.parametrization.symbolic import TripleExponent
from impact.parametrization.symbolic import TripleExponentGeneral
from impact.vis import Plots


class RunTheSolutin(unittest.TestCase):

    # @unittest.skip('')
    def test_visualize_the_results(self):
        with open('config/input.json') as f:
            data = json.load(f)

        models = {
            'config/datafit.json': Symbolic,
            'config/triple-exponent.json': TripleExponent,
            # 'config/triple-exponent-general.json': TripleExponentGeneral,
        }

        # TODO: Move conigs to the amplitude definitions
        #
        for config, algo in models.iteritems():
            dataset = DataSet(data['data'][-1])
            print algo.name
            visualisator = Plots()
            visualisator.draw_results(algo(), dataset, config)

            # analysis = ImpactAnalysis(algo(), config)
            # pandas = analysis.run(dataset)
            # pandas.to_csv(algo.name + ".csv")
            # visualisator = Plots()
            # visualisator.draw_results(model,
            # totem7tev, 'config/triple-exponent.json')
            # visualisator.draw_results(algo, dataset, config)
            # print pandas
