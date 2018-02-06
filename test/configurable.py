import unittest
import json
import hashlib

from impact.impactanalysis import ImpactAnalysis
from impact.datapoint import DataSet

class Configurable(unittest.TestCase):

    def setUp(self):
        with open('config/test.json') as f:
            self.data = json.load(f)

        p = lambda x: self.data[x]
        self.infile = p('infile')
        self.dataset = DataSet(self.infile, self.data)

        hsum = self.data['hsum256']
        self.checkInputFile(hsum)
        self.longMessage = True


    def checkInputFile(self, hsum):
        hashsum = hashlib.sha256()
        with open(self.infile, 'rb') as f:
                data = f.read()
        hashsum.update(data)

        msg = "You are using a wrong file to test your data. Your hash sums don't coincide."\
              "\n\nActual:  {}\nNominal: {}".format(hashsum.hexdigest(), hsum)
        self.assertEqual(hashsum.hexdigest(), hsum, msg)
