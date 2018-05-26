#!/usr/bin/python

from impactanalysis import ImpactAnalysis
from datapoint import DataSet
import json
import sys


def main():
    with open(sys.argv[1]) as f:
        data = json.load(f)

    for p in data['data']:
        analysis = ImpactAnalysis()
        dataset = DataSet(p)
        values, errors = analysis.run(dataset)


if __name__ == '__main__':
    main()
