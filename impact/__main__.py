#!/usr/bin/python

from impactanalysis import ImpactAnalysis
import json
import sys

def main():
    with open(sys.argv[1]) as f:
        data = json.load(f)

    infile = data['infile']
    for p in data['data']:
        # TODO: Wrap this up in options
        analysis = ImpactAnalysis(infile, p)
        values, errors = analysis.run()

if __name__ == '__main__':
    main()