#!/usr/bin/python

from impact.impactanalysis import ImpactAnalysis
import json
import sys

def main():
	with open(sys.argv[1]) as f:
		data = json.load(f)

	infile = data['infile']
	for p in data['data']:
		analysis = ImpactAnalysis(infile, p["PROCESS"], p["ENERGY"], p["SIGMA"], p["RHO"], p["DSIGMA"], p["DRHO"], 100)
		values, errors = analysis.run()

if __name__ == '__main__':
	main()