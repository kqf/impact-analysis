#!/usr/bin/python

from impact.ErrorEstimator import ErrorEstimator
import json
import sys

def main():
	with open(sys.argv[1]) as f:
		data = json.load(f)


	for p in data['data']:
		err = ErrorEstimator(p["PROCESS"], p["ENERGY"], p["SIGMA"], p["RHO"], p["DSIGMA"], p["DRHO"], 100)
		values, errors = err.main()

if __name__ == '__main__':
	main()