.PHONY: test

test:
	# $(ROOT_PYTHON) -m unittest test.test_impact_amplitude
	$(ROOT_PYTHON) -m unittest discover 