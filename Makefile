.PHONY: test all

all:
	$(ROOT_PYTHON) run.py config/input.json

test:
	$(ROOT_PYTHON) -m unittest discover 	
