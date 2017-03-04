.PHONY: test all clean

all:
	$(ROOT_PYTHON) run.py config/input.json

test:
	$(ROOT_PYTHON) -m unittest discover 	


clean:
	rm -f *.eps *.dat results/*.eps results/*.dat	
