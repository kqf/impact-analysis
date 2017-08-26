.PHONY: test all clean

all:
	$(ROOT_PYTHON) impact config/input.json

test:
	$(ROOT_PYTHON) -m unittest discover 	


clean:
	rm -f *.eps *.dat 
	rm -f output/*.eps output/*.dat	
