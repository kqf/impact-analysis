.PHONY: test all clean

all:
	python impact config/input.json

test:
	python -m unittest discover 	


clean:
	rm -f *.eps *.dat 
	rm -f output/*.eps output/*.dat	
