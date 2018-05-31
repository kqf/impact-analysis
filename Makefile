.PHONY: test all clean

all:
	python -m unittest analysis.run

test:
	python -m unittest discover 	


clean:
	rm -f *.eps *.dat 
	rm -f output/*.eps output/*.dat	
