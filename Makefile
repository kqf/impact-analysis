.PHONY: test all clean

all:
	-mkdir -p output
	python -m unittest analysis.run
	-mv *.eps output
	-mv *.pdf output
	-mv *.csv output
	-mv *.tex output

test:
	python -m unittest discover 	


clean:
	rm -f *.eps *.dat 
	rm -f output/*.eps output/*.dat	
