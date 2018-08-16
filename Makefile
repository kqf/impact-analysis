.PHONY: test all clean

lower:
	python -m unittest analysis.low

main:
	python -m unittest analysis.run

test:
	python -m unittest discover 	


clean:
	rm -f *.eps *.dat 
	rm -f output/*.eps output/*.dat	
