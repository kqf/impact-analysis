.PHONY: test all clean

main:
	python -m unittest analysis.run

lower:
	python -m unittest analysis.low

test:
	python -m unittest discover 	


clean:
	rm -f *.eps *.dat 
	rm -f output/*.eps output/*.dat	
