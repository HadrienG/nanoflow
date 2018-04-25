clean:
	rm -rf work/
	rm -rf .nextflow.log*
	cd test && rm -rf work && rm -rf .nextflow.log* && rm -rf results
