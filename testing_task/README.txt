Requirements:

	python2.7 (http://www.python.org/download/releases/2.7/)
	PyCogent (http://pycogent.org)


Please ensure the quality of the following code:

	lib/annotation.py
	lib/annotate.py


Code description:

Given coordinates and a chromosome, this program for looks up it's annotation and annotates an input file.
	- Input:
		- Tab-delimited file: master_file_unannotated.txt
		- GTF formatted file with genome annotations: **NOTE** unpack hg19_annotations.tgz
	- Output: 
		- Annotated file of gene name that input position overlaps: master_file_annotated.txt


The following command generates annotated output.

$ python lib/annotate.py --gtf_file hg19_annotations.gtf --coordinate_file master_file_unannotated.txt -o master_file_annotated.txt