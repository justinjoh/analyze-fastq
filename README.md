### analyze-fastq
Made for the "Cpf1-CRISPR-gate-sequence" -- research work in Lambert Lab at Cornell University
##### Python 2 or 3 works
Requires Python modules:
* numpy
* matplotlib

Usage:  
Put output_SW.py in directory containing fastq files

Run from command line "python output_SW.py" to show list of fastq files to select from  
Or enter "python output_SW.py somefilenameR1.fastq somefilenameR2.fastq" to specify the files

Will display a "heatmap" displaying agreement of the SWSWSWS regions
