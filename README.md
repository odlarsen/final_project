# User Manual
This program is designed to take in 4 sequence alignments of genes within a bacteria type in fasta format and MIC antibiotic resistance values in a csv file.  It will find mutations that have the most correlation with antibiotic resistance.  Download the python code and sequence alignments for the test data.  There are two sets of test data, one for Staphylococcus aureus that begin with Sa and one for Staphylococcus haemolyticus that start with Sh.  

# Getting Started
The paths within the code indicating file location should be altered to their location on your machine.  You also must change the paths of the output files.  To run this code, open it in an IDE and run it or open a command prompt with python installed and use "python mic_alleles.py" to run the code.  Biopython, csv, itertools, pandas, and scipy must be installed on your machine.

In order to run this code for Staphylococcus aureus, switch out the names of Staphylococcus haemolyticus files.  The outfile names are also changable to what ever is desired.

# Output
Two output files will be written.  One will contain the number of unique alleles and the alleles and the other will contain a chi-squared value for which mutation(s) have an effect on MIC value, the positions of mutations relative to each other, the genes where the mutation(s) are located, and the amino acid present at each mutation for every sequence.

# Runtime Considerations
Runtime is dependent on the amount of mutations and sequences within an alignment.
