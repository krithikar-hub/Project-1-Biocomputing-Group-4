Pipeline:-
•	Obtained .fasta file using API
The nucleotide sequences were downloaded from NCBI using Python and the `requests` library.
•	Obtained genebank file using API and reported the accession id
The genebank files were downloaded from NCBI using Python and the `requests` library.
•	Obtained protein sequences
The protein sequences were downloaded from NCBI using Python and the `requests` library.
•	Extracted fasta headers
grep "^>" DIO2.fasta DIO2_protein.fasta > headers_DIO2.txt
head -n 1 MT_ATP6_CDS.fasta > headers.txt
•	Extracted CDS
CDS was extracted using the position
•	Translated the CDS
The CDS was translated using the standard genetic code in Python 
•	Cross-checked my translated protein with the deposited protein
The CDS was cross-checked with deposited protein using for loop and if-else conditioning
•	Use the translate function in Biopython and cross-check with my translated protein
The CDS was translated and cross-checked using Biopython.
•	Saved the results in .txt file
•	Git clone and git pull the files
git clone <repository-url> 
git add . 
git commit -m "Add your_gene results" 
git push
•	Created a summary table for all 3 genes
Saved the relevant informations in .txt files and defined a function to call all the genes to create the summary table

Members	Category	Genes
Rohini Chakraborty	A (Selenoprotein)	DIO2
Krithika Rajagopal	B (Mitochondrial)	MT-ATP6
Naresh Abhange  	C (Nuclear control)	PGK1
