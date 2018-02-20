# ARCS (Analysis of Relationship based on Chromosomal Segments)

ARCS is a nonparametric algorithm that reconstructs chromosomal segments that have 0, 1, or 2 alleles identical-by-descent.
ARCS was implemented in C++11 with only standard libraries. To compile:
```
g++ -std=c++11 ARCS.cpp -o ARCS
```

ARCS takes a number of input files:
* A comma separated file with SNP genotyping results (example_snp_genotypes)
* A list of highly polymorphic (ie. large minor allele frequency) SNP markers to examine (example_variable_sites)
* A tab separated file with features of SNP markers (example_snp_features) with the following fields:
[chromosome]  [name of the marker]  [position in cm]  [position in bp]
* A list of all samples in the dataset (example_sample_list)
* A tab separated file specifying the subject pairs to analyze (example_subject_pairs) with the following fields:
[subject 1]  [subject 2]  [a double value to indicate relationship. Note: only used for benchmarking; any double values can be used]

To run ARCS:
```
./ARCS example_snp_genotypes example_variable_sites example_snp_features example_sample_list example_subject_pairs [output file name]
```

Output:
A tab limited file with each line being a reconstructed segment. Each line has the following fields:
[subject1 subject2]  [Length of the segment]  [chromosome]  [Mean number of IBS alleles]  [relationship (or any double value specified in the last field of the subject pairs file)]  [start bp of the segment]
