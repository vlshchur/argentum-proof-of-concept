# stable
1. Brief file description.
analyzeARG.py - read trees in Newick format (indeed MACS output, ignoring SITEs), computes the minimal branch size density (average or pairwise)

argentum.cpp - reads a "binary" file and builds a tcPBWT-ARG

binaryTree.py - generates random binary trees under a well-mixture assumption and counts minimal branch sizes

macs-edit.pl - mainly used to convert MACS-file to a "binary" file (sites represented as strings of 0 and 1)

vcf2bin.pl - read vcf and outputs all lines where a value of AA field is known as a string of 0 (ancestral allele) and 1 (derived allele)


2. Argentum: preparing and processing your data.

2.1. Argentum supports only a specific following file format. The first line of a data file is arbitrary and is not processed. So it can be used as a description. Next lines correspond to sites and contatain an equal number of 0 and 1, where 0 corresponds to an ancestral allele and 1 corresponds to a derived allele. The processing ends when either EOF or an empty line is reached. So empty line can be followed with any additional comments which will be ignored by the program.

2.2 The scripts to recode vcf and macs data are provided. Be aware that AA field should be filled in vcf!

Recode macs: ./macs-edit.pl [filename] -bin    Output file "out-macs.txt"

Recode vcf:  ./vcf2bin.pl [filename]           Output file "out-bin.txt"
