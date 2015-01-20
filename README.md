# stable
aa-repair.pl - read vcf and outputs all lines where a value of AA field is known as a string of 0 (ancestral allele) and 1 (derived allele)
analyzeARG.py - read trees in Newick format (indeed MACS output, ignoring SITEs), computes the minimal branch size density (average or pairwise)
binaryTree.py - generates random binary trees under a well-mixture assumption and counts minimal branch sizes
macs-edit.pl - mainly used to convert MACS-file to a "binary" file (sites represented as strings of 0 and 1)
pbwt-tree-reconstructor.cpp - reads a "binary" file and builds a tcPBWT-ARG
vcf2bin.pl - read vcf and outputs all lines where a value of AA field is known as a string of 0 (ancestral allele) and 1 (derived allele)
