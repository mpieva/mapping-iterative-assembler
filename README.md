mapping-iterative-assembler
===========================

The basic idea of this program is to align DNA sequencing fragments
(shotgun or targeted resequencing) to a reference, then call a
consensus.  Then the consensus is used as new reference and the process
is repeated until convergence.  Since it was originally designed to be
used on ancient DNA, it supports a position specific substitution
matrix, which improves both alignment and consensus calling on
chemically damaged aDNA.

MIA has been used to assemble a number of Neandertal and early modern
human mitochondria.   Occasionally it has been used on smallish nuclear
regions, but it will probably not scale to a genome wide analysis.


contamination-checker
=====================

This program takes the output of MIA and tries to make sure an assembled
mitochondrion is free from contamination.  It works by looking for
positions in the called consensus where it differs from a panel of known
human mitochondria, then classifies each read as either belonging to the
sample or a putative contaminant.

