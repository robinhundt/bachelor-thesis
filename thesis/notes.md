# Notes on thesis
Containing notes and ramblings on this thesis


## Abstract

## Introduction

- maybe this could include reasoning that majority of time in Dialign is spend constructing pairwise alignments and only about 1-2% is spent adding diagonals to closure
    - finding enough suitable diagonals **faster** and greedily aligning them could drastically improve runtime

## Basics

### MSA
- basis for numerous further analysis such as "inferring phylogenetic relationships, homology search of functional elements, classification of proteins, designing detection markers" \cite{Russell2016}

- describe multiple sequence alignments and their importance


## Prior Work
- Spaced Word Matches
- Dialign
    - Gabios-Lib (as part of Dialign 2.2)

## Algorithm

- more detailed description of gabios lib 

## Implementation
- reimplementation of finding spaced word matches and GABIOS LIB
- implementation in Rust

## Evaluation

### Dataset(s)
#### Balibase
- use version 3 -> why?
    - comparability
- describe balibase and it's characteristics -> flaws in core block definition?


### MAFFT

- Describe process of evaluation
    - which dataset(s) have been used a

## Conclusion

### Further Work




# TODOs

What is the objective/scoring function????
Change output so that unaligned pos are lowercase

## IMPORTANT pdf/a validation
The thesis has to be in PDF/A format for archiving. This template automatically generates your file according to the PDF/A standard. It can be that this is violated by changing the document, for example by inserting non-compliant graphics. Please check, whether your final version of the document is PDF/A compliant, e.g., by using a freely available tool like \emph{veraPDF}\footnote{https://verapdf.org/home/}.


# Reading
- use maximal exact matches (MEMs) that are locally longest common subse- quences without the uniqueness condition of MUM. MUMs and MEMs can be efficiently detected by a suffix tree, suffix array, or simple k-mer table [82] Crochemore M, Hancart C, Lecroq T (2007) Algorithms on strings. Cambridge University Press, Cambridge


# Ideas
- Extending found spaced word matches
- provide multiple pattern sets with different weight ratios, first construct diagonals for pattern sets with high weights -> add diagonals above $THRESHOLD but retain those below -> construct diagonals with lower weight ratio (this should be faster since earlier diagonals provide anchors and search space is smaller) -> merge new diagonals with earlier ones -> add to alignment
    - is threshold needed? Pattern with high enough weight provide indirect threshold anyway