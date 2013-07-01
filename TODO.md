#TODO List

## Fixes:
* Things that are not working (and should have tests written for them):
  * Plotting sequencing analysis
    * DNA objects need sequence feature model first

## Convert to oo API:
* Make sequences (DNA, RNA, Peptide) inherit from MutableSequence ABC
* Remove use of sequence.utils outside of sequence module, if possible
  * Functions that currently require use of strings for DNA should be addressed as follows:
    1. Replace with DNA objects - if they're too slow, optimize DNA object methods.
    2. Remove from module (put in an extension)

## Organization:
* Remove features that don't contribute to core workflow?:
  * OrthoSeq
  * NUPACK
  * Gene splitting

## Wishlist:
* Add Breslauer, SantaLucia98, and Sugimoto methods
 See doi: 10.1093/bioinformatics/bti066 for good comparison
   "Comparison of different melting temperature calculation
    methods for short DNA sequences"


* Choose best primer from database for sequencing
* Repeat detection for OrthoSeq, OligoAssembly
* Check cloning primers for self-self, pair structure using Nupack
* Check Gibson object for whether the desired overlap is a high-probability structure
* Add RNAfold / RNAStructure / Vienna as drop-in NUPACK replacements
* Drop emboss dependency by using pycogent / pycogent approach for needleman-wunsch
  * pycogent method (pure python) is significantly slower than emboss
  * try numpy and / or c
* Implement santa lucia and sugimoto params in Tm
