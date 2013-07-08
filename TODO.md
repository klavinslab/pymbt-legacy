#TODO List

## Fixes:
* Things that are not working (and should have tests written for them):
  * Plotting sequencing analysis
    * DNA objects need sequence feature model first

## OO API:
* Make sequences (DNA, RNA, Peptide) inherit from MutableSequence ABC

## Organization:
* Remove features that don't contribute to core workflow?:
  * Gene splitting

## Wishlist:
* Add Breslauer, SantaLucia98, and Sugimoto methods
 See doi: 10.1093/bioinformatics/bti066 for good comparison
   "Comparison of different melting temperature calculation
    methods for short DNA sequences"


* Choose best primer from database for sequencing
* Repeat detection for OligoAssembly
* Check cloning primers for self-self, pair structure using Nupack
* Check Gibson object for whether the desired overlap is a high-probability structure
* Add RNAfold / RNAStructure / Vienna as drop-in NUPACK replacements
* Drop emboss dependency by using pycogent / pycogent approach for needleman-wunsch
  * pycogent method (pure python) is significantly slower than emboss
  * try numpy or cython
