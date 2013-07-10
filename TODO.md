#TODO List

## Fixes:
* Things that are not working (and should have tests written for them):
    * Plotting sequencing analysis
        * DNA objects need sequence feature model first

## OO API:
* Make sequences (DNA, RNA, Peptide) inherit from MutableSequence ABC

## Organization:

## Wishlist:
* High Priority
    * Drop emboss dependency by using pycogent / pycogent approach for needleman-wunsch
        * pycogent method (pure python) is significantly slower than emboss
        * try numpy or cython
    * Add ViennaRNA as partial NUPACK replacement
    * Check cloning primers for self-self, pair structure using Nupack
    * Make Gibson reaction object
        * Check Gibson object for whether the desired overlap is a high-probability structure
* Low Priority
    * Repeat detection for OligoAssembly
