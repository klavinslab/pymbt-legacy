# TODO List

## Fixes:
* Features should be updated on \_\_getattr\_\_

## OO API:
* Make sequences (DNA, RNA, Peptide) inherit from MutableSequence ABC

## Organization:

## Wishlist:
* High Priority
    * Analysis:
        * Read/write sequencing alignments (Sanger methods)
            * Need a format!
        * ViennaRNA
            * Equivalent (to Nupack) functions for:
                * complexes
                * concentrations
            * Return same format for pairs as Nupack does
            * Add tests
    * Design:
    * Reaction:
        * Gibson
    * SeqIO:
    * Sequence:
* Low Priority
    * Analysis:
        * Drop emboss dependency by using pycogent / pycogent approach for needleman-wunsch
            * pycogent method (pure python) is significantly slower than emboss
            * optimize pycogent method or try numpy or cython
        * Nupack's method outputs could be cleaner. They're usually dicts and
          some return a small fraction of the useful output of nupack
    * Data:
        * Scientific database access object(s)
            * NEB API for restriction enzymes?
            * SGD, EcoCYC
    * Design:
        * Repeat detection for OligoAssembly
        * Automatic sequencing primer design
    * Reaction:
        * More complex restriction digest (checks for all possible cut orders,
        if ambiguities are possible)
    * SeqIO:
    * Sequence:
        * special genome object
        * comparisons - current behaviors to question:
            * != when features are different
            * != when comparing to a string
    * Visualization:
        * New module? How do you visualize a sequence, reaction, design, etc?
