# TODO List

## Fixes:
* Features should be updated on \_\_getitem\_\_
* Features should be renamed when truncated

## OO API:
* Make sequences (DNA, RNA, Peptide) inherit from MutableSequence ABC

## Organization:

## Wishlist:
* High Priority
    * In general, replace .run stuff with methods
        * classes to replace with functions:
            *Gibson? 
    * Design:
        * Gibson
    * Reaction:
        * Remove .run methods. New idea: split into reaction objects/types and
          reaction functions. Reaction functions do reaction on 'init' and
          return reaction objects/types.
    * Sequence:
        * Features should update whenever sequence changes
* Medium Priority
    * Analysis:
        * Read/write sequencing alignments (Sanger methods)
            * Need a format!
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
    * Sequence:
        * special genome object
        * comparisons - current behaviors to question:
            * != when features are different
            * != when comparing to a string
    * Visualization:
        * New module? How do you visualize a sequence, reaction, design, etc?
