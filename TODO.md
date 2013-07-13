# TODO List

## Fixes:

## OO API:
* Make sequences (DNA, RNA, Peptide) inherit from MutableSequence ABC

## Organization:

## Wishlist:
* High Priority
    * Analysis:
        * Read/write sequencing alignments (Sanger methods)
        * Add ViennaRNA as partial NUPACK replacement
        * Check cloning primers for self-self, pair structure using Nupack
    * Data:
        * Scientific database access object(s)
    * Design:
        * Given a sequence with a subsequence to a (Sanger) sequence,
          design primers to give full coverage of the subsequence.
    * Reaction:
        * PCR
        * Gibson
            * Check for structures / alternative binding
        * BASIC Restriction digest
    * SeqIO:
        * write function
    * Sequence:
        * write methods per-object, rather than in seqio?
        * 'read' function defined for sequence module rather than seqio (to
          avoid circular deps)?
    * Visualization:
        * New module? How do you visualize a sequence, reaction, design, etc?
* Low Priority
    * Analysis:
        * Drop emboss dependency by using pycogent / pycogent approach for needleman-wunsch
            * pycogent method (pure python) is significantly slower than emboss
            * optimize pycogent method or try numpy or cython
        * Nupack's method outputs could be cleaner. They're usually dicts and
          some return a small fraction of the useful output of nupack
    * Data:
    * Design:
        * Repeat detection for OligoAssembly
    * Reaction:
        * More complex restriction digest (checks for all possible cut orders,
        if ambiguities are possible)
    * SeqIO:
    * Sequence:
        * special genome object
        * comparisons - current behaviors to question:
            * != when features are different
            * != when comparing to a string
