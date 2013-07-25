# TODO List

## Fixes:
## Organization:
## Wishlist:
* High Priority
    * Multi-functional class API:
        * sequence are inputted on init
        * calls to methods:
            * should return output and not write public attributes if calculations are fast
            * should also write public attribute(s) if calculations are slow
        * classes that violate these rules:
          * StructureWindows
          * OligoAssembly
          * Gibson (reaction)
    * Design:
        * Gibson
    * Reaction:
        * Ambiguous Gibson Detection
    * Sequence:
* Medium Priority
    * Analysis:
        * Read/write sequencing alignments (Sanger methods)
            * Need a format!
    * Sequence:
        * handle enzymes that cut outside of recognition site
        * molecular weight method for gels with digested DNA, gapped DNA, etc.
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
