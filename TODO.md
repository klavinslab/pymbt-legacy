# TODO List

## Fixes:
## Organization:
## Wishlist:
* High Priority
    * Sanger report() on Yaoyu's pGP4-GalZif-ZifAVNY doesn't include a T in the sequence near the insertion
    * Oligo assembly reaction
    * Recognize all of genbank feature types, read/write them
* Medium Priority
    * Read/write sequencing alignments
        * Need a format! There are a lot. How to choose? Make my own? Just pick one that does the job for now.
    * Sequence:
        * handle enzymes that cut outside of recognition site
        * molecular weight method for gels with digested DNA, gapped DNA, etc.
* Low Priority
    * Drop emboss dependency by using pycogent / pycogent approach for needleman-wunsch
        * pycogent method (pure python) is significantly slower than emboss
        * optimize pycogent method or try numpy or cython
    * Nupack's method outputs could be cleaner. They're usually dicts and some return a small fraction of the useful output of nupack
    * Scientific database access object(s)
        * SGD, EcoCYC
    * Repeat detection for OligoAssembly
    * Automatic sequencing primer design
    * More complex restriction digest (checks for all possible cut orders, if ambiguities are possible)
    * Visualizations for a sequence, reaction, design, etc?
