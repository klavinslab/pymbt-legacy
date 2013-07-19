# TODO List

## Fixes:
* Features should be updated whenever sequence changes:
    * on \_\_getitem\_\_
* Features should be renamed when truncated

## OO API:
* Make sequences (DNA, RNA, Peptide) inherit from MutableSequence ABC

## Organization:

## Wishlist:
* High Priority
    * In general, replace .run classes with functions 
        * classes to replace with functions:
            *Gibson? 
    * Classes / etc that have inconsistent APIs:
        * Nupack
        * Vienna
        * StructureWindows
        * OligoAssembly
        * Gibson (reaction)
    * Design:
        * Rename 'design_primer' to 'primer' and 'design_primer_pcr' to 'primers'
        * Gibson
    * Reaction:
        * Ambiguous Gibson Detection
    * Sequence:
        * May want a base class for DNA/RNA. Things common to both:
            * Both have init sequence that gets checked for alphabet
            * Both have option to disable checks
            * reverse complement method is identical but for ss/ds branch
            * both have locate method which could share some but not all code
            * __getitem__ identical but for topology change
            * __delitem__ should be identical but isn't
            * __setitem__ should be identical exept for rev comp of bottom strand for DNA
            * __repr__ is identical except DNA also reports topology
            * __str__ is identical
            * __len__ is identical
            * __add__ is much more complex for DNA than RNA due to possibility of ds gaps
            * __radd__ should be identical
            * __mul__ should be identical except for topology branch
            * __eq__ is identical
            * __ne__ is identical
            * __contains__ is identical
        * Peptide shares these with RNA (basically everything!): 
            * locate method is identical to RNA's
            * copy method is identical to RNA's
            * __getitem__ is identical except Peptide ignores bottom strand
              (RNA currently has no purpose for a bottom strand either!!!!!!)
            * __delitem__ is identical except for bottom strand
            * __setitem__ is identical
            * __repr__ is identical except for bottom strand / words
            * __str__ is identical
            * __len__ is identical
            * __add__ is identical
            * __radd__ should be identical
            * __mul__ is identical
            * __eq__ is identical
            * __ne__ is identical
            * __str__ is identical

* Medium Priority
    * Tests!
        * design primer with overhangs
        * gibsons
        * test for feature extraction working
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
