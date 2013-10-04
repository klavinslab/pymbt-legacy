
Primer Design
=============

One of the first things anyone learns in a molecular biology lab is how
to design primers. The exact strategies vary a lot and are sometimes
polymerase-specific. ``pymbt`` uses the Klavins' lab approach of
targeting a specific melting temperature (Tm) and nothing else, with the
exact Tm targeted being between 65°C and 72°C, the choice being personal
preference. ``pymbt`` currently defaults to 72°C on the Phusion
(modified Breslauer 1986) Tm calculator.

``pymbt.design_primer`` is a function that takes in a ``sequence.DNA``
object and rapidly finds the 5' subsequence that is closest to the
desired Tm (within a user-definable error range). If the entire sequence
would make a primer with too low of a Tm, a descriptive error is
produced.

For this tutorial, let's design primers that will amplify the gene EYFP.

.. code:: python

    from pymbt import design, seqio, sequence
First we read in a plasmid from Havens et al. 2012 and isolate the EYFP
sequence.

.. code:: python

    plasmid = seqio.read_dna("./files_for_docs/maps/pGP4G-EYFP.ape")
    eyfp = plasmid.extract("EYFP")
    print len(eyfp)
    eyfp

.. parsed-literal::

    717




.. parsed-literal::

    linear dsDNA:
    atggtgagcaagggcgaggagctgttcaccggggtggtgc ... cgccgccgggatcactctcggcatggacgagctgtacaag
    taccactcgttcccgctcctcgacaagtggccccaccacg ... gcggcggccctagtgagagccgtacctgctcgacatgttc



Designing primers is straightforward - you just call
``design.design_primer`` with a ``sequence.DNA`` object as the input.

.. code:: python

    # Forward and reverse, one at a time using design_primer()
    forward = design.design_primer(eyfp)
    reverse = design.design_primer(eyfp.reverse_complement())
    # Both at once using design_primers()
    forward, reverse = design.design_primers(eyfp)
    # design_primer has many options, including adding overhangs
    custom_forward = design.design_primer(eyfp, tm=65, min_len=12, 
                                          tm_undershoot=1, tm_overshoot=1, 
                                          end_gc=True, tm_parameters="santalucia98", 
                                          overhang=sequence.DNA("GGGGGATCGAT"))
    print forward
    print
    print custom_forward

.. parsed-literal::

    atggtgagcaagggcgaggag
    
    gggggatcgatatggtgagcaagggcgaggagctgttcac

