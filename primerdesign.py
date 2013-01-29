from tm_calc import calc_tm
from numpy import array as narray

# Should just return a TmRecord, doesn't need a new class.

def designPrimer(s,tm=72,minstart=10,tm_errorminus=1,tm_errorplus=3,endGC=True):
    # Check Tm to see if it's already too low
    source_tm = calc_tm(s)

    if source_tm < tm-tm_errorminus:
        return("Tm of full-length sequence is already lower than desired Tm and error parameters allow.")

    primer_tm = 0
    bases = minstart
    s = s[0:70] # Trim down max length to increase efficiency
    primers = []
    tms = []

    # Storing primer sequence separate from Tm may increase efficiency
    # First, generate all primers up to Set Tm + Allowed Error. Even though it's iterating, it's faster than generating
    # all reasonably-sized oligos and finding their tms
    while primer_tm <= (tm+tm_errorplus) and bases <= len(s):
        # Increase oligo length
        primers.append(s[0:bases])
        tms.append(calc_tm(primers[-1]))
        primer_tm = tms[-1]

        bases += 1
    if primer_tm <= tm+tm_errorplus:
        return("Could not generate an oligo with high enough Tm")
    # Trim to primers matching minimum tm
    tm_ind = narray(tms)>=tm-tm_errorminus
    primers = narray(primers)[tm_ind].tolist()
    tms = narray(tms)[tm_ind].tolist()

    primerlist = zip(primers,tms)

    gc_list = []
    if endGC:
        for i,x in enumerate(primerlist):
            primer_current = x[0]
            # TODO:
            # Fix this '.upper()' hack - use REs to check end of primer
            if primer_current.upper().endswith('G') or primer_current.upper().endswith('C'):
                gc_list.append(x)
        if len(gc_list) == 0:
            return("No primers in this range ending in G or C could be generated")
        primerlist = gc_list


    # which primer(s) are closest to the desired tm?
    diffs = []
    for i in range(len(primerlist)):
        diffs.append(abs(72-primerlist[i][1]))

    min_ind = [i for i, x in enumerate(diffs) if x==min(diffs)]
    best_primer = primerlist[min_ind[0]]
    best_primer = calc_tm(best_primer[0])

    return(best_primer)

if __name__ == "__main__":
    s = raw_input('Input Sequence:')
    print "Your primer:"
    print ""
    out = designPrimer(s)
#    assert Tm_staluc('CAGTCAGTACGTACGTGTACTGCCGTA') == 59.865612727457972

