from intermine.webservice import Service
import pymbt
try:
    import requests
except ImportError:
    print "requests module could not be imported, fetching sequences will " + \
          "not work."
    pass


def fetch_yeast_locus_sequence(locus_name):
    """Acquire a sequence from SGD http://www.yeastgenome.org
        given a locus systematic name, get +/-1kb
    """
    service = Service("http://yeastmine.yeastgenome.org/yeastmine/service")

    # Get a new query on the class (table) you will be querying:
    query = service.new_query("Gene")

    # The view specifies the output columns
    query.add_view(
        "secondaryIdentifier", "symbol", "length", "flankingRegions.direction",
        "flankingRegions.sequence.length", "flankingRegions.sequence.residues"
    )

    # Uncomment and edit the line below (the default) to select a custom sort
    # order:

    # You can edit the constraint values below
    query.add_constraint("flankingRegions.direction", "=", "both", code="C")
    query.add_constraint("Gene", "LOOKUP", locus_name, "S. cerevisiae",
                         code="B")
    query.add_constraint("flankingRegions.distance", "=", "1.0kb", code="A")

    # Uncomment and edit the code below to specify your own custom logic:
    query.set_logic("A and B and C")

    for row in query.rows():
        # FIXME: Use logger module instead
        # print row["secondaryIdentifier"]
        # print row["symbol"], row["length"]
        # print row["flankingRegions.direction"]
        # print row["flankingRegions.sequence.length"]
        # print row["flankingRegions.sequence.residues"]

        # this is stupid but with my example I get duplicated records, I am
        # just taking the last one now....
        # https://pythonhosted.org/intermine/intermine.query.Query-class.html
        seq = pymbt.sequence.DNA(row["flankingRegions.sequence.residues"])

    return seq


def get_yeast_sequence(chromosome, start, end, rev=False):
    """Acquire a sequence from SGD http://www.yeastgenome.org
    :param chromosome: Yeast chromosome.
    :type chromosome: int
    :param start: A biostart.
    :type start: int
    :param end: A bioend.
    :type end: int
    :param rev: Get the reverse complement.
    :type rev: bool
    :returns: A DNA sequence.
    :rtype: pymbt.DNA

    """
    if start != end:
        if rev:
            reversed = "-REV"
        param_url = str(chromosome) + "&beg=" + str(start) + "&end=" + \
            str(end) + "&rev=" + reversed
        url = "http://www.yeastgenome.org/cgi-bin/getSeq?map=a2map&chr=" + \
            param_url
        res = requests.get(url)
        # ok... sadely, I contacted SGD and they haven;t implemented this so
        # I have to parse their yeastgenome page, but
        # it is easy between the raw sequence is between <pre> tags!

        # warning that's for the first < so we need +5!
        begin_index = res.text.index("<pre>")

        end_index = res.text.index("</pre>")
        sequence = res.text[begin_index + 5:end_index]
        sequence = sequence.replace('\n', '').replace('\r', '')

    else:
        sequence = ""
    return pymbt.sequence.DNA(sequence)
