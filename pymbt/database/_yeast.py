"""Yeast database query functions."""
from intermine.webservice import Service
import pymbt
try:
    import requests

except ImportError:
    print "requests module could not be imported, fetching sequences will " + \
          "not work."


class SGD(object):
    def __init__(self):
        baseurl = "http://yeastmine.yeastgenome.org/yeastmine"
        self.service = Service(baseurl + "/service")
        self.getsequrl = "http://www.yeastgenome.org/cgi-bin/getSeq"

    def sequence_from_name(self, locus_name, upstream=0, downstream=0,
                           reverse_complement=False):
        """Acquire a sequence from SGD http://www.yeastgenome.org.

        :param locus_name: Common name or systematic name for the locus (e.g.
                           ACT1 or YFL039C).
        :type locus_name: str
        :param flanking_size: The length of flanking DNA (on each side).
        :type flanking_size: int

        """
        params = {"map": "a2map",
                  "seqname": locus_name,
                  "flankl": upstream,
                  "flankr": downstream}
        if reverse_complement:
            params["rev"] = "-REV"

        res = requests.get(self.getsequrl, params=params)
        # OK... sadly, I contacted SGD and they haven;t implemented this so
        # I have to parse their yeastgenome page, but
        # it is easy between the raw sequence is between <pre> tags!

        seq_text = res.text.split("<pre>")[1].split("</pre>")[0]
        sequence = seq_text.replace("\n", "").replace("\r", "")
        return pymbt.DNA(sequence)

    def sequence_from_coordinates(self, chromosome, start, end,
                                  reverse_complement=False):
        """Acquire a sequence from SGD http://www.yeastgenome.org.

        :param chromosome: Yeast chromosome.
        :type chromosome: int
        :param start: Start location on the chromosome.
        :type start: int
        :param end: End location on the chromosome.
        :type end: int
        :param reverse_complement: Get the reverse complement.
        :type revervse_complement: bool
        :returns: A DNA sequence.
        :rtype: pymbt.DNA

        """
        if start == end:
            raise ValueError("start and end must be different")

        params = {"map": "a2map",
                  "chr": chromosome,
                  "beg": start,
                  "end": end}
        if reverse_complement:
            params["rev"] = "-REV"

        res = requests.get(self.getsequrl, params=params)

        seq_text = res.text.split("<pre>")[1].split("</pre>")[0]
        sequence = seq_text.replace("\n", "").replace("\r", "")
        return pymbt.DNA(sequence)

    def location_from_name(self, locus_name):
        """Acquire the location of a gene from SGD http://www.yeastgenome.org
        :param locus_name: Common name or systematic name for the locus (e.g.
                           ACT1 or YFL039C).
        :type locus_name: str
        :returns location: [int: chromosome, int:biostart, int:bioend,
                            int:strand]
        :rtype location: list

        """
        query = self.service.new_query("Gene")

        # Add views
        query.add_view("chromosome.primaryIdentifier",
                       "chromosomeLocation.start", "chromosomeLocation.end",
                       "chromosomeLocation.strand")

        # Add constraints
        query.add_constraint("Gene", "LOOKUP", locus_name, code="A")
        query.add_constraint("organism.shortName", "=", "S. cerevisiae",
                             code="B")

        # Set logic
        query.set_logic("A and B")

        # Process results
        chromosomes = {"chrI": 1,
                       "chrII": 2,
                       "chrIII": 3,
                       "chrIV": 4,
                       "chrV": 5,
                       "chrVI": 6,
                       "chrVII": 7,
                       "chrVIII": 8,
                       "chrIX": 9,
                       "chrX": 10,
                       "chrXI": 11,
                       "chrXII": 12,
                       "chrXIII": 13,
                       "chrXIV": 14,
                       "chrXV": 15,
                       "chrXVI": 16}
        first_result = query.rows().next()
        chr_number = chromosomes[first_result["chromosome.primaryIdentifier"]]

        return [chr_number,
                int(first_result["chromosomeLocation.start"]),
                int(first_result["chromosomeLocation.end"]),
                int(first_result["chromosomeLocation.strand"])]

    def sysname_from_name(self, locus_name):
        """Retrieve systematic yeast gene name from the common name.

        :param locus_name: Common name or systematic name for the locus (e.g.
                           ACT1 or YFL039C).
        :type locus_name: str
        :returns: Systematic name for yeast gene (e.g. YOR128C).
        :rtype: str

        """
        query = self.service.new_query("Gene")

        # Add views
        query.add_view("primaryIdentifier", "secondaryIdentifier", "symbol",
                       "name", "sgdAlias", "crossReferences.identifier",
                       "crossReferences.source.name")

        # Add constraints
        query.add_constraint("organism.shortName", "=", "S. cerevisiae",
                             code="B")
        query.add_constraint("Gene", "LOOKUP", locus_name, code="A")

        # Set logic
        query.set_logic("A and B")

        first_result = query.rows().next()
        sysname = first_result["secondaryIdentifier"]
        return sysname


def get_yeast_promoter_ypa(locus_name):
    """Retrieve promoter from Yeast Promoter Atlas
    (http://ypa.csbb.ntu.edu.tw).

    :param locus_name: Common name or systematic name for the locus (e.g.
                       ACT1 or YFL039C).
    :type locus_name: str
    :returns: Double-stranded DNA sequence of the promoter.
    :rtype: pymbt.DNA

    """
    sgd = SGD()
    loc = sgd.location_from_name(locus_name)
    gid = sgd.sysname_from_name(locus_name)

    ypa_baseurl = "http://ypa.csbb.ntu.edu.tw/do"
    params = {"act": "download",
              "nucle": "InVitro",
              "right": str(loc[2]),
              "left": str(loc[1]),
              "gene": str(gid),
              "chr": str(loc[0])}

    response = requests.get(ypa_baseurl, params=params)
    text = response.text
    # FASTA records are just name-sequence pairs split up by > e.g.
    # >my_dna_name
    # GACGATA
    # TODO: most of this is redundant, as we just want the 2nd record
    record_split = text.split(">")
    record_split.pop(0)
    parsed = []
    for record in record_split:
        parts = record.split("\n")
        sequence = pymbt.DNA("".join(parts[1:]))
        sequence.name = parts[0]
        parsed.append(sequence)

    return parsed[1]
