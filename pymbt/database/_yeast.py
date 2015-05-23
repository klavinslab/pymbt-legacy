from intermine.webservice import Service
import pymbt
# TODO: Use httplib instead if we only need to do one requests-style function
try:
    import requests
except ImportError:
    print 'requests module could not be imported, fetching sequences will ' + \
          'not work.'
    pass



def fetch_yeast_locus_sequence(locus_name, flanking_size=0):
    '''Acquire a sequence from SGD http://www.yeastgenome.org.

    :param locus_name: Common name or systematic name for the locus (e.g. ACT1
                       or YFL039C).
    :type locus_name: str
    :param flanking_size: The length of flanking DNA (on each side) to return
    :type flanking_size: int

    '''
    service = Service('http://yeastmine.yeastgenome.org/yeastmine/service')

    # Get a new query on the class (table) you will be querying:
    query = service.new_query('Gene')

    # The view specifies the output columns
    # secondaryIdentifier: the systematic name (e.g. YFL039C)
    # symbol: short name (e.g. ACT1)
    # length: sequence length
    # flankingRegions.direction: Upstream or downstream (or both) of locus
    # flankingRegions.sequence.length: length of the flanking regions
    # flankingRegions.sequence.residues: sequence of the flanking regions
    query.add_view('secondaryIdentifier', 'symbol', 'length',
                   'flankingRegions.direction',
                   'flankingRegions.sequence.length',
                   'flankingRegions.sequence.residues')

    # You can edit the constraint values below
    query.add_constraint('flankingRegions.direction', '=', 'both',
                         code='A')
    query.add_constraint('Gene', 'LOOKUP', locus_name, 'S. cerevisiae',
                         code='B')
    query.add_constraint('flankingRegions.distance', '=',
                         '{:.1f}kb'.format(flanking_size / 1000.),
                         code='C')
    # Uncomment and edit the code below to specify your own custom logic:
    query.set_logic('A and B and C')

    # TODO: What to do when there's more than one result?
    first_result = query.rows().next()
    # FIXME: Use logger module instead
    # print first_result['secondaryIdentifier']
    # print first_result['symbol'], row['length']
    # print first_result['flankingRegions.direction']
    # print first_result['flankingRegions.sequence.length']
    # print first_result['flankingRegions.sequence.residues']

    seq = pymbt.DNA(first_result['flankingRegions.sequence.residues'])
    # TODO: add more metadata

    return seq


def get_yeast_sequence(chromosome, start, end, reverse_complement=False):
    '''Acquire a sequence from SGD http://www.yeastgenome.org
    :param chromosome: Yeast chromosome.
    :type chromosome: int
    :param start: A biostart.
    :type start: int
    :param end: A bioend.
    :type end: int
    :param reverse_complement: Get the reverse complement.
    :type revervse_complement: bool
    :returns: A DNA sequence.
    :rtype: pymbt.DNA

    '''
    if start != end:
        if reverse_complement:
            rev_option = '-REV'
        else:
            rev_option = ''
        param_url = '&chr=' + str(chromosome) + '&beg=' + str(start) + \
                    '&end=' + str(end) + '&rev=' + rev_option
        url = 'http://www.yeastgenome.org/cgi-bin/getSeq?map=a2map' + \
            param_url
        res = requests.get(url)
        # ok... sadely, I contacted SGD and they haven;t implemented this so
        # I have to parse their yeastgenome page, but
        # it is easy between the raw sequence is between <pre> tags!

        # warning that's for the first < so we need +5!
        begin_index = res.text.index('<pre>')

        end_index = res.text.index('</pre>')
        sequence = res.text[begin_index + 5:end_index]
        sequence = sequence.replace('\n', '').replace('\r', '')
    else:
        sequence = ''

    return pymbt.DNA(sequence)

def get_yeast_gene_location(gene_name):
    """Acquire the location of a gene from SGD http://www.yeastgenome.org
    :param gene_name: Name of the gene.
    :type gene_name: string
    :returns location: [int: chromosome, int:biostart, int:bioend, int:strand]
    :rtype location: list
    
    """
    service = Service("http://yeastmine.yeastgenome.org/yeastmine/service")

    # Get a new query on the class (table) you will be querying:
    query = service.new_query("Gene")
    
    # The view specifies the output columns
    query.add_view("primaryIdentifier", "secondaryIdentifier", "symbol", 
                   "name", "organism.shortName", 
                   "chromosome.primaryIdentifier",
                   "chromosomeLocation.start", "chromosomeLocation.end",
                   "chromosomeLocation.strand")
    
    # Uncomment and edit the line below (the default) to select a custom sort order:
    # query.add_sort_order("Gene.primaryIdentifier", "ASC")
    
    # You can edit the constraint values below
    query.add_constraint("organism.shortName", "=", "S. cerevisiae", code = "B")
    query.add_constraint("Gene", "LOOKUP", gene_name, code = "A")
    
    # Uncomment and edit the code below to specify your own custom logic:
    # query.set_logic("A and B")
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
                   "chrXVI": 16
                   }
    

    return  [chromosomes[first_result["chromosome.primaryIdentifier"]], 
             first_result["chromosomeLocation.start"], 
             first_result["chromosomeLocation.end"], 
             int(first_result["chromosomeLocation.strand"])]