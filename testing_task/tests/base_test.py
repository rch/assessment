import os, sys, six, unittest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../lib'))
from collections import namedtuple
from parse.gtf import GtfParser, gtf_to_gene
from annotation import Feature
import annotate

Opts = namedtuple('Opts', ['coordinate_file','outfile','gtf_file',
                           'chromosome_indices','coordinate_indices'])
 
class BaseTest(unittest.TestCase):
    
    def __init__(self, *args, **kwargs): 
        super(BaseTest, self).__init__(*args, **kwargs)
    
    def setUp(self):
        opts = {'coordinate_file': 'master_file_unannotated.txt',
            'outfile': 'master_file_annotated.txt',
            'gtf_file': 'hg19_annotations.gtf',
            'chromosome_indices': '1, 1',
            'coordinate_indices': '5, 6'}
        self.opts = Opts(**opts)
        self.chromosomes, self.gene_spans = gtf_to_gene(open(self.opts.gtf_file, 'U')) 
        self.coordinate_indices = map(int, self.opts.coordinate_indices.strip().split(', '))
        self.chromosome_indices = map(int, self.opts.chromosome_indices.strip().split(', '))
        self.master_file_path = self.opts.coordinate_file
        self.out_path = self.opts.outfile
    
    def test_annotate_master_file(self):
        """Not a great test. QueryTest is more granular.
        """
        annotate.annotate_master_file(master_file_path=self.opts.coordinate_file,
                         out_path=self.opts.outfile,
                         chromosomes=self.chromosomes,
                         gene_spans=self.gene_spans,
                         coordinate_indices=self.coordinate_indices,
                         chromosome_indices=self.chromosome_indices)
        last_line = "M00517:73:000000000-A3H68:1:2114:15894:27217    chr21   -   1   80  9827137 9827214 Novel:Novel:n/a:n/a Novel:Novel:n/a:n/a"
        with open(self.opts.outfile) as f:
            for line in f: pass
        self.assertTrue(line, last_line)
