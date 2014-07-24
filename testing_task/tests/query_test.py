import os, sys, six, unittest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../lib'))
from annotation import Feature
import annotate

 
class QueryTest(unittest.TestCase):
    
    def __init__(self, *args, **kwargs): 
        super(QueryTest, self).__init__(*args, **kwargs)
    
    def setUp(self):
        self.query_feature = Feature(start=20704380, end=20704380, chromosome='chr12') 
    
    def test_find_gene_overlaps(self):
        """Find gene overlaps for a known query.
        """
        chrom = self.query_feature.Chromosome
        test_feature = Feature(start=20522179, end=20837042) # should overlap
        query_feature, gene_features = annotate.find_gene_overlaps(chrom=chrom,
            start=self.query_feature.Start,
            end=self.query_feature.End,
            gene_spans={chrom: [test_feature]}
        )
        self.assertTrue(gene_features)
    
    def test_annotation_from_gene_ids(self):
        """Check annotation from gene IDs with a known query and gene ID.
        """
        chrom = self.query_feature.Chromosome
        gene_id = 'PDE3A'
        test_feature = Feature(gene_id=gene_id, start=20523179, end=20709594)
        annotation = annotate.annotation_from_gene_ids(self.query_feature,
                                              chromosomes={chrom:{gene_id:[test_feature]}},
                                              chrom=chrom,
                                              gene_features=
                                              [Feature(gene_id=gene_id, start=20522179, end=20837042)])  
        self.assertTrue(isinstance(annotation, six.string_types))
        self.assertTrue(annotation.find(test_feature.GeneId))
        