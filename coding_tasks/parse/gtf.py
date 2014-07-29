#!/usr/bin/env python
# file archer.primers.parse.gtf.py

"""Parses GTF file format.
"""
from collections import defaultdict

from annotation import Chromosome, Feature


__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Development"


def split_attributes(line):
    """Returns dict of attribute label to value.
    """
    res = {}
    for attribute in line.strip().split(';'):
        if attribute:
            k, v = attribute.strip().split()
            res[k] = v.strip('"').strip("'")
    return res


def gtf_parser(lines, add_gtf_line=False):
    """Returns a dict of Chromosome objects from parsing a GTF file.
    """
    res = defaultdict(list)
    for line in MinimalGtfParser(lines):
        attrs = split_attributes(line[8])
        if add_gtf_line:
            gtf_line = '\t'.join(line).strip()
        else:
            gtf_line = None
        feature = Feature(chromosome=line[0],
                          gene_id=attrs['gene_id'],
                          feature_type=line[2],
                          start=int(line[3]),
                          end=int(line[4]),
                          strand=line[6],
                          exon_number=attrs.get('exon_number', None),
                          transcript_id=attrs.get('transcript_id', None),
                          gtf_line=gtf_line)
                
        res[line[0]].append(feature)
        
    chromosomes = {}
    for chr_name, features in res.items():
        chromosomes[chr_name] = Chromosome(name=chr_name,
                                           features=sorted(features))
    return chromosomes


def GtfParser(lines, add_gtf_line=False):
    """
    DEPRECATED:: Please use gtf_parser
    PEP8 Complinance
    """
    return gtf_parser(lines, add_gtf_line)


def minimal_gtf_parser(lines):
    """Returns a populated Genome object from parsing a GTF file.
    """
    for line in lines:
        # comments and blank lines
        if "#" in line:
            (line, comments) = line.split("#", 1)
        else:
            comments = ''
        line = line.strip()
        if not line:
            continue
        
        # parse columns
        cols = line.split('\t')
        if len(cols) == 8:
            cols.append('')
        assert len(cols) == 9, line
        (seqname, source, feature, start, end, score,
         strand, frame, attributes) = cols
        
        yield [seqname, source, feature, start, end, score,
               strand, frame, attributes, comments]


def MinimalGtfParser(lines):
    """
    DEPRECATED:: Please use minimal_gtf_parser
    PEP8 Compliance
    """
    return minimal_gtf_parser(lines)


def gtf_to_gene(lines, add_gtf_line=False, pass_back_genes_by_names=False):
    """Returns a populated Genome object from parsing a GTF file.
    """
    # chromosomes is a dict of dict of lists storing features by chromosome and
    # gene_id.
    chromosomes = defaultdict(dict)
    #gene_spans is a dict of lists storing gene spans keyed by chromosome
    gene_spans = defaultdict(list)
    
    #gene_indices is a dict of dict storing start and end indices for each gene.
    gene_indices = defaultdict(dict)
    
    # This is a dict of lists storing features keyed by genes
    gene_by_names = defaultdict(list)

    #Iterate through lines in GTF file.
    for line in MinimalGtfParser(lines):
        attrs = split_attributes(line[8])
        #Add GTF line from file if requested
        if add_gtf_line:
            gtf_line = '\t'.join(line)
        else:
            gtf_line = None
        #Cache gene_id and chromosome
        gene_id = attrs['gene_id']
        chrom = line[0]
        
        #Build feature object
        feature = Feature(chromosome=chrom,
                          gene_id=gene_id,
                          feature_type=line[2],
                          start=int(line[3]),
                          end=int(line[4]),
                          strand=line[6],
                          exon_number=attrs.get('exon_number'),
                          transcript_id=attrs.get('transcript_id'),
                          gtf_line=gtf_line)
        
        #Store gene indices
        min_idx = min(feature.Start, feature.End)
        max_idx = max(feature.Start, feature.End)
        
        #If gene hasn't been encountered, add it.
        gene_chrom = gene_id+'|'+chrom
        if gene_chrom not in gene_indices:
            gene_indices[gene_chrom]['MinIdx'] = min_idx
            gene_indices[gene_chrom]['MaxIdx'] = max_idx
            #gene_indices[gene_chrom]['Chromosome'] = chrom
        #Else accumulate current min and max indices.
        else:
            if min_idx < gene_indices[gene_chrom]['MinIdx']:
                gene_indices[gene_chrom]['MinIdx'] = min_idx
            if max_idx > gene_indices[gene_chrom]['MaxIdx']:
                gene_indices[gene_chrom]['MaxIdx'] = max_idx
            #gene_indices[gene_id]['Chromosome'] = chrom
        
        #Add current feature to chromosomes keyed by chromosome and gene_id
        if gene_id not in chromosomes[chrom]:
            chromosomes[chrom][gene_id] = [feature]
        else:
            chromosomes[chrom][gene_id].append(feature)
    #Walk through chromosomes and sort features
    for chrom, genes in chromosomes.items():
        for gene, features in genes.items():
            chromosomes[chrom][gene] = sorted(features)
            if pass_back_genes_by_names:
                gene_by_names[gene] = features

    #Walk through gene indices to build gene_spans for faster indexing
    for gene_chrom, gene_data in gene_indices.items():
        gene_id, chrom = gene_chrom.split("|")
        #chrom=gene_data['Chromosome']
        start = gene_data['MinIdx']
        end = gene_data['MaxIdx']
        #append new Feature keyed by chromosome.
        gene_spans[chrom].append(Feature(chromosome=chrom,
                                         gene_id=gene_id,
                                         start=start,
                                         end=end))
            
    #Sort gene_spans for use with bisect functionality
    for chrom, features in gene_spans.items():
        gene_spans[chrom] = sorted(features)

    if pass_back_genes_by_names:
        return chromosomes, gene_spans, gene_by_names

    return chromosomes, gene_spans

