import os, argparse
from collections import defaultdict, namedtuple
from operator import itemgetter
from Bio import SeqIO

Entry = namedtuple('Entry',['seqname',
                            'source',
                            'feature',
                            'start',
                            'end',
                            'score',
                            'strand',
                            'frame',
                            'attribute'])
                            
def parse(handle): 
    """ Parse GTF input
        seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix.
        source - name of the program that generated this feature, or the data source (database or project name)
        feature - feature type name, e.g. Gene, Variation, Similarity
        start - Start position of the feature, with sequence numbering starting at 1.
        end - End position of the feature, with sequence numbering starting at 1.
        score - A floating point value.
        strand - defined as + (forward) or - (reverse).
        frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
        attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.
    """
    for line in handle:
        vals = [int(v) if v.isdigit() else v for v in line.strip().split('\t')]
        yield Entry(*vals)

def main():
    parser = argparse.ArgumentParser(description='reports percent of sequences greater than 30 nt long')
    parser.add_argument('positions', type=str, nargs=1, help='positions file')
    parser.add_argument('annotations', type=str, nargs=1, help='annotations file')
    parser.add_argument('output', type=str, nargs=1, help='output file')
    args = parser.parse_args()
    
    positions = defaultdict(list)
    with open(args.positions[0]) as f:
        for num, line in enumerate(f):
            try:
                chrom, position = line.strip().split('\t')
            except ValueError:
                print line
            positions[chrom].append(int(position))
    
    annotations = defaultdict(list)
    with open(args.annotations[0]) as f:
        for entry in parse(f):
            annotations[entry.seqname].append(entry)
    
    results = defaultdict(list)
    for chrom, entries in positions.iteritems():
        entries.sort()
        annotations[chrom].sort(key=itemgetter(Entry._fields.index('start')))
        for pos in entries:
            for num, annotation in enumerate(annotations[chrom]):
                if annotation.end >= pos and annotation.start <= pos:
                    results[(chrom, pos)].append(annotation)
    
    with open(args.output[0],'w') as f:
        for key, entries in results.iteritems():
            print key[0], key[1], ','.join(set([v.attribute for v in entries]))
            f.write('{}\n'.format(key[0], key[1], ','.join(set([v.attribute for v in entries]))))

if __name__ == '__main__':
    main()
    print('\a')    
