import os, argparse
from collections import Counter
from Bio import SeqIO


def parse(handle): 
        return SeqIO.parse(handle, "fasta")


def main():
    parser = argparse.ArgumentParser(description='reports percent of sequences greater than 30 nt long')
    parser.add_argument('fasta', type=str, nargs=1, help='fasta file')
    args = parser.parse_args()
    index = Counter()
    with open(args.fasta[0]) as f:
        for num, entry in enumerate(parse(f)):
            index[str(entry.seq)] += 1
    for seq, num in index.most_common(10):
        print '>', num
        print seq

if __name__ == '__main__':
    main()
