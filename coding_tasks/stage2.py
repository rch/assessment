import os, argparse
from collections import Counter
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description='find 10 most frequent sequences')
    parser.add_argument('fasta', type=str, nargs=1, help='fasta file')
    args = parser.parse_args()
    
    index = Counter()
    with open(args.fasta[0]) as f:
        for entry in SeqIO.parse(f, "fasta"):
            index[str(entry.seq)] += 1
            
    with open('most_common.txt','w') as f:
        for seq, num in index.most_common(10):
            f.write('> {}\n{}\n'.format(num, seq))


if __name__ == '__main__':
    main()
