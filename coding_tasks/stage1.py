import os, argparse
from itertools import count
from fractions import Fraction
from decimal import Decimal
from Bio import SeqIO

def fastq_files(top):
    for entry in os.walk(top):
        for filename in entry[2]:
            if os.path.splitext(filename)[1] == '.fastq':
                yield os.path.join(entry[0], filename)
            

def parse(handle): 
        return SeqIO.parse(handle, "fastq")


def main():
    parser = argparse.ArgumentParser(description='reports percent of sequences greater than 30 nt long')
    parser.add_argument('directory', type=str, nargs='?', default='.', help='top-level directory')
    args = parser.parse_args()
    for filename in fastq_files(args.directory):
        with open(filename) as f:
            i = count()
            for num, entry in enumerate(parse(f)):
                if len(entry.seq) > 30:
                    i.next()
            print filename, "{:.4f}%".format(float(i.next())/num * 100) 


if __name__ == '__main__':
    main()
