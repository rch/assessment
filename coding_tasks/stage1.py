import os, argparse
from Bio import SeqIO


def fastq_files(top):
    for entry in os.walk(top):
        for filename in entry[2]:
            if os.path.splitext(filename)[1] == '.fastq':
                yield os.path.join(entry[0], filename)
            

def main():
    parser = argparse.ArgumentParser(description='reports percent of sequences greater than 30 nt long')
    parser.add_argument('directory', type=str, nargs='?', default='.', help='top-level directory')
    args = parser.parse_args()
    
    with open('percent_long.txt','w') as output:
        for filename in fastq_files(args.directory):
            with open(filename) as f:
                i = 0
                for num, entry in enumerate(SeqIO.parse(f, "fastq")):
                    if len(entry.seq) > 30:
                        i += 1
                output.write('{}\t{:.4f}%\n'.format(filename, float(i)/num * 100)) 


if __name__ == '__main__':
    main()
