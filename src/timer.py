import time
import argparse
import numpy as np
from parsers import parse_fasta, parse_fastq
from fm import fm_index
import csv

def main():
    argparser = argparse.ArgumentParser(
        description="FM-index exact pattern matching",
        usage="\n\tfm -p genome\n\tfm genome reads"
    )
    argparser.add_argument(
        "-p", action="store_true",
        help="preprocess the genome. output is written to [genome].bin"
    )
    argparser.add_argument(
        "genome",
        help="Simple-FASTA file containing the genome.",
        type=argparse.FileType('r')
    )
    argparser.add_argument(
        "reads", nargs="?",
        help="Simple-FASTQ file containing the reads.",
        type=argparse.FileType('r')
    )
    argparser.add_argument(
        "output", nargs="?",
        type=str
    )
    args = argparser.parse_args()

    genome=parse_fasta(args.genome)
    reads=parse_fastq(args.reads)

    time1, time2 =time_fm(genome, reads)
    with open("{}genome.csv".format(args.output),"w") as f:
        write=csv.writer(f)
        write.writerow(time1)
    with open("{}reads.csv".format(args.output),"w") as f:
        write=csv.writer(f)
        write.writerow(time2)

def time_fm(genome, reads):
    out1=[]
    out2=[]
    for chr in genome:
        print(chr)
        t=time.process_time()
        current=fm_index(genome[chr])
        out1.append(time.process_time()-t)
        t=time.process_time()
        for read in reads:
            current.search(reads[read])
        out2.append((time.process_time()-t)/len(reads))
    return out1, out2

if __name__ == '__main__':
    main()
