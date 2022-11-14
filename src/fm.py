from __future__ import annotations
import argparse
from dataclasses import dataclass
import sys


def main():
    argparser = argparse.ArgumentParser(
        description="FM-index exact pattern matching",
        usage="\n\tfm -p genome\n\tfm genome reads"
    )
    argparser.add_argument(
        "-p", action="store_true",
        help="preprocess the genome."
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
    args = argparser.parse_args()

    if args.p:
        print(f"Preprocess {args.genome}")
    else:
        # here we need the optional argument reads
        if args.reads is None:
            argparser.print_help()
            sys.exit(1)
        print(f"Search {args.genome} for {args.reads}")


@dataclass
class rotating_string:
    x: str

    def __getitem__(self, i) -> str:
        return self.x[(i) % len(self.x)]

    def _len__(self) -> int:
        return len(self.x)


def sa_construct(x: str) -> list[int]:
    ''' 
    Simple SA construction algorithm 

    If we get a more efficient version running, this will be updated, but right now we just need a working version.
    '''
    suf = [x[i:] for i, _ in enumerate(x)]
    suf = sorted(suf)
    return [len(x) - len(i) for i in suf]


def bwt(x: str) -> str:
    "returns the last column of a BWT from string x"
    sa = sa_construct(x)
    x = rotating_string(x)

    l = [x[-1-i] for i in sa]
    return ''.join(l)


def compress(x: str) -> str:
    ''' compresses string x using run length encoding '''
    rle = []
    run = 1
    for i, c in enumerate(x):

        if c == x[i+1]


if __name__ == '__main__':
    main()
