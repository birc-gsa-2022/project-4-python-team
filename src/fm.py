from __future__ import annotations
import argparse
from dataclasses import dataclass
import sys
from io import TextIOWrapper


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

    def __getitem__(self, i: int) -> str:
        return self.x[(i) % len(self.x)]

    def __len__(self) -> int:
        return len(self.x)


def sa_construct(x: str) -> list[int]:
    '''
    Simple SA construction algorithm

    If we get a more efficient version running, this will be updated, but right now we just need a working version.
    '''
    suf = [x[i:] for i, _ in enumerate(x)]
    suf = sorted(suf)
    return [len(x) - len(i) for i in suf]


def fm_search(x: rotating_string, p: str, sa: list[int]) -> list[int]:
    # make c and o table (this should be moved to pre processing)
    alpha = sorted(set(x.x))
    c = {a: 0 for a in alpha}
    for a in x.x:
        c[a] += 1
    accsum = 0
    for bucket in c:
        c[bucket], accsum = accsum, accsum + c[bucket]

    # should probably stop using dict for this, and just map them down to indices in a list...
    O = [{a: 0 for a in c} for _ in range(len(x)+1)]
    for i, o in enumerate(O):
        if i != 0 and i < len(x):
            o[x[sa[i]]] = O[i-1][x[sa[i]]] + (x[sa[i]] == x[sa[i-1]])

    print('\n'.join([str(foo) for foo in O]))

    l, r = 0, len(x)
    for a in p[::-1]:
        if l == r:
            break


def bwt(x: str) -> str:
    "returns the last column of a BWT from string x"
    x += "$"
    sa = sa_construct(x)
    rx = rotating_string(x)

    l = [rx[i-1] for i in sa]
    return ''.join(l)


def rle(x: str) -> str:
    ''' compresses string x using run length encoding '''
    rle = ''
    run = 1
    end = len(x)
    i = 1
    while i < end:
        if x[i] == x[i-1]:
            run += 1
            i += 1
        else:
            rle += str(run) + x[i-1]
            run = 1
            i += 1
    else:
        rle += str(run)+x[i-1]
    return rle


def compress(x: str) -> str:
    l = bwt(x)
    return rle(l)


def decompress(x: str) -> str:
    decom = ''
    for i, c in zip(x[0::2], x[1::2]):
        decom += int(i) * c
    return decom


def ith_occurence(i: int, c: str, x: str) -> int:
    for j, char in enumerate(x):
        if char == c:
            i -= 1
            if i < 0:
                return j


def reverse_bwt(l: str) -> str:
    # yes, we know that we can radix the first alphabet by bit represenation, but hey this is python.
    alpha = sorted(set(l))
    buckets = {c: 0 for c in alpha}
    f = [str()]*len(l)

    # bucket sort the l column
    for c in l:
        buckets[c] += 1
    accsum = 0
    for bucket in buckets:
        buckets[bucket], accsum = accsum, accsum + buckets[bucket]
    for c in l:
        f[buckets[c]] = c
        buckets[c] += 1

    # use buckets to repreesent ranks
    rank = {c: 0 for c in alpha}
    ranked_f: list[tuple[str, int]] = []
    for char in f:
        ranked_f.append((char, rank[char]))
        rank[char] += 1

    out = ""
    # find sentinel
    first = ith_occurence(0, '$', l)
    out += ranked_f[first][0]
    rank = ranked_f[first][1]

    i = 1
    while i < len(l):
        next = ith_occurence(rank, out[i-1], l)
        out += ranked_f[next][0]
        rank = ranked_f[next][1]
        i += 1
    return out


# For pre processing, i figured it might be worthwile to save the SA given, that it makes the FM index lookups easier.
def pre_process(x: str, outfile: str):
    t = compress(x)
    sa = sa_construct(x+"$")

    with open(outfile, "w") as file:
        file.write(t+'\n')
        file.write('\t'.join([str(i) for i in sa]))


def read_preprocessed_genome(infile: TextIOWrapper) -> tuple[str, list[int], rotating_string]:
    t = decompress(infile.readline())
    sa = infile.readline().split("\t")
    sa = [int(i) for i in sa]
    genome = rotating_string(reverse_bwt(t))
    return t, sa, genome


if __name__ == '__main__':
    # main()
    x = "mississippi"
    sa = sa_construct(x+"$")
    print(sa)
    t = bwt(x)
    print(t)
    print(rle(t))
    print(reverse_bwt(t))

    comp = compress(t)
    print(comp)
    print(decompress(comp))
    fm_search(rotating_string(x+"$"), "isi", sa)
