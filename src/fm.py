from __future__ import annotations
import argparse
from dataclasses import dataclass
from parsers import parse_fasta, parse_fastq
import sys
from io import TextIOWrapper, BufferedReader, BufferedWriter
import pickle


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
    args = argparser.parse_args()

    if args.p:
        print(f"Preprocess {args.genome.name}")
        genome = parse_fasta(args.genome)
        processed_genome = {chr: fm_index(genome[chr]) for chr in genome}
        with open(f"{args.genome.name}.bin", "wb") as file:
            pickle.dump(processed_genome, file)

    else:
        # here we need the optional argument reads
        if args.reads is None:
            print("argument: reads is not optional")
            argparser.print_help()
            sys.exit(1)

        try:
            with open(f"{args.genome.name}.bin", "rb") as file:
                genome = pickle.load(file)
        except FileNotFoundError:
            print(
                f"{args.genome.name} has not been preprocessed, you should probably do that... see how here:")
            argparser.print_help()
            sys.exit(1)

        reads = parse_fastq(args.reads)

        print(f"Search {args.genome.name} for {args.reads.name}")
        out = []
        for chr in genome:
            for read in reads:
                hits = genome[chr].search(reads[read])
                for hit in hits:
                    out.append(
                        f'{read}\t{chr}\t{hit+1}\t{len(reads[read])}M\t{reads[read]}')
        out.sort()
        print('\n'.join(out))


@dataclass
class rotating_string:
    x: str

    def __getitem__(self, i: int) -> str:
        return self.x[(i) % len(self.x)]

    def __len__(self) -> int:
        return len(self.x)


class fm_index:
    rx: rotating_string
    sa: list[int]
    o: list[list[int]]
    c: dict[str, int]
    a_map: dict[str, int]

    def __init__(self, x: str) -> None:
        '''
        constructs our fm index from a str
        '''
        self.rx, self.sa, self.o, self.c, self.a_map = pre_process(x)

    def search(self, p: str) -> list[int]:
        if not p or self.rx.x == '$':
            return []
        l, r = 0, len(self.rx)
        for a in p[::-1]:
            if l == r:
                break
            l = self.c[a] + self.o[l][self.a_map[a]]
            r = self.c[a] + self.o[r][self.a_map[a]]
        return self.sa[l:r]

    def __repr__(self) -> str:
        return f"{self.rx.x}\n{self.sa}\n{self.o}\n{self.c}"


def sa_construct(x: str) -> list[int]:
    '''
    Simple SA construction algorithm

    If we get a more efficient version running, this will be updated, but right now we just need a working version.
    '''
    suf = [x[i:] for i, _ in enumerate(x)]
    suf = sorted(suf)
    return [len(x) - len(i) for i in suf]


def fm_preprocess(x: rotating_string, sa: list[int]) -> tuple[list[list[int]], dict[str, int], dict[str, int]]:
    # make c and o table
    alpha = sorted(set(x.x))
    a_map = {c: i for i, c in enumerate(alpha)}

    # C table, leaving it as a dict, so we can
    c = {a: 0 for a in alpha}
    for a in x.x:
        c[a] += 1
    accsum = 0
    for bucket in c:
        c[bucket], accsum = accsum, accsum + c[bucket]

    # should probably stop using dict for this, and just map them down to indices in a list...
    O = [[0 for _ in alpha] for _ in range(len(x)+1)]
    for i in range(1, len(x)+1):
        for a in a_map:
            O[i][a_map[a]] = O[i-1][a_map[a]] + \
                (a_map[a] == a_map[x[len(x) + sa[i-1]-1]])

    return O, c, a_map


def pre_process(x: str):
    sa = sa_construct(x+'$')
    rx = rotating_string(x+'$')
    o, c, a_map = fm_preprocess(rx, sa)
    return rx, sa, o, c, a_map


def dump_fm(index: fm_index, file: BufferedWriter) -> None:
    pickle.dump(obj=index, file=file)
    return None


def load_fm(file: BufferedReader) -> fm_index:
    return pickle.load(file)


# probably deprecated by now.
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


if __name__ == '__main__':
    main()
    #x = "mississippi"
    #fm = fm_index(x)
    # print(repr(fm))
    # with open("test.bin", "wb") as file:
    #    dump_fm(fm, file)
    # with open("test.bin", "rb") as file:
    #    loaded = load_fm(file)
    # print(repr(loaded))
