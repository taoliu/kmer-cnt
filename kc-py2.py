#!/usr/bin/env python

import sys
import gzip

base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)

def count_kmer(h, k, seq):
    l = len(seq)
    if l < k: return
    for i in range(l - k + 1):
        kmer_for = seq[i:(i+k)]
        if 'N' in kmer_for: continue
        kmer_rev = kmer_for.translate(comp_tab)[::-1]
        if kmer_for < kmer_rev: kmer = kmer_for
        else: kmer = kmer_rev
        if kmer in h:
            h[kmer] += 1
        else: h[kmer] = 1
    return h

def count_stdin(k):
    counter = {}
    seq = []
    for line in sys.stdin:
        if line[0] == '>':
            if len(seq) > 0:
                count_kmer(counter, k, ''.join(seq))
                seq = []
            else:
                seq.append(line[:-1])
    if len(seq) > 0:
        count_kmer(counter, k, ''.join(seq).upper())
    return counter

def count_gzfile(k, filename):
    counter = {}
    seq = []
    fh = gzip.open(filename, "rt")
    for line in fh:
        if line[0] == '>':
            if len(seq) > 0:
                count_kmer(counter, k, ''.join(seq))
                seq = []
            else:
                seq.append(line[:-1])
    if len(seq) > 0:
        count_kmer(counter, k, ''.join(seq).upper())
    fh.close()
    return counter

def print_hist(counter):
    hist = [0] * 256
    for kmer in counter:
        cnt = counter[kmer]
        if cnt > 255: cnt = 255
        hist[cnt] += 1
    for i in range(1, 256):
        print("{}\t{}".format(i, hist[i]))

def main():
    if len(sys.argv) < 3:
        sys.stderr.write("need 2 paras: %s <n kmer> <gzip filename>\n" % sys.argv[0])
        sys.exit(1)
    nkmer = int(sys.argv[1])
    filename =  str(sys.argv[2])
    counter = count_gzfile( nkmer, filename )
    print_hist(counter)

if __name__ == '__main__':
    main()
