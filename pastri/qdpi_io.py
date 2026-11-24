import gzip
import heapq
from contextlib import ExitStack
from itertools import groupby
from operator import itemgetter

import numpy as np

EMPTY = np.array([], dtype=int)

def to_int(data):
    """
    String to numpy
    """
    if data == '.':
        return EMPTY
    return np.fromstring(data, dtype=int, sep=",")


def _iter_loci(fh, file_idx, subset=None):
    """
    Lazily parse one qdpi file.
    Yields (sort_key, file_idx, locus, fields).
    """
    for raw in fh:
        line = raw.decode().rstrip('\n')
        if not line:
            continue
        data = line.split('\t')
        locus = tuple(data[:3])
        if subset and locus not in subset:
            continue
        yield (data[0], int(data[1]), int(data[2])), file_idx, locus, data


def stream_qdpi(in_fhs, subset=None, with_unphased=False):
    """
    Given a set of qdpi files, stream/transform/join.
    Yields per-locus (chrom,start,end) and h1/h2 lists of arrays.
    Each list has one entry per input file, in input order. Files that
    don't contain the locus get `missing` in their slot.
    """
    with ExitStack() as stack:
        files = [stack.enter_context(gzip.open(f)) for f in in_fhs]
        n = len(files)
        streams = [_iter_loci(fh, i, subset) for i, fh in enumerate(files)]
        merged = heapq.merge(*streams, key=itemgetter(0))

        for _, group in groupby(merged, key=itemgetter(0)):
            h0_parts = [EMPTY] * n
            h1_parts = [EMPTY] * n
            h2_parts = [EMPTY] * n
            locus = None
            for _key, i, locus, data in group:
                if with_unphased:
                    h0_parts[i] = to_int(data[4])
                h1_parts[i] = to_int(data[5])
                h2_parts[i] = to_int(data[6])

            if with_unphased:
                yield locus, h0_parts, h1_parts, h2_parts
            else:
                yield locus, h1_parts, h2_parts
