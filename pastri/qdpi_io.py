import gzip
import numpy as np

EMPTY = np.array([], dtype=int)

def to_int(data):
    """
    String to numpy
    """
    if data == '.':
        return EMPTY
    return np.fromstring(data, dtype=int, sep=",")

def stream_qdpi(in_fhs, subset=None):
    """
    Given a set of qdpi files, stream/transform/join
    Yields per-locus (chrom,start,end) and h1/h2 lists of arrays
    """
    files = [gzip.open(_) for _ in in_fhs]
    while True:
        # Join each qdpi bed file locus
        try:
            lines = [next(_).decode() for _ in files]
        except StopIteration:
            break

        locus = None
        h1_parts = []
        h2_parts = []
        for line in lines:
            data = line.strip().split('\t')
            if subset and tuple(data[:3]) not in subset:
                break
            locus = data[:3]
            h1_parts.append(to_int(data[5]))
            h2_parts.append(to_int(data[6]))
        if locus is not None:
            yield tuple(locus), h1_parts, h2_parts

