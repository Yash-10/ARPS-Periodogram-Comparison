import random
import lightkurve
from lightkurve import LightCurve


def rand_parts(seq, n, l):
    """Extracts `n` random contiguous sub-sequences from `seq`, each of length `l`.

    Args:
        seq (list): List containing the elements.
        n (int): No. of subsequences to extract.
        l (int): No. of elements in each subsequence.

    Returns:
        list: The desired result.

    """
    random.seed(1)
    n = int(n)
    l = int(l)
    indices = range(int(len(seq) - (l - 1) * n))
    result = []
    offset = 0
    for i in sorted(random.sample(indices, int(n))):
        i += offset
        result.append(seq[i:i+l])
        offset += l - 1
    return result


def fold_lc(lc, t, period):
    Lc = LightCurve(flux=lc, time=t)
    folded = Lc.fold(period=period, epoch_time=0, epoch_phase=0)
    return folded
