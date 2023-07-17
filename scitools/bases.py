from collections import Counter
from itertools import product
from typing import Callable

import primer3
from Levenshtein import distance

_revtrans = str.maketrans("ATCGatcg", "TAGCtagc")


def reverse_complement(seq: str) -> str:
    return seq[::-1].translate(_revtrans)


def calculate_base_diversity(seq: list[str]) -> list[dict[str, float]]:
    length = len(seq[0])
    if not all(len(x) == length for x in seq):
        raise ValueError("Sequences must be the same length")
    return [{k: v / len(seq) for k, v in Counter([x[i] for x in seq]).items()} for i in range(length)]


def min_distance(seqs: list[str], dist: Callable[[str, str], int] = distance) -> float:
    out = float("inf")
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            out = min(out, dist(seqs[i], seqs[j]))
    return out


def gen_sample_space(template: str) -> list[str]:
    mapping = {
        "A": "A",
        "T": "T",
        "G": "G",
        "C": "C",
        "N": "ATGC",
        "W": "AT",
        "S": "GC",
        "M": "AC",
        "K": "GT",
        "R": "AG",
        "Y": "CT",
        "B": "CGT",
        "D": "AGT",
        "H": "ACT",
        "V": "ACG",
    }
    iters = [mapping[x] for x in template]
    base = ["".join(x) for x in product(*iters)]
    return base


def gc_content(seq: str):
    return (seq.count("G") + seq.count("C") + seq.count("g") + seq.count("c")) / len(seq)


def screen_umi_interactions(
    seqs: list[str],
    umis: list[str],
    conditions: dict[str, float],
    f: Callable[[str, str], str] = lambda x, y: x,
) -> list[list[float]]:
    idxs = []
    for seq in seqs:
        o = [primer3.calc_hairpin_tm(f(seq, umi), **conditions) for umi in umis]
        idxs.append(o)
    return idxs
