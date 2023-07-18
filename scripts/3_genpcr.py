# %%
import random
from typing import Any, Callable

import polars as pl
from Levenshtein import distance, hamming
from nupack import Model, SetSpec, Strand, Tube, tube_analysis

from scitools.bases import gen_sample_space, reverse_complement
from scitools.constants import ADT_LINKER, NEXTERA2, OVERHANG, P5, P7, SMALLRNA, TRUSEQ1
from scitools.sci import gen_set, hairpin
from scitools.utils import repeat

my_model = Model(material="dna", celsius=63, sodium=0.3, magnesium=0.003)

existing_p5 = pl.read_csv("scripts/sci-v2-p5.tsv", separator="\t", has_header=False)["column_1"]
ok = []
for seq in existing_p5:
    for x in ok:
        if distance(x, seq) < 4:
            break
    else:
        ok.append(seq)
print(len(ok))
# %%


def gen_p7(idx: str) -> str:
    return P7 + idx.lower() + NEXTERA2


def gen_p7adt(idx: str) -> str:
    return ADT_LINKER + idx.lower() + SMALLRNA[:19]


def gen_p5(idx: str):
    return P5 + idx.lower() + TRUSEQ1[:21]


def gen_lig(idx: str):
    # T is presumably a buffer to prevent digestion
    return (
        reverse_complement(OVERHANG)
        + idx.lower()
        + "T"
        + "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
        + reverse_complement(idx).lower()
    )


def gen_rt(idx: str):
    return "5/Phos/" + OVERHANG + idx.lower() + "N" * 9 + "T" * 28 + "VN"


# %%

p7_screen = dict(
    nextera=lambda seq: ((mfe := hairpin(gen_p7(seq), my_model)), mfe >= -0.06),
    smallrna=lambda seq: ((mfe := hairpin(gen_p7adt(seq), my_model)), mfe >= -0.06),
)

p5_screen = dict(p5=lambda seq: ((mfe := hairpin(gen_p5(seq), my_model)), mfe >= -0.37))


def gen_idxs(seed: int, k: int = 10, screens: dict[str, Callable[[str], tuple[Any, bool]]] = {}, **kwargs):
    base = gen_sample_space(repeat("N", k))
    random.seed(seed)
    random.shuffle(base)
    return [x[0] for x in gen_set(base, "".join(random.choices("ATGC", k=k)), screens=screens, **kwargs)]


df = pl.DataFrame(
    {
        # "rt": gen_idxs(0, no_gg_start=False),
        # "lig": gen_idxs(1, no_four_gc=False, gc_range=(0.4, 0.7)),
        # "p7adt": gen_idxs(2, screens={"p7": p7_screen["smallrna"]}),
        # "p7": gen_idxs(3, screens={"p7": p7_screen["nextera"]}),
        "p5": gen_idxs(4, screens={"p5": p5_screen["p5"]}, out=[(x, None) for x in ok]),
    }
)

# %%
df.write_csv("sciv3idx.csv")

# %%
seqed = df.with_columns(
    rt=pl.col("rt").apply(gen_rt),
    lig=pl.col("lig").apply(gen_lig),
    p7adt=pl.col("p7adt").apply(gen_p7adt),
    p7=pl.col("p7").apply(gen_p7),
    p5=pl.col("p5").apply(gen_p5),
)
# %%
seqed.write_csv("sciv3.csv")
# %%
