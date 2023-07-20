# %%
from typing import Callable, Sequence

import numpy as np
import polars as pl
from Levenshtein import distance

from scitools.bases import min_distance, wells_iter
from scitools.construct.generator import gen_p5, gen_p7, gen_p7adt, gen_rt
from scitools.construct.sci import generate_n

mapping = {"A": 0, "C": 1, "G": 2, "T": 3}


def extract_idx(df: pl.DataFrame, start: int, length: int):
    return df["column_3"].str.replace_all(" ", "").str.slice(start, length).to_list()


def max_retention(data: Sequence[str], start: int, thresh: int = 4):
    out = [start]
    for i in range(len(data)):
        if i == start:
            continue
        for x in out:
            if distance(data[x], data[i]) < thresh:
                break
        else:
            out.append(i)
    return out


def gen_fix(df: pl.DataFrame, f: Callable):
    data = np.array(df["idx"].to_list())
    counts = {start: len(max_retention(data, start)) for start in range(len(data))}
    best = max(counts, key=counts.get)
    print(best, counts[best])

    idxs_good = sorted(max_retention(data, best))
    idxs_bad = sorted(set(range(96)) - set(idxs_good))

    new = generate_n(0, out=[(x, None) for x in np.array(data)[idxs_good]], max_gen=96)
    assert len(set(data[idxs_good]) & set(new[len(idxs_good) :])) == 0
    new_df = (
        df[idxs_bad]
        .clone()
        .replace("idx", pl.Series(new[len(idxs_good) :]))
        .with_columns(seq=pl.col("idx").apply(f))
    )
    assert min_distance((df[idxs_good]["idx"].to_list() + new_df["idx"].to_list()), distance) >= 4

    return new_df


p7 = pl.read_csv("scirnaseq/sci-v2-p7.tsv", separator="\t", has_header=False, new_columns=["idx", "well"])

# %%
new_p7 = generate_n(9, out=[(x, None) for x in p7["idx"]], max_gen=192)
p7_adt = pl.DataFrame(dict(idx=new_p7[96:], well=list(wells_iter()))).with_columns(
    seq=pl.col("idx").apply(gen_p7adt),
    name="sciv3-P7adt-" + pl.col("well"),
)
p7_adt.write_csv("fixmess/sciv3-P7adt.csv")

# P7
wells = [f"{r}{c:02d}" for r in "ABCDEFGH" for c in range(4, 13)]
new_p7 = pl.DataFrame(dict(idx=new_p7[24:96], well=wells)).with_columns(
    seq=pl.col("idx").apply(gen_p7), name="sciv3-P7-" + pl.col("well")
)
# Move to P5 to get that to 24.
new_p7[:-2].write_csv("fixmess/sciv3-P7.csv")

# %%
p5 = pl.read_csv("scirnaseq/sci-v2-p5.tsv", separator="\t", has_header=False, new_columns=["idx", "well"])
new_p5 = gen_fix(p5, gen_p5).with_columns(name="sciv3-P5-" + pl.col("well"))

pl.concat([new_p5, new_p7[-2:]]).write_csv("fixmess/sciv3-P5.csv")

# %%
new_rt = (
    pl.DataFrame({"idx": (generate_n(ord("R"), max_gen=96))})
    .with_columns(
        seq=pl.col("idx").apply(gen_rt),
        well=pl.lit(list(wells_iter())),  # type: ignore
    )
    .with_columns(name="sciv3-RT-" + pl.col("well"))
)
new_rt.write_csv("fixmess/sciv3-RT.csv")


# %%
