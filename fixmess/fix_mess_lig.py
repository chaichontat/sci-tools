# %%
from typing import Any, Callable, Collection, Sequence

import numpy as np
import polars as pl
from Levenshtein import distance

from scitools.bases import min_distance
from scitools.construct.generator import gen_lig
from scitools.construct.sci import generate_n

pl.Config.set_fmt_str_lengths(65)


lig = pl.read_csv(
    "scirnaseq/lig_plate.tsv",
    separator="\t",
    has_header=False,
    new_columns=["well", "name", "seq"],
).with_columns(
    col=pl.col("well").str.slice(1, 23).cast(pl.UInt16),
    row=pl.col("well").str.slice(0, 1),
    seq=pl.col("seq").str.replace_all(" ", ""),
)

existing = lig.filter(pl.col("col").apply(lambda x: np.bitwise_and(x, 1)) != 0).with_columns(
    idx=pl.col("seq").str.slice(7, 10)
)

dump = lig.filter(pl.col("col").apply(lambda x: np.bitwise_and(x, 1)) == 0).with_columns(
    idx=pl.col("seq").str.slice(7, 10)
)


# %%
new = generate_n(ord("l"), out=[(x, None) for x in existing["idx"]], max_gen=96)
assert len(set(existing["idx"]) & set(new[48:])) == 0

new_df = (
    dump.clone()
    .replace("idx", pl.Series(new[48:]))
    .with_columns(seq=pl.col("idx").apply(gen_lig), name="sciv3-lig-" + pl.col("well"))
)
new_df.write_csv("fixmess/sciv3-lig.csv")

# %%
combi = pl.concat([existing, new_df])
assert min_distance(combi["idx"], distance) >= 4
assert combi["seq"].apply(len).unique().len() == 1

# %%
