# %%
import gzip
import sys
from itertools import product
from pathlib import Path

import click
import polars as pl
import pyfastx


def extract_idx(file: Path, length_limit: int):
    df = pl.read_csv(
        file, has_header=False, new_columns=["idx", "well"], separator="\t"
    )
    return dict(zip(df["well"], df["idx"].str.slice(0, length=length_limit)))


LENGTH_LIMIT = {"p5": 6, "p7": 6, "lig": 8, "rt": 8}


# %%
@click.command()
@click.argument("input", type=click.Path(exists=True))
@click.option(
    "--idxs",
    type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=Path),
    default=Path.cwd(),
)
@click.option("--lane", type=int)
def run(input: str, lane: int, idxs: Path = Path.cwd()):
    p5, p7, lig, rt = [
        extract_idx(idxs / name, lim)
        for name, lim in zip(
            ["finalP5.tsv", "finalP7.tsv", "finalLig.tsv", "finalRT.tsv"],
            [8, 8, 6, 6],
        )
    ]

    Is = "I" * 28
    Ns = "N" * 28
    fs = zip(*[
        gzip.open(Path(input).parent / f"L{lane:03d}_{name}.tsv.gz", "rt")
        for name in ["rt", "lig", "p7", "p5"]
    ])

    for i, ((name, seq, qual), (_rt, _lig, _P7, _P5)) in enumerate(
        zip(
            pyfastx.Fastq(
                input,
                build_index=False,
                full_name=True,
            ),
            fs,
        )
    ):
        umi, umiqual = seq[15:24], qual[15:24]
        if any("unk" in x for x in [_rt, _lig, _P7, _P5]):
            sys.stdout.write(f"@{name}\n{Ns}{umi}\n+\n{Is}{umiqual}\n\n")
        else:
            well, qual = zip(*[y.strip().split("\t") for y in [_rt, _lig, _P5, _P7]])
            qual = [qual[0][:6], qual[1][:6], qual[2][:8], qual[3][:8]]

            seq = [x[y] for x, y in zip([rt, lig, p5, p7], well)]
            seq = "".join(seq)

            sys.stdout.write(f"@{name}\n{seq}{umi}\n+\n{''.join(qual)}{umiqual}\n\n")

        if i % 1000000 == 0:
            sys.stderr.write(f"{i}\n")


if __name__ == "__main__":
    run()
# %%


def gen_whitelist(idxs: Path, output: Path):
    p5, p7, lig, rt = [
        extract_idx(idxs / name, lim).values()
        for name, lim in zip(
            ["finalP5.tsv", "finalP7.tsv", "finalLig.tsv", "finalRT.tsv"],
            [8, 8, 6, 6],
        )
    ]
    p7p5 = [x + y for x, y in zip(p7, p5)]
    whitelist = ["".join(combi) for combi in product(rt, lig, p7p5)]
    with open(output, "w") as f:
        f.write("\n".join(whitelist))


# %%
# gen_whitelist(Path("/fast/2024"), Path("/fast/2024/whitelist.txt"))
# %%
