# %%
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from functools import wraps
from pathlib import Path
from typing import Literal

import click
from loguru import logger

RUST_PATH = ""
# if not (RUST_PATH / "splitter").exists() and (RUST_PATH / "remap").exists():
#     raise RuntimeError

INDEX_PATH = Path("/home/chaichontat/working/mousenuc/")

kfiles = {
    k: INDEX_PATH / v
    for k, v in dict(
        i="index.idx", g="t2g.txt", c1="cdna.txt", c2="nascent.txt"
    ).items()
}
assert all(v.exists() for v in kfiles.values())


@dataclass
class Param:
    file: str
    read: Literal["I1", "I2", "R1", "R2"]
    index: int
    tol: int


def check_done(f):
    @wraps(f)
    def inner(path: Path):
        if path.exists():
            return
        f(path)
        path.touch()

    return inner


@click.command()
@click.argument(
    "path",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path),
)
@click.argument("lanes", type=str)
@click.option(
    "--idxs",
    type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=Path),
    default=Path.cwd(),
)
@click.option("--output", type=str)
@click.option("--threads", type=int, default=16)
def run(path: Path, lanes: str, output: str, idxs: Path, threads: int = 16):
    lanes = list(map(int, lanes.split(",")))

    indices = " ".join([f"-{k} {v}" for k, v in kfiles.items()])
    files = [
        f"{path}/L00{i}_combi.fastq.gz {path}/Undetermined_S0_L00{i}_R2_001.fastq.gz"
        for i in lanes
    ]

    subprocess.run(
        f'kb count --h5ad -x \'0,0,28:0,28,37:1,0,0\' -t {threads} --workflow nac --strand forward'
        f' {indices}'
        f' {" ".join(files)} -o {path}/{output} --overwrite --verbose'
        f' -w {idxs}/whitelist.txt',
        shell=True,
        check=True,
    )


if __name__ == "__main__":
    run()
# %%
# import polars as pl

# lig = pl.read_csv("finalLig.tsv", separator="\t", has_header=False)

# mapping = dict(
#     zip(
#         [x for x in lig["column_2"]],
#         [x for x in Path("shortlong.tsv").read_text().splitlines()],
#     )
# )

# with open(path) as f:
#     with open(out_path, "w") as g:
#         for i, line in enumerate(f):
#             g.write(mapping[line])
#             if i % 1000000 == 0:
#                 print(i)


# # %%

# %%
