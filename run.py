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

RUST_PATH = Path("/home/chaichontat/rawSeq/")
if not (RUST_PATH / "splitter").exists() and (RUST_PATH / "remap").exists():
    raise RuntimeError

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
    "i1", type=click.Path(exists=True, file_okay=True, dir_okay=False, path_type=Path)
)
@click.option(
    "--idxs",
    type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=Path),
    default=Path.cwd(),
)
@click.option("--rt", type=int, default=1)
@click.option("--lig", type=int, default=1)
@click.option("--nokallisto", is_flag=True)
def run(i1: Path, idxs: Path, rt: int = 1, lig: int = 1, nokallisto: bool = False):
    PARAMS = {
        "rt": Param("finalRT.tsv", "R1", 23, rt),
        "p5": Param("finalP5.tsv", "I2", 0, 1),
        "p7": Param("finalP7.tsv", "I1", 0, 1),
        "lig": Param("finalLig.tsv", "R1", 0, lig),
    }
    path = i1.parent
    lane = re.search(r"Undetermined_S0_(L[0-9]+).*", i1.name).group(1)

    done_files = [
        path / f"{lane}_1idxpcr.DONE",
        path / f"{lane}_2remap.DONE",
        path / f"{lane}_3rt.DONE",
        path / f"{lane}_4fakeidx.DONE",
        path / f"{lane}_5kallisto.DONE",
    ]

    # with check_done(path / f"{lane}_1idxpcr.DONE"):
    if not (done_files[0]).exists():
        with ThreadPoolExecutor() as exc:
            for k, v in PARAMS.items():
                if k == "rt":
                    continue
                name = re.sub(r"(I1|I2|R1|R2)", v.read, i1.name)
                exc.submit(
                    subprocess.run,
                    f"pigz -dc {path/name} | {RUST_PATH}splitter - {idxs.resolve()}/{v.file} --index {v.index} --tol {v.tol} | pigz > {path}/{lane}_{k}.tsv.gz",
                    shell=True,
                )
        done_files[0].touch()

    k, v = "rt", PARAMS["rt"]
    name = re.sub(r"(I1|I2|R1|R2)", v.read, i1.name)

    if not (done_files[1]).exists():
        logger.info(f"Remapping. Reading from {lane}_lig.tsv.gz.")
        subprocess.run(
            f"pigz -dc {path}/{lane}_lig.tsv.gz| {RUST_PATH}remap {idxs}/shortlong.tsv | pigz > {path}/{lane}_ligmode.tsv.gz",
            shell=True,
            env=os.environ,
            check=True,
        )
        done_files[1].touch()

    if not (done_files[2]).exists():
        subprocess.run(
            f'fish -c "conda activate seq && pigz -dc {path/name} | {RUST_PATH}splitter - {idxs.resolve()}/{v.file} --index {v.index} --tol {v.tol}'
            f' -s (pigz -dc {path}/{lane}_ligmode.tsv.gz | psub) | pigz > {path}/{lane}_rt.tsv.gz"',
            shell=True,
            env=os.environ,
            check=True,
        )
        done_files[2].touch()

    if not (done_files[3]).exists():
        subprocess.run(
            f"python makefakeindex.py --idxs {idxs} --lane {int(lane[1:])} {path/name} | pigz > {path}/{lane}_combi.fastq.gz",
            shell=True,
            check=True,
        )
        done_files[3].touch()

    if not (done_files[4]).exists() and not nokallisto:
        indices = " ".join([f"-{k} {v}" for k, v in kfiles.items()])

        subprocess.run(
            f'kb count --h5ad -x \'0,0,28:0,28,37:1,0,0\' -t 16 --workflow nac --strand forward'
            f' {indices}'
            f' {path}/{lane}_combi.fastq.gz {path/re.sub(r"(I1|I2|R1|R2)", "R2", i1.name)} -o {path/lane} --overwrite --verbose'
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
