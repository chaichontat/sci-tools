# %%
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import click
from loguru import logger


@dataclass
class Param:
    file: str
    read: Literal["I1", "I2", "R1", "R2"]
    index: int
    tol: int


@click.command()
@click.argument("i1", type=click.Path(exists=True, file_okay=True, dir_okay=False, path_type=Path))
@click.option(
    "--idxs",
    type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=Path),
    default=Path.cwd(),
    description="Path to barcode files.",
)
@click.option(
    "--rustpath",
    type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=Path),
    description="Path to rust executables.",
)
def run(i1: Path, idxs: Path, rustpath: Path):
    PARAMS = {
        "rt": Param("finalRT.tsv", "R1", 23, 1),
        "p5": Param("finalP5.tsv", "I2", 0, 1),
        "p7": Param("finalP7.tsv", "I1", 0, 1),
        "lig": Param("finalLig.tsv", "R1", 0, 1),
    }
    path = i1.parent
    lane = re.search(r"(L[0-9]+)", i1.name).group(1)

    if not (rustpath / "splitter").exists() and (rustpath / "remap").exists():
        raise FileNotFoundError("Rust executables not found. Please download from github.com/chaichontat/")

    done_files = [
        path / f"{lane}_1idxpcr.DONE",
        path / f"{lane}_2remap.DONE",
        path / f"{lane}_3rt.DONE",
        path / f"{lane}_4fakeidx.DONE",
    ]

    if not (done_files[0]).exists():
        with ThreadPoolExecutor() as exc:
            for k, v in PARAMS.items():
                if k == "rt":
                    continue
                name = re.sub(r"(I1|I2|R1|R2)", v.read, i1.name)
                exc.submit(
                    subprocess.run,
                    f"pigz -dc {path/name} | {rustpath/'splitter'} - {idxs.resolve()}/{v.file} --index {v.index} --tol {v.tol} | pigz > {path}/{lane}_{k}.tsv.gz",
                    shell=True,
                )
        done_files[0].touch()

    k, v = "rt", PARAMS["rt"]
    name = re.sub(r"(I1|I2|R1|R2)", v.read, i1.name)

    if not (done_files[1]).exists():
        logger.info(f"Remapping. Reading from {lane}_lig.tsv.gz.")
        subprocess.run(
            f"pigz -dc {path}/{lane}_lig.tsv.gz | {rustpath/'remap'} {idxs}/shortlong.tsv | pigz > {path}/{lane}_ligmode.tsv.gz",
            shell=True,
            env=os.environ,
            check=True,
        )
        done_files[1].touch()

    if not (done_files[2]).exists():
        shell = os.environ.get("SHELL", "/bin/sh")
        if shell.endswith("fish"):
            substitution = f"(pigz -dc {path}/{lane}_ligmode.tsv.gz | psub)"
        else:
            substitution = f"<(pigz -dc {path}/{lane}_ligmode.tsv.gz)"

        subprocess.run(
            f'{shell} -c "pigz -dc {path/name} | {rustpath/"splitter"} - {idxs.resolve()}/{v.file} --index {v.index} --tol {v.tol}'
            f' -s {substitution} | pigz > {path}/{lane}_rt.tsv.gz"',
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
