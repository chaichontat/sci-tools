import subprocess
from pathlib import Path

import click


@click.command()
@click.argument(
    "path",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path),
)
@click.option("--lanes", type=str)
@click.option(
    "--idxs",
    type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=Path),
    default=Path.cwd(),
)
@click.option("--output", type=str, required=True)
@click.option("--threads", "-t", type=int, default=16)
def run(path: Path, lanes: str, output: str, idxs: Path, threads: int = 16):
    kfiles = {
        k: idxs / v for k, v in dict(i="index.idx", g="t2g.txt", c1="cdna.txt", c2="nascent.txt").items()
    }
    assert all(v.exists() for v in kfiles.values())
    indices = " ".join([f"-{k} {v}" for k, v in kfiles.items()])

    laness = list(map(int, lanes.split(",")))
    files = [f"{path}/L00{i}_combi.fastq.gz {path}/Undetermined_S0_L00{i}_R2_001.fastq.gz" for i in laness]

    subprocess.run(
        f"kb count --h5ad -x '0,0,28:0,28,37:1,0,0' -t {threads} --workflow nac"
        f" {indices} --strand forward"
        f' {" ".join(files)} -o {path}/{output} --overwrite --verbose'
        f" -w NONE",
        shell=True,
        check=True,
    )


if __name__ == "__main__":
    run()
