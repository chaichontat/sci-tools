import pandas as pd


def gen_plate(name: str, seqs: list[str]) -> pd.DataFrame:
    wells = [f"{row}{col:02d}" for row in "ABCDEFGH" for col in range(1, 13)]
    return pd.DataFrame(
        {
            "Well Position": wells,
            "Name": [name + w for w in wells],
            "Sequence": seqs,
        }
    )


def gen_tube(name: str, seq: str, conc: str = "25nm", purify: str = "STD") -> str:
    return f"{name}\t{conc}\t{purify}\t{seq}"
