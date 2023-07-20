# %%
from pathlib import Path

import polars as pl
import xlsxwriter

dfs = {f.stem: pl.read_csv(f) for f in Path("fixmess").glob("sci*.csv")}

with xlsxwriter.Workbook("fixmess/fix.xlsx") as wb:
    for name, df in dfs.items():
        print(name)
        df.select(**{"Well Position": "well", "Name": "name", "Sequence": "seq"}).write_excel(
            wb, worksheet=name
        )

# %%
