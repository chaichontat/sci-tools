#%%
import random
from binascii import hexlify
from collections import Counter
from itertools import product
from typing import Callable

import colorama
import pandas as pd
import primer3
import seaborn as sns
from hexhamming import hamming_distance_string
from Levenshtein import distance
from nupack import Model, SetSpec, Strand, Tube, tube_analysis

#%%


def gen_rt():
    bridges = pd.read_csv("./ps_bridges.tsv", sep="\t")
    htvn = "".join(["H"] * 12) + "".join(["T"] * 23) + "VN"
    rtplate = gen_plate(
        "sciv2-RT-", ["/5Phos/CTCACTG" + bridges.iloc[s[1], 1][:10].lower() + htvn for s in selected]
    )
