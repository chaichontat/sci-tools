# %%
import random
from Levenshtein import hamming

from nupack import Model, SetSpec, Strand, Tube, tube_analysis

from scitools.bases import gen_sample_space
from scitools.constants import NEXTERA2, P5, P7, SMALLRNA, TRUSEQ1
from scitools.sci import gen_set, hairpin
from scitools.utils import repeat

my_model = Model(material="dna", celsius=63, sodium=0.3, magnesium=0.003)


def gen_p7(seq: str, type_: str) -> str:
    if type_ == "nextera":
        return P7 + seq.lower() + NEXTERA2
    if type_ == "smallrna":
        return P7 + seq.lower() + SMALLRNA[:19]
    raise ValueError("Invalid primer type")


def gen_p5(seq: str):
    return P5 + seq + TRUSEQ1[:21]


p7_screen = dict(
    nextera=lambda seq: ((mfe := hairpin(gen_p7(seq, "nextera"), my_model)), mfe >= -0.06),
    smallrna=lambda seq: ((mfe := hairpin(gen_p7(seq, "smallrna"), my_model)), mfe >= -0.06),
)


p5_screen = dict(p5=lambda seq: ((mfe := hairpin(gen_p5(seq), my_model)), mfe >= -0.37))

base = gen_sample_space(repeat("N", 10))

random.seed(42)
random.shuffle(base)
# %%
random.seed(10)
p7out = gen_set(base, "".join(random.choices("ATGC", k=10)), screens=p7_screen)
p7out = sorted(p7out, key=lambda x: -x[1]["nextera"])
# %%
base = gen_sample_space(repeat("N", 10))
random.seed(50)
random.shuffle(base)
random.seed(50)
p5out = gen_set(base, "".join(random.choices("ATGC", k=10)), screens=p5_screen, gc_range=(0.4, 0.6))


# %%
chosen = []
it = iter(p5out)
while len(chosen) < 96:
    cand = next(it)
    for s in chosen:
        if hamming(cand[0], s[0]) < 4:
            break
    else:
        chosen.append(cand)
# %%
