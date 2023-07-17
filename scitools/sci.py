from typing import Callable, Iterable, TypeVar

from Levenshtein import distance
from nupack import Model, SetSpec, Strand, Tube, tube_analysis

from scitools.bases import gc_content

try:
    profile
except NameError:
    profile = lambda x: x


def hairpin(seq: str, model: Model, conc: float = 1e-7) -> float:
    s = Strand(seq.upper(), "s")
    t1 = Tube(strands={s: conc}, name="t1", complexes=SetSpec(max_size=2))
    tube_result = tube_analysis(tubes=[t1], compute=["mfe"], model=model)
    return tube_result["(s)"].mfe_stack


T = TypeVar("T")


@profile
def gen_set(
    samples: Iterable[str],
    starter: str,
    *,
    screens: dict[str, Callable[[str], tuple[T, bool]]],
    max_gen: int = 800,
    min_dist: int = 4,
    gc_range: tuple[float, float] = (0.3, 0.6),
) -> list[tuple[str, dict[str, T]]]:
    out = [
        (
            starter,
            {name: f(starter)[0] for name, f in screens.items()},
        )
    ]
    for i, x in enumerate(samples):
        if len(out) > max_gen:
            break

        if i % 5000 == 0:
            print(i, len(out))

        seq = "".join(x)
        if seq.endswith("CC"):
            continue
        if "GGGG" in seq or "CCCC" in seq or "AAAA" in seq or "TTTT" in seq:
            continue
        if gc_content(seq) < gc_range[0] or gc_content(seq) > gc_range[1]:
            continue
        for x in out:
            if distance(seq, x[0]) < min_dist:
                break
        else:
            res: dict[str, T] = {}
            for name, f in screens.items():
                if not (r := f(seq))[1]:
                    break
                res[name] = r[0]
            else:
                out.append((seq, res))
    return out
