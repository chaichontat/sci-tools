from scitools.bases import reverse_complement
from scitools.constants import ADT_LINKER, NEXTERA2, OVERHANG, P5, P7, SMALLRNA, TRUSEQ1


def gen_p7(idx: str) -> str:
    return P7 + idx.lower() + NEXTERA2


def gen_p7adt(idx: str) -> str:
    return ADT_LINKER + idx.lower() + SMALLRNA[:19]


def gen_p5(idx: str):
    return P5 + idx.lower() + TRUSEQ1[:21]


def gen_lig(idx: str):
    # T is presumably a buffer to prevent digestion
    return (
        reverse_complement(OVERHANG)
        + idx.lower()
        + "T"
        + "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
        + reverse_complement(idx).lower()
    )


def gen_rt(idx: str):
    return "/5Phos/" + OVERHANG + idx.lower() + "N" * 9 + "T" * 28 + "VN"
