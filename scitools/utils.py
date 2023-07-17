from binascii import hexlify

import colorama
from hexhamming import hamming_distance_string


def printc(seq: str):
    for c in seq:
        if c == "A" or c == "a":
            print(colorama.Fore.GREEN + c, end="")
        elif c == "T" or c == "t":
            print(colorama.Fore.RED + c, end="")
        elif c == "C" or c == "c":
            print(colorama.Fore.BLUE + c, end="")
        elif c == "G" or c == "g":
            print(colorama.Fore.YELLOW + c, end="")
        else:
            print(colorama.Fore.WHITE + c, end="")
    print(colorama.Fore.RESET)


def hamming(a: bytes, b: bytes) -> int:
    return hamming_distance_string(hexlify(a), hexlify(b)) >> 1


def repeat(seq: str, n: int) -> str:
    return "".join([seq] * n)
