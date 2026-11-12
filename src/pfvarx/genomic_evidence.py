"""Shared interval and alignment conventions."""

from collections import defaultdict
from pathlib import Path

FLAG_FILTER = 0x4 | 0x100 | 0x200 | 0x400 | 0x800


def read_bed(path: Path | None) -> dict[str, list[tuple[int, int]]]:
    out: dict[str, list[tuple[int, int]]] = defaultdict(list)
    if path is None or not path.is_file():
        return out
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith(("#", "track")):
            continue
        c = line.split("\t")
        out[c[0]].append((int(c[1]), int(c[2])))
    for k in out:
        out[k].sort()
    return out


def in_intervals(iv: list[tuple[int, int]], pos: int) -> bool:
    lo, hi = 0, len(iv)
    while lo < hi:
        mid = (lo + hi) // 2
        if iv[mid][1] <= pos:
            lo = mid + 1
        else:
            hi = mid
    return lo < len(iv) and iv[lo][0] <= pos < iv[lo][1]
