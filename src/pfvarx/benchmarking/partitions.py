"""Immutable benchmarking / held-out partition of truth samples, locked before tuning."""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import random
from pathlib import Path
from typing import Any

from pfvarx.io import read_json, write_json


def lock(root: Path, holdout_fraction: float, seed: int) -> dict[str, Any]:
    truth = read_json(root / "inputs" / "truth" / "truth.json")
    path = root / "inputs" / "partitions" / "partition.json"
    if path.is_file():
        raise ValueError(f"partition already locked: {path} (partitions are immutable)")
    ids = sorted(truth["samples"])
    rng = random.Random(seed)
    k = max(1, round(holdout_fraction * len(ids)))
    held = sorted(rng.sample(ids, k))
    part = {"schema": "pfvarx.partition/1", "locked_utc": dt.datetime.now(dt.UTC).isoformat(),
            "seed": seed, "benchmarking": [i for i in ids if i not in held], "held_out": held}
    part["sha256"] = hashlib.sha256(json.dumps(part, sort_keys=True).encode()).hexdigest()
    write_json(path, part)
    return part


def load(root: Path) -> dict[str, Any]:
    path = root / "inputs" / "partitions" / "partition.json"
    if not path.is_file():
        raise ValueError("no locked partition: run 'pfvarx benchmark-partition' before evaluation")
    part = read_json(path)
    body = {k: v for k, v in part.items() if k != "sha256"}
    if hashlib.sha256(json.dumps(body, sort_keys=True).encode()).hexdigest() != part["sha256"]:
        raise ValueError("partition file was modified after locking (hash mismatch)")
    return part
