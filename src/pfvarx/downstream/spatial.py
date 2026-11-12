"""DA007 spatial and temporal summary of a source outcome (default: DA002 markers).

Eligibility is inherited from the source analysis: only samples analysed there, with
valid coordinates inside the declared country (Natural Earth 1:110m polygons, bundled)
and a valid collection date, are included. Outputs are aggregated to country x month
and to a 0.5-degree grid (no raw coordinates in shared outputs). A GAM (mgcv, REML)
of marker prevalence on time is fitted only when at least 10 callable samples exist;
otherwise the analysis is descriptive. Preferential sampling is a stated limitation.
"""

from __future__ import annotations

import datetime as dt
import json
import shutil
import subprocess
from pathlib import Path
from typing import Any

from pfvarx.constants import Reason, Terminal
from pfvarx.downstream.common import Context, write_manifest
from pfvarx.io import read_tsv, write_json, write_tsv

GEO = Path(__file__).resolve().parents[3] / "config" / "geo" / "ne_110m_countries.geojson"


def _point_in_ring(x: float, y: float, ring: list[list[float]]) -> bool:
    inside = False
    j = len(ring) - 1
    for i in range(len(ring)):
        xi, yi = ring[i][0], ring[i][1]
        xj, yj = ring[j][0], ring[j][1]
        if (yi > y) != (yj > y) and x < (xj - xi) * (y - yi) / ((yj - yi) or 1e-12) + xi:
            inside = not inside
        j = i
    return inside


def in_country(lon: float, lat: float, country: str, geo: dict[str, Any]) -> bool | None:
    for f in geo["features"]:
        if f["properties"]["name"].lower() != country.lower():
            continue
        g = f["geometry"]
        polys = g["coordinates"] if g["type"] == "MultiPolygon" else [g["coordinates"]]
        return any(_point_in_ring(lon, lat, poly[0]) and not any(_point_in_ring(lon, lat, h) for h in poly[1:])
                   for poly in polys)
    return None  # country not in layer: cannot validate


def run(ctx: Context, threads: int = 1) -> None:
    out = ctx.out("DA007")
    source = ctx.eff["config"]["analyses"]["spatial_outcome_source"]
    src_manifest = ctx.results / "downstream" / source / "status_manifest.json"
    if not src_manifest.is_file() or json.loads(src_manifest.read_text())["terminal_state"] not in (
            "COMPLETED", "COMPLETED_WITH_EXCLUSIONS"):
        write_manifest(ctx, "DA007", Terminal.BLOCKED_MISSING_PREREQUISITE.value,
                       reason_codes=[Reason.SOURCE_OUTCOME_INCOMPLETE.value], requested=ctx.samples, eligible=[],
                       excluded={}, analysed=[], outputs=[], caveats=[], method={"source": source})
        return
    if source != "DA002":
        write_manifest(ctx, "DA007", Terminal.SKIPPED_INELIGIBLE.value, reason_codes=[Reason.METHOD_NOT_VALIDATED.value],
                       requested=ctx.samples, eligible=[], excluded={}, analysed=[], outputs=[],
                       caveats=[f"source {source} has no registered spatial adapter"], method={"source": source})
        return
    geo = json.loads(GEO.read_text()) if GEO.is_file() else {"features": []}
    calls = read_tsv(ctx.results / "downstream" / "DA002" / "marker_calls.tsv")
    analysed = sorted({c["sample_id"] for c in calls})
    valid: dict[str, dict[str, Any]] = {}
    excluded: dict[str, list[str]] = {}
    for sid in analysed:
        m = ctx.metadata.get(sid, {})
        try:
            lat, lon = float(m.get("latitude", "")), float(m.get("longitude", ""))
            ok = in_country(lon, lat, m.get("country", ""), geo)
        except ValueError:
            ok = False
        if ok is False:
            excluded.setdefault(Reason.METADATA_GEOGRAPHY_INVALID.value, []).append(sid)
            continue
        try:
            date = dt.date.fromisoformat(m.get("collection_date", ""))
        except ValueError:
            excluded.setdefault(Reason.METADATA_DATE_INVALID.value, []).append(sid)
            continue
        valid[sid] = {"country": m.get("country"), "month": date.strftime("%Y-%m"),
                      "cell": f"{round(lat * 2) / 2:.1f},{round(lon * 2) / 2:.1f}",
                      "geo_validated": ok is True}
    min_n = int(ctx.params["da007_min_samples_per_unit"])
    agg: dict[tuple[str, str, str], list[str]] = {}
    for c in calls:
        if c["sample_id"] not in valid:
            continue
        v = valid[c["sample_id"]]
        agg.setdefault((c["marker_id"], v["country"], v["month"]), []).append(c["state"])
    rows = []
    for (marker, country, month), states in sorted(agg.items()):
        callable_ = [s for s in states if s != "UNKNOWN"]
        k = sum(1 for s in callable_ if s in ("MUTANT", "MIXED"))
        rows.append({"marker_id": marker, "country": country, "month": month, "samples": len(states),
                     "callable": len(callable_), "mutant_or_mixed": k,
                     "prevalence": f"{k / len(callable_):.4f}" if len(callable_) >= min_n else "",
                     "status": "SUMMARISED" if len(callable_) >= min_n else "INSUFFICIENT_SAMPLES"})
    write_tsv(out / "country_month_summary.tsv", rows,
              ["marker_id", "country", "month", "samples", "callable", "mutant_or_mixed", "prevalence", "status"])
    cells: dict[tuple[str, str], list[str]] = {}
    for c in calls:
        if c["sample_id"] in valid:
            cells.setdefault((c["marker_id"], valid[c["sample_id"]]["cell"]), []).append(c["state"])
    grid = [{"marker_id": m, "cell": cell, "lat": float(cell.split(",")[0]), "lon": float(cell.split(",")[1]),
             "callable": sum(1 for s in st if s != "UNKNOWN"),
             "mutant_or_mixed": sum(1 for s in st if s in ("MUTANT", "MIXED"))} for (m, cell), st in sorted(cells.items())]
    write_json(out / "grid_summary.json", grid)
    model_status = "NOT_FITTED_DESCRIPTIVE_ONLY"
    outputs = [out / "country_month_summary.tsv", out / "grid_summary.json"]
    if len(valid) >= 10 and shutil.which("Rscript"):
        fit_output = out / "gam_fit.tsv"
        fit_output.unlink(missing_ok=True)
        flat = out / "gam_input.tsv"
        recs = [{"sample_id": c["sample_id"], "marker_id": c["marker_id"],
                 "y": 1 if c["state"] in ("MUTANT", "MIXED") else 0,
                 "t": (dt.date.fromisoformat(valid[c["sample_id"]]["month"] + "-01") - dt.date(2000, 1, 1)).days}
                for c in calls if c["sample_id"] in valid and c["state"] != "UNKNOWN"]
        write_tsv(flat, recs, ["sample_id", "marker_id", "y", "t"])
        r = subprocess.run(["Rscript", str(ctx.root / "workflow" / "scripts" / "da007_gam.R"), str(flat),
                            str(fit_output)], capture_output=True, text=True)
        (out / "gam.log").write_text(r.stdout + r.stderr)
        if r.returncode != 0 or not fit_output.is_file():
            model_status = "FIT_FAILED"
        elif not read_tsv(fit_output):
            model_status = "INSUFFICIENT_MODEL_SUPPORT"
        else:
            model_status = "FITTED"
            outputs.append(fit_output)
    write_manifest(ctx, "DA007", (Terminal.COMPLETED_WITH_EXCLUSIONS if excluded else Terminal.COMPLETED).value,
                   reason_codes=[], requested=ctx.samples, eligible=sorted(valid), excluded=excluded,
                   analysed=sorted(valid), outputs=outputs,
                   caveats=["Descriptive aggregation; preferential sampling is not modelled.",
                            "DA007 summaries use a 0.5-degree grid. The general dashboard can contain sample coordinates.",
                            "Eligibility inherited from the source analysis " + source + "."],
                   method={"source": source, "aggregation": "country x month; 0.5-degree grid", "model": model_status,
                           "geo_layer": "Natural Earth 1:110m countries (public domain)"},
                   results_summary={"units": len(rows)})
