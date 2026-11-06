#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Genera Object Classifiers JSON de QuPath desde thresholds GMM."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import time
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

# =============================================================================
# CONFIGURACION EDITABLE
# =============================================================================
THRESHOLDS_CSV = Path(
    "/mnt/imgserver/IA/Projects/2026/2026_2_10_lpoluha/qupath_3 _tune (copy)/output/ref_dist_std/thresholds_by_marker_all_regions.csv"
)
OUTPUT_DIR = Path(
    "/mnt/imgserver/IA/Projects/2026/2026_2_10_lpoluha/qupath_3 _tune (copy)/output/ref_dist_std/qupath_object_classifiers"
)

# None = exportar todos. Ejemplos: ["__ALL__"] o ["90170"]
GROUPS_TO_EXPORT = None
REGIONS_TO_EXPORT = None

# IMPORTANTE: QuPath debe encontrar estos nombres exactos en sus mediciones.
# Correspondencia marcador biologico -> canal de la imagen
MEASUREMENT_BY_MARKER = {
    "S100B": "AF751: Cell: Mean",  # S100 beta
    "K14": "AF488: Cell: Mean",
    "MDK": "AF555: Cell: Mean",
    "PTN": "AF647: Cell: Mean",
    "DAPI": "AF405: Cell: Mean",
}

POPULATIONS = [
    #"S100B+",
    # "K14+",
    #"S100B-K14-",
    #"S100B+MDK+PTN-",
    "PTN+",
    "S100B+PTN+MDK-",
    "S100B+MDK+PTN+",
    #"K14+MDK+PTN-",
    "K14+PTN+MDK-",
    "K14+MDK+PTN+",
    #"S100B-K14-MDK+PTN-",
    "S100B-K14-PTN+MDK-",
    "S100B-K14-MDK+PTN+",
    #"S100B+MDK+",
    "S100B+PTN+",
    "S100B+MDK+PTN+",
]

FILTER = "DETECTIONS_ALL"
INDENT = 2
TOKEN_PATTERN = re.compile(r"(S100B|K14|MDK|PTN|DAPI)([+-])", re.IGNORECASE)


def unique(items: Iterable[str]) -> List[str]:
    return list(dict.fromkeys(items))


def safe_filename(value: object) -> str:
    text = re.sub(r"[^A-Za-z0-9._+-]+", "_", str(value)).strip("_")
    return text or "NA"


def parse_population(population: str) -> List[Tuple[str, bool]]:
    conditions: List[Tuple[str, bool]] = []
    rebuilt = ""
    canonical = {k.upper(): k for k in MEASUREMENT_BY_MARKER}
    for match in TOKEN_PATTERN.finditer(population):
        marker_raw, sign = match.groups()
        marker = canonical[marker_raw.upper()]
        conditions.append((marker, sign == "+"))
        rebuilt += marker + sign
    if not conditions or rebuilt.upper() != population.upper():
        raise ValueError(f"Poblacion no valida: {population}")
    return conditions


def stable_color(name: str) -> List[int]:
    digest = hashlib.sha256(name.encode("utf-8")).digest()
    return [55 + digest[i] % 181 for i in range(3)]


def selected(value: str, selection: Sequence[str] | None) -> bool:
    return selection is None or value in selection


def load_thresholds(path: Path) -> Dict[Tuple[str, str, str], float]:
    if not path.exists():
        raise FileNotFoundError(f"No existe THRESHOLDS_CSV: {path}")

    lookup: Dict[Tuple[str, str, str], float] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"analysis_region", "gmm_group", "marker", "status", "threshold_raw"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise KeyError(f"Faltan columnas: {sorted(missing)}")

        for row in reader:
            if row["status"].strip().lower() != "ok":
                continue
            region = row["analysis_region"].strip()
            group = row["gmm_group"].strip()
            marker = row["marker"].strip()
            if not selected(region, REGIONS_TO_EXPORT) or not selected(group, GROUPS_TO_EXPORT):
                continue
            if marker not in MEASUREMENT_BY_MARKER:
                continue
            try:
                threshold = float(row["threshold_raw"])
            except (TypeError, ValueError):
                continue
            if not math.isfinite(threshold):
                continue
            key = (region, group, marker)
            if key in lookup and not math.isclose(lookup[key], threshold, abs_tol=1e-12):
                raise ValueError(f"Threshold duplicado y diferente para {key}")
            lookup[key] = threshold
    return lookup


def simple_classifier(marker: str, positive: bool, threshold: float, timestamp: int) -> dict:
    class_name = f"{marker}{'+' if positive else '-'}"
    path_class = {"name": class_name, "color": stable_color(class_name)}
    function = {
        "classifier_fun": "ClassifyByMeasurementFunction",
        "measurement": MEASUREMENT_BY_MARKER[marker],
        "threshold": float(threshold),
    }
    if positive:
        # Misma logica que Python: intensidad > threshold
        function["pathClassAbove"] = path_class
    else:
        # Misma logica que Python: intensidad <= threshold
        function["pathClassBelow"] = path_class
        function["pathClassEquals"] = path_class

    return {
        "object_classifier_type": "SimpleClassifier",
        "function": function,
        "pathClasses": [path_class],
        "filter": FILTER,
        "timestamp": timestamp,
    }


def main() -> None:
    lookup = load_thresholds(THRESHOLDS_CSV)
    populations = unique(POPULATIONS)
    parsed = {p: parse_population(p) for p in populations}
    region_groups = sorted({(region, group) for region, group, _ in lookup})
    if not region_groups:
        raise RuntimeError("No hay thresholds exportables.")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = []
    base_timestamp = int(time.time() * 1000)
    generated = 0

    for region, group in region_groups:
        available = {
            marker: threshold
            for (r, g, marker), threshold in lookup.items()
            if r == region and g == group
        }
        target_dir = OUTPUT_DIR / safe_filename(region) / safe_filename(group)
        target_dir.mkdir(parents=True, exist_ok=True)

        for population, conditions in parsed.items():
            markers = unique(marker for marker, _ in conditions)
            missing = [marker for marker in markers if marker not in available]
            if missing:
                manifest.append({
                    "analysis_region": region,
                    "gmm_group": group,
                    "population": population,
                    "status": "skipped_missing_thresholds",
                    "missing_markers": ";".join(missing),
                    "json_path": "",
                })
                continue

            classifiers = [
                simple_classifier(marker, positive, available[marker], base_timestamp + generated * 20 + i)
                for i, (marker, positive) in enumerate(conditions)
            ]
            output = {
                "object_classifier_type": "CompositeClassifier",
                "classifiers": classifiers,
            }
            output_path = target_dir / f"{safe_filename(population)}.json"
            output_path.write_text(
                json.dumps(output, indent=INDENT, ensure_ascii=False) + "\n",
                encoding="utf-8",
            )
            generated += 1
            manifest.append({
                "analysis_region": region,
                "gmm_group": group,
                "population": population,
                "status": "generated",
                "missing_markers": "",
                "json_path": str(output_path),
            })

    manifest_path = OUTPUT_DIR / "manifest_qupath_object_classifiers.csv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        fields = ["analysis_region", "gmm_group", "population", "status", "missing_markers", "json_path"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(manifest)

    print(f"Clasificadores generados: {generated}")
    print(f"Salida: {OUTPUT_DIR}")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
