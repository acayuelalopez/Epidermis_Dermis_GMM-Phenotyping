#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fromQuPathIntensities_phenotyping_FAST.py

Versión acelerada del phenotyping desde exports single-cell de QuPath.

Principales cambios de velocidad frente al script original:
1) Cada archivo *_single_cell se lee UNA sola vez, y solo con las columnas necesarias.
2) Las columnas de intensidad se renombran a <MARKER>__value para evitar búsquedas repetidas.
3) Las GMM se calculan sobre datos ya cargados y filtrados por Parent, sin releer CSVs por marcador.
4) El cálculo de silhouette usa una muestra configurable pequeña.
5) Los plots pesados se pueden desactivar o submuestrear.

Este script está pensado para tus datos:
- single-cell: columna Image = "2025_11_18__0309.czi - Scene #1"
- metadata.xlsx: columna Image Name = "2025_11_18__0309"
- región: columna Parent = "Annotation (dermis)" / "Annotation (epidermis)"
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score

# =============================================================================
# CONFIGURACIÓN
# =============================================================================

INPUT_DIR = "/mnt/imgserver/IA/Projects/2026/2026_2_10_lpoluha/qupath_0_tune_process/measurements_cellpose_process"
OUTPUT_DIR = "/mnt/imgserver/IA/Projects/2026/2026_2_10_lpoluha/qupath_0_tune_process/output/ref_dist_std"

METADATA_XLSX = "/mnt/imgserver/IA/Projects/2026/2026_2_10_lpoluha/misc/info_1_batch_process.xlsx"
METADATA_SHEET = None
METADATA_IMAGE_COL = "Image"

POPULATIONS_TXT = "/mnt/imgserver/IA/Projects/2026/2026_2_10_lpoluha/misc/populations.txt"

GMM_GROUPING_COLUMN = "Ref"
PLOT_GROUPINGS = [GMM_GROUPING_COLUMN, "TYPE", "DIAGNOSIS", "Ref", "COMPONENT_TYPE", "Image Name", "SIMPLIFIED_NAME"]
# "DIAGNOSIS_2.0"
FILE_PATTERN = "*_single_cell.*"

QUPATH_IMAGE_COL = "Image"
PARENT_COL = "Classification"
X_COL = "Centroid X µm"
Y_COL = "Centroid Y µm"

REGIONS = None  # None = detectar automáticamente

# Regiones jerárquicas. Para la región principal se incluirán también todas
# las células clasificadas dentro de las subregiones indicadas.
# Las subregiones siguen generando además sus propios resultados independientes.
HIERARCHICAL_REGIONS = {
    "Annotation (Lesion)": [
        "Annotation (Lesion_dermis)",
        "Annotation (Lesion_epidermis)",
    ],
}

PHENOTYPE_COL = "phenotype"
POPULATION_COL = "population"

TRANSFORM = "log1p"  # "log1p" o "none"

# Metodo por defecto para calcular el threshold de positividad.
# Valores permitidos: "GMM" o "Single Distribution".
DEFAULT_THRESHOLD_METHOD = "Single Distribution"

# Solo se usa cuando el metodo es "Single Distribution".
# 0 = media; 1 = media + 1*SD; 2 = media + 2*SD.
DEFAULT_SINGLE_DISTRIBUTION_STD_MULTIPLIER = 1.0

# Excepciones por marcador. Si un marcador no aparece aquí, usa los valores
# DEFAULT indicados arriba. Los nombres no distinguen mayúsculas/minúsculas.
MARKER_THRESHOLD_SETTINGS = {
     "S100B": {"method": "GMM"},
     "MDK": {"method": "GMM"},
     "K14": {"method": "GMM"},
    "PTN": {"method": "Single Distribution", "std_multiplier": 1.0},
}

# Compartimento por defecto para el cálculo de los GMM.
# Valores permitidos: "Cell", "Nucleus", "Cytoplasm" o "Membrane".
DEFAULT_INTENSITY_COMPARTMENT = "Cell"

# Excepciones por marcador. Si un marcador no aparece aquí, usa el valor anterior.
# Los nombres de marcador y compartimento no distinguen mayúsculas/minúsculas.
MARKER_INTENSITY_COMPARTMENTS = {
    #"S100B": "Nucleus",
    #"MDK": "Cytoplasm",
    # "K14": "Membrane",
    # "PTN": "Cell",
}

# True = incluir no clasificadas en plots; False = ocultarlas solo en plots.
# Los CSV siempre conservan todas las células.
PLOT_UNCLASSIFIED_CELLS = True
UNCLASSIFIED_PHENOTYPE_LABELS = {
    "Negative", "Unknown", "UNKNOWN", "Unclassified", ""
}
K_MIN = 2
K_MAX = 8
BIC_FRACTION = 0.01
ELBOW_MARGIN = 2
RANDOM_STATE = 10

# -----------------------------
# OPCIONES DE VELOCIDAD
# -----------------------------
# Baja esto si sigue lento. 30000-50000 suele ser mucho más rápido que 200000.
MAX_GMM_CELLS = 500000000

# silhouette es lo más caro; 3000-10000 suele ser suficiente para selección de k.
SILHOUETTE_SAMPLE_SIZE = 10000

# Si quieres ir mucho más rápido, pon False.
ENABLE_GMM_PLOTS = True
ENABLE_PHENOTYPE_PLOTS = True
ENABLE_SPATIAL_PLOTS = True

# Submuestreo para spatial plots por Scene/Image. 0 = usar todas.
SPATIAL_MAX_POINTS_PER_SCENE = 0  # 20000

# Columna QC opcional. Si no quieres QC, dejar "".
QC_COL = ""


# =============================================================================
# UTILIDADES
# =============================================================================

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def detect_sep(path: Path) -> str:
    with open(path, "rb") as fh:
        head = fh.read(4096).decode("utf-8", errors="ignore")
    return "\t" if head.count("\t") > head.count(",") else ","


def safe_filename(s: object) -> str:
    s = str(s).replace("/", "_").replace("\\", "_")
    return re.sub(r"[^A-Za-z0-9._\-]+", "_", s).strip("_") or "NA"


def normalize_for_match(s: object) -> str:
    s = str(s).strip().lower()
    for ext in [".ome.tiff", ".ome.tif", ".czi", ".tiff", ".tif"]:
        s = s.replace(ext, "")
    s = re.sub(r"\s+", "", s)
    return s


def file_candidates(input_dir: Path, pattern: str) -> List[Path]:
    allowed = {".csv", ".tsv", ".txt"}
    return sorted([p for p in input_dir.glob(pattern) if p.suffix.lower() in allowed])


def read_header(path: Path) -> List[str]:
    return pd.read_csv(path, sep=detect_sep(path), nrows=0, engine="python").columns.tolist()


def find_col(df: pd.DataFrame, target: str, required: bool = True) -> Optional[str]:
    for c in df.columns:
        if str(c).strip().lower() == target.lower():
            return c
    for c in df.columns:
        if target.lower() in str(c).strip().lower():
            return c
    if required:
        raise KeyError(f"No se encontró la columna '{target}'")
    return None


def derive_image_name_from_qupath_image(image_value: object) -> str:
    s = str(image_value).strip()
    s = re.sub(r"\.czi.*$", "", s, flags=re.IGNORECASE)
    s = re.sub(r"\s+-\s+Scene\s+#?\d+.*$", "", s, flags=re.IGNORECASE)
    return s.strip()


# =============================================================================
# POPULATIONS / MARKERS
# =============================================================================

# Marcadores permitidos en populations.txt y en las columnas del CSV.
#
# El orden más largo primero evita coincidencias parciales si en el futuro
# introduces marcadores cuyos nombres sean prefijos de otros.
ALLOWED_MARKERS = [
    "S100B",
    "DAPI",
    "K14",
    "MDK",
    "PTN",
]

# DAPI está presente en el CSV, pero normalmente no debe utilizarse
# para definir poblaciones positivas o negativas.
#
# Si quieres utilizar DAPI en populations.txt, elimínalo de este conjunto.
EXCLUDED_PHENOTYPE_MARKERS = {
    "DAPI",
}


def read_populations(pop_txt: Path) -> List[str]:
    """
    Lee las poblaciones, elimina líneas vacías, comentarios y duplicados.
    """
    populations = []

    for line in pop_txt.read_text(
            encoding="utf-8",
            errors="ignore",
    ).splitlines():

        population = line.strip()

        if not population:
            continue

        if population.startswith("#"):
            continue

        populations.append(population)

    output = []
    seen = set()

    for population in populations:

        if population not in seen:
            output.append(population)
            seen.add(population)

    return output


def population_pattern() -> re.Pattern:
    """
    Construye un patrón para encontrar expresiones como:

        S100B+
        K14-
        MDK+
        PTN-

    Los nombres más largos se comprueban primero.
    """
    markers_sorted = sorted(
        ALLOWED_MARKERS,
        key=len,
        reverse=True,
    )

    marker_expression = "|".join(
        re.escape(marker)
        for marker in markers_sorted
    )

    return re.compile(
        rf"({marker_expression})([+-])",
        flags=re.IGNORECASE,
    )


POPULATION_PATTERN = population_pattern()


def parse_population(population: str) -> List[Tuple[str, bool]]:
    """
    Convierte una población en condiciones marcador/estado.

    Ejemplos
    --------
    S100B+
        -> [("S100B", True)]

    S100B-K14-
        -> [("S100B", False), ("K14", False)]

    S100B+MDK+PTN-
        -> [
            ("S100B", True),
            ("MDK", True),
            ("PTN", False),
        ]
    """
    conditions = []

    for match in POPULATION_PATTERN.finditer(str(population)):

        marker = match.group(1).upper()
        sign = match.group(2)

        canonical_marker = next(
            (
                allowed_marker
                for allowed_marker in ALLOWED_MARKERS
                if allowed_marker.upper() == marker
            ),
            marker,
        )

        if canonical_marker in EXCLUDED_PHENOTYPE_MARKERS:
            continue

        conditions.append(
            (
                canonical_marker,
                sign == "+",
            )
        )

    return conditions


def extract_markers(population: str) -> List[str]:
    """
    Devuelve los marcadores encontrados en una población,
    independientemente de que sean positivos o negativos.
    """
    return [
        marker
        for marker, expected_positive
        in parse_population(population)
    ]


def build_marker_list_from_populations(
        populations: List[str],
) -> List[str]:
    """
    Construye la lista única de marcadores utilizados
    en populations.txt.
    """
    markers = []

    for population in populations:
        markers.extend(
            extract_markers(population)
        )

    output = []
    seen = set()

    for marker in markers:

        if marker not in seen:
            output.append(marker)
            seen.add(marker)

    return output


def validate_populations(
        populations: List[str],
) -> None:
    """
    Verifica que todas las poblaciones tengan al menos una
    condición válida.
    """
    invalid_populations = []

    for population in populations:

        conditions = parse_population(population)

        if not conditions:
            invalid_populations.append(population)

    if invalid_populations:
        invalid_text = "\n".join(
            f"  - {population}"
            for population in invalid_populations
        )

        raise ValueError(
            "Hay poblaciones que no contienen marcadores válidos:\n"
            f"{invalid_text}\n\n"
            "Marcadores permitidos: "
            + ", ".join(ALLOWED_MARKERS)
        )


def assign_population_from_bools(
        populations: List[str],
        bools: Dict[str, np.ndarray],
) -> np.ndarray:
    """
    Asigna poblaciones respetando los signos + y -.

    Una condición positiva requiere:
        bools[marker] == True

    Una condición negativa requiere:
        bools[marker] == False

    Las poblaciones con más condiciones tienen prioridad.
    """
    if not bools:
        return np.array([], dtype=str)

    number_of_cells = len(
        next(iter(bools.values()))
    )

    if number_of_cells == 0:
        return np.array([], dtype=str)

    parsed_populations = []

    for original_index, population in enumerate(populations):

        conditions = parse_population(population)

        if conditions:
            parsed_populations.append(
                (
                    population,
                    conditions,
                    original_index,
                )
            )

    # Prioridad:
    # 1. Mayor número de condiciones.
    # 2. Orden original de populations.txt.
    parsed_populations.sort(
        key=lambda item: (
            -len(item[1]),
            item[2],
        )
    )

    output = np.array(
        ["Negative"] * number_of_cells,
        dtype=object,
    )

    for population, conditions, original_index in parsed_populations:

        population_mask = np.ones(
            number_of_cells,
            dtype=bool,
        )

        for marker, expected_positive in conditions:

            marker_values = bools.get(marker)

            if marker_values is None:
                marker_values = np.zeros(
                    number_of_cells,
                    dtype=bool,
                )

            if expected_positive:

                population_mask &= marker_values

            else:

                population_mask &= ~marker_values

        assignable_mask = (
                (output == "Negative") &
                population_mask
        )

        output[assignable_mask] = population

    return output.astype(str)


# =============================================================================
# INTENSIDADES
# =============================================================================
VALID_INTENSITY_COMPARTMENTS = {
    "cell": "Cell",
    "nucleus": "Nucleus",
    "cytoplasm": "Cytoplasm",
    "membrane": "Membrane",
}

def normalize_intensity_compartment(compartment: str) -> str:
    key = str(compartment).strip().lower()
    if key not in VALID_INTENSITY_COMPARTMENTS:
        allowed = ", ".join(VALID_INTENSITY_COMPARTMENTS.values())
        raise ValueError(
            f"Compartimento no válido: {compartment!r}. Permitidos: {allowed}"
        )
    return VALID_INTENSITY_COMPARTMENTS[key]

def compartment_for_marker(marker: str) -> str:
    overrides = {
        str(name).strip().lower(): value
        for name, value in MARKER_INTENSITY_COMPARTMENTS.items()
    }
    selected = overrides.get(
        str(marker).strip().lower(), DEFAULT_INTENSITY_COMPARTMENT
    )
    return normalize_intensity_compartment(selected)

def find_compartment_mean_column(
        columns: Iterable[str], marker: str, compartment: str
) -> Optional[str]:
    """Busca la columna Mean del marcador y compartimento solicitados."""
    columns = list(columns)
    compartment = normalize_intensity_compartment(compartment)
    exact = f"{marker}: {compartment}: Mean"
    if exact in columns:
        return exact

    marker_low = str(marker).strip().lower()
    compartment_low = compartment.lower()
    candidates = []
    for column in columns:
        column_low = str(column).strip().lower()
        if marker_low in column_low and compartment_low in column_low and "mean" in column_low:
            candidates.append(column)
    if not candidates:
        return None

    other_compartments = [
        value.lower() for value in VALID_INTENSITY_COMPARTMENTS.values()
        if value.lower() != compartment_low
    ]
    def score(column: str) -> int:
        column_low = str(column).strip().lower()
        result = 100 if column_low == exact.lower() else 0
        if f": {compartment_low}: mean" in column_low:
            result += 20
        elif f"{compartment_low}: mean" in column_low:
            result += 10
        result -= 20 * sum(x in column_low for x in other_compartments)
        result -= 5 * sum(x in column_low for x in ["std", "median", "min", "max"])
        return result
    candidates.sort(key=lambda column: (-score(column), str(column)))
    return candidates[0]

def transform_values(x: np.ndarray, mode: str) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if mode == "none":
        return x
    if mode == "log1p":
        return np.log1p(np.clip(x, a_min=0, a_max=None))
    raise ValueError(f"Transform no soportada: {mode}")


def inverse_transform_value(x: float, mode: str) -> float:
    return float(np.expm1(x)) if mode == "log1p" else float(x)


# =============================================================================
# METADATA
# =============================================================================

def load_metadata_xlsx(path: Path, sheet: Optional[str], image_col: str) -> pd.DataFrame:
    meta = pd.read_excel(path, sheet_name=sheet if sheet else 0, engine="openpyxl")
    col_img = find_col(meta, image_col, required=True)
    meta = meta.copy()
    if col_img != METADATA_IMAGE_COL:
        meta[METADATA_IMAGE_COL] = meta[col_img]
    meta[METADATA_IMAGE_COL] = meta[METADATA_IMAGE_COL].astype(str).str.strip()
    meta["__image_name_norm__"] = meta[METADATA_IMAGE_COL].map(normalize_for_match)
    return meta


def read_qupath_image_value(path: Path) -> str:
    try:
        df0 = pd.read_csv(path, sep=detect_sep(path), nrows=1, engine="python")
        if QUPATH_IMAGE_COL in df0.columns and len(df0):
            return str(df0.iloc[0][QUPATH_IMAGE_COL])
    except Exception:
        pass
    return ""


def match_metadata_for_file(path: Path, meta: pd.DataFrame) -> pd.Series:
    filename_norm = normalize_for_match(path.name)
    qupath_image = read_qupath_image_value(path)
    qupath_norm = normalize_for_match(qupath_image)
    target = filename_norm + " " + qupath_norm

    matches = []
    for idx, row in meta.iterrows():
        im = row["__image_name_norm__"]
        if im and im in target:
            matches.append((len(im), idx))
    if matches:
        matches.sort(reverse=True)
        return meta.loc[matches[0][1]]

    fallback = pd.Series({c: np.nan for c in meta.columns})
    fallback[METADATA_IMAGE_COL] = derive_image_name_from_qupath_image(qupath_image) if qupath_image else Path(
        path.name).stem.replace("_single_cell", "")
    fallback["__image_name_norm__"] = normalize_for_match(fallback[METADATA_IMAGE_COL])
    return fallback


def build_file_metadata(files: List[Path], meta: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for f in files:
        m = match_metadata_for_file(f, meta)
        row = m.to_dict()
        row["source_file"] = f.name
        row["source_path"] = str(f)
        row["qupath_image_from_file"] = read_qupath_image_value(f)
        rows.append(row)
    return pd.DataFrame(rows)


# =============================================================================
# LECTURA ACELERADA: CADA ARCHIVO UNA SOLA VEZ
# =============================================================================

def load_all_cells_fast(files: List[Path], file_meta: pd.DataFrame, markers: List[str]) -> pd.DataFrame:
    fmap = file_meta.set_index("source_file")
    chunks = []

    for i, f in enumerate(files, start=1):
        print(f"[LOAD {i}/{len(files)}] {f.name}")
        header = read_header(f)
        marker_compartments = {
            marker: compartment_for_marker(marker) for marker in markers
        }
        marker_cols = {
            marker: find_compartment_mean_column(
                header, marker, marker_compartments[marker]
            )
            for marker in markers
        }
        marker_cols = {m: c for m, c in marker_cols.items() if c is not None}
        print(f"\n[INTENSITY COLUMNS] File: {f.name}")

        for marker in markers:
            selected_column = marker_cols.get(marker)

            compartment = marker_compartments[marker]
            if selected_column is None:
                print(f"  WARNING: {marker} -> {compartment}: Mean column not found")
            else:
                print(f"  {marker} [{compartment}] -> {selected_column}")
        base_cols = [QUPATH_IMAGE_COL, PARENT_COL, X_COL, Y_COL]
        if QC_COL:
            base_cols.append(QC_COL)
        usecols = [c for c in base_cols if c in header] + list(marker_cols.values())
        usecols = list(dict.fromkeys(usecols))

        df = pd.read_csv(f, sep=detect_sep(f), usecols=usecols, engine="c")
        df["source_file"] = f.name

        # Renombrar columnas de intensidades a nombres estándar
        rename = {col: f"{m}__value" for m, col in marker_cols.items()}
        df = df.rename(columns=rename)

        # Asegurar columnas de marcadores aunque falten en un archivo
        for m in markers:
            vcol = f"{m}__value"
            if vcol not in df.columns:
                df[vcol] = np.nan
            else:
                df[vcol] = pd.to_numeric(df[vcol], errors="coerce").astype("float32")

        if X_COL in df.columns:
            df[X_COL] = pd.to_numeric(df[X_COL], errors="coerce").astype("float32")
        if Y_COL in df.columns:
            df[Y_COL] = pd.to_numeric(df[Y_COL], errors="coerce").astype("float32")

        # Añadir metadata sin sobrescribir Image de QuPath
        meta_row = fmap.loc[f.name]
        for col, val in meta_row.items():
            if col in ["source_file", "source_path", "__image_name_norm__"]:
                continue
            if col == QUPATH_IMAGE_COL:
                df[f"metadata__{col}"] = val
            else:
                df[col] = val

        if METADATA_IMAGE_COL not in df.columns or df[METADATA_IMAGE_COL].isna().all():
            if QUPATH_IMAGE_COL in df.columns:
                df[METADATA_IMAGE_COL] = df[QUPATH_IMAGE_COL].map(derive_image_name_from_qupath_image)
            else:
                df[METADATA_IMAGE_COL] = Path(f.name).stem.replace("_single_cell", "")

        chunks.append(df)

    if not chunks:
        return pd.DataFrame()
    return pd.concat(chunks, ignore_index=True, copy=False)


VALID_THRESHOLD_METHODS = {
    "gmm": "GMM",
    "single distribution": "Single Distribution",
    "single_distribution": "Single Distribution",
    "single-distribution": "Single Distribution",
}

def normalize_threshold_method(method: str) -> str:
    key = str(method).strip().lower()
    if key not in VALID_THRESHOLD_METHODS:
        allowed = ", ".join(sorted(set(VALID_THRESHOLD_METHODS.values())))
        raise ValueError(
            f"Método de threshold no válido: {method!r}. Permitidos: {allowed}"
        )
    return VALID_THRESHOLD_METHODS[key]

def threshold_settings_for_marker(marker: str) -> Dict[str, object]:
    settings_by_marker = {
        str(name).strip().lower(): settings
        for name, settings in MARKER_THRESHOLD_SETTINGS.items()
    }
    marker_settings = settings_by_marker.get(str(marker).strip().lower(), {})
    if isinstance(marker_settings, str):
        marker_settings = {"method": marker_settings}
    method = normalize_threshold_method(
        marker_settings.get("method", DEFAULT_THRESHOLD_METHOD)
    )
    multiplier = float(
        marker_settings.get(
            "std_multiplier", DEFAULT_SINGLE_DISTRIBUTION_STD_MULTIPLIER
        )
    )
    if not np.isfinite(multiplier) or multiplier < 0:
        raise ValueError(
            f"std_multiplier no válido para {marker}: {multiplier}. "
            "Debe ser un número finito >= 0."
        )
    return {"method": method, "std_multiplier": multiplier}

def single_distribution_threshold_one_marker(
        values_raw: np.ndarray, std_multiplier: float
) -> Dict[str, object]:
    """Threshold = media + std_multiplier * desviación estándar en escala raw."""
    x_raw = np.asarray(values_raw, dtype=float)
    x_raw = x_raw[np.isfinite(x_raw)]
    if x_raw.size < 2:
        return {
            "status": "too_few_values",
            "n": int(x_raw.size),
            "n_raw_before_sample": int(x_raw.size),
            "threshold_raw": np.nan,
            "threshold_method": "Single Distribution",
            "std_multiplier": float(std_multiplier),
        }
    mean_raw = float(np.mean(x_raw))
    std_raw = float(np.std(x_raw, ddof=1))
    threshold_raw = float(mean_raw + float(std_multiplier) * std_raw)
    return {
        "status": "ok",
        "n": int(x_raw.size),
        "n_raw_before_sample": int(x_raw.size),
        "threshold_raw": threshold_raw,
        "threshold_transformed": np.nan,
        "raw_mean": mean_raw,
        "raw_std": std_raw,
        "std_multiplier": float(std_multiplier),
        "threshold_method": "Single Distribution",
        "values_raw": x_raw,
    }

def threshold_one_marker(values_raw: np.ndarray, marker: str) -> Dict[str, object]:
    settings = threshold_settings_for_marker(marker)
    if settings["method"] == "GMM":
        result = gmm_threshold_one_marker(values_raw)
        result["threshold_method"] = "GMM"
        result["std_multiplier"] = np.nan
        return result
    return single_distribution_threshold_one_marker(
        values_raw, settings["std_multiplier"]
    )

# =============================================================================
# GMM
# =============================================================================

def _fit_gmm_1d(X: np.ndarray, n_components: int, random_state: int = 10) -> GaussianMixture:
    gmm = GaussianMixture(n_components=n_components, covariance_type="full", random_state=random_state)
    gmm.fit(X)
    return gmm


def _ssq_to_nearest_mean(X: np.ndarray, means: np.ndarray) -> float:
    dist2 = (X - means.reshape(1, -1)) ** 2
    return float(np.sum(np.min(dist2, axis=1)))


def _safe_silhouette_1d(X: np.ndarray, labels: np.ndarray, sample_size: int, random_state: int) -> float:
    if len(np.unique(labels)) < 2:
        return float("nan")
    try:
        ss = min(sample_size, X.shape[0]) if X.shape[0] > sample_size else None
        return float(silhouette_score(X, labels, sample_size=ss, random_state=random_state))
    except Exception:
        return float("nan")


def evaluate_gmm_k_sweep(values: np.ndarray) -> pd.DataFrame:
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < max(10, K_MIN * 3):
        return pd.DataFrame(columns=["k", "ssq", "silhouette", "bic"])

    X = v.reshape(-1, 1)
    rows = []
    max_allowed = min(K_MAX, max(K_MIN, v.size // 3))
    for k in range(K_MIN, max_allowed + 1):
        try:
            gmm = _fit_gmm_1d(X, k, RANDOM_STATE)
            labels = gmm.predict(X)
            rows.append({
                "k": int(k),
                "ssq": _ssq_to_nearest_mean(X, gmm.means_.reshape(-1)),
                "silhouette": _safe_silhouette_1d(X, labels, SILHOUETTE_SAMPLE_SIZE, RANDOM_STATE),
                "bic": float(gmm.bic(X)),
            })
        except Exception:
            continue
    return pd.DataFrame(rows)


def detect_elbow_k(k_values: np.ndarray, costs: np.ndarray) -> Optional[int]:
    k = np.asarray(k_values, dtype=float)
    y = np.asarray(costs, dtype=float)
    mask = np.isfinite(k) & np.isfinite(y)
    k, y = k[mask], y[mask]
    if len(k) < 3:
        return int(k[0]) if len(k) else None

    p1 = np.array([k[0], y[0]], dtype=float)
    p2 = np.array([k[-1], y[-1]], dtype=float)
    denom = np.linalg.norm(p2 - p1)
    if denom == 0:
        return int(k[0])

    dists = []
    for xi, yi in zip(k, y):
        p = np.array([xi, yi], dtype=float)
        d = abs((p2[0] - p1[0]) * (p1[1] - p[1]) - (p1[0] - p[0]) * (p2[1] - p1[1])) / denom
        dists.append(d)
    return int(k[int(np.argmax(dists))])


def choose_optimal_k(metrics: pd.DataFrame) -> Tuple[Optional[int], Optional[int]]:
    if metrics.empty:
        return None, None
    m = metrics.dropna(subset=["k", "ssq", "bic"]).copy()
    if m.empty:
        return None, None
    elbow_k = detect_elbow_k(m["k"].to_numpy(), m["ssq"].to_numpy())
    if elbow_k is None:
        elbow_k = int(m.loc[m["bic"].idxmin(), "k"])

    allowed = m[(m["k"] >= elbow_k) & (m["k"] <= elbow_k + int(ELBOW_MARGIN))].copy()
    if allowed.empty:
        allowed = m.copy()
    bic_min = allowed["bic"].min()
    bic_delta = abs(bic_min) * float(BIC_FRACTION)
    candidates = allowed[allowed["bic"] <= bic_min + bic_delta].copy()
    if candidates.empty:
        candidates = allowed[allowed["bic"] == bic_min].copy()

    if candidates["silhouette"].notna().any():
        candidates = candidates.sort_values(["silhouette", "k"], ascending=[False, True])
    else:
        candidates = candidates.sort_values(["bic", "k"], ascending=[True, True])
    return int(candidates.iloc[0]["k"]), int(elbow_k)


def subsample_values(values: np.ndarray, max_cells: int, random_state: int) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if max_cells and max_cells > 0 and values.size > max_cells:
        rng = np.random.default_rng(random_state)
        idx = rng.choice(values.size, size=max_cells, replace=False)
        return values[idx]
    return values


def gmm_threshold_one_marker(values_raw: np.ndarray) -> Dict[str, object]:
    x_raw = np.asarray(values_raw, dtype=float)
    x_raw = x_raw[np.isfinite(x_raw)]
    if TRANSFORM == "log1p":
        x_raw = x_raw[x_raw >= 0]

    x_fit_raw = subsample_values(x_raw, MAX_GMM_CELLS, RANDOM_STATE)
    x_t = transform_values(x_fit_raw, TRANSFORM)
    if x_t.size < 20:
        return {"status": "too_few_values", "n": int(x_t.size), "threshold_raw": np.nan}

    metrics = evaluate_gmm_k_sweep(x_t)
    k_best, elbow_k = choose_optimal_k(metrics)
    if k_best is None:
        return {"status": "gmm_failed", "n": int(x_t.size), "threshold_raw": np.nan, "metrics": metrics}

    X = x_t.reshape(-1, 1)
    gmm = _fit_gmm_1d(X, k_best, RANDOM_STATE)
    means = np.sort(gmm.means_.reshape(-1))
    threshold_t = float(np.mean(means))
    threshold_raw = inverse_transform_value(threshold_t, TRANSFORM)

    return {
        "status": "ok",
        "n": int(x_t.size),
        "n_raw_before_sample": int(x_raw.size),
        "k_best": int(k_best),
        "elbow_k": int(elbow_k) if elbow_k is not None else np.nan,
        "threshold_transformed": threshold_t,
        "threshold_raw": threshold_raw,
        "means_transformed": ";".join(f"{z:.6g}" for z in means),
        "weights": ";".join(f"{z:.6g}" for z in gmm.weights_.reshape(-1)),
        "metrics": metrics,
        "gmm": gmm,
        "values_transformed": x_t,
    }


# =============================================================================
# PLOTS
# =============================================================================

def save_single_distribution_plot(
        values_raw: np.ndarray, mean_raw: float, std_raw: float,
        threshold_raw: float, std_multiplier: float, outpath: Path,
        title: str, xlabel: str
) -> None:
    if values_raw is None or len(values_raw) == 0:
        return
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(values_raw, bins=80, alpha=0.75, label="Distribución de intensidad")
    ax.axvline(mean_raw, linestyle="-", linewidth=1.5, label=f"Media = {mean_raw:.4g}")
    ax.axvline(
        threshold_raw, linestyle="--", linewidth=2,
        label=f"Threshold = media + {std_multiplier:g}×SD = {threshold_raw:.4g}"
    )
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Número de células")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def save_gmm_evaluation_plot(metrics: pd.DataFrame, outpath: Path, title: str, chosen_k, elbow_k) -> None:
    if not ENABLE_GMM_PLOTS or metrics is None or metrics.empty:
        return
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, col, ylabel in zip(axes, ["ssq", "silhouette", "bic"], ["SSQ", "Silhouette", "BIC"]):
        ax.plot(metrics["k"], metrics[col], marker="o")
        if elbow_k is not None and np.isfinite(elbow_k):
            ax.axvline(elbow_k, linestyle="--", linewidth=1, label="elbow_k")
        if chosen_k is not None and np.isfinite(chosen_k):
            ax.axvline(chosen_k, linestyle=":", linewidth=1.5, label="chosen_k")
        ax.set_xlabel("k")
        ax.set_ylabel(ylabel)
        ax.legend(fontsize=7)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def save_hist_gmm_plot(values_t: np.ndarray, gmm: GaussianMixture, threshold_t: float, outpath: Path, title: str,
                       xlabel: str) -> None:
    if not ENABLE_GMM_PLOTS or values_t is None or gmm is None or len(values_t) == 0:
        return
    x = np.asarray(values_t, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(x, bins=80, density=True, alpha=0.45)
    grid = np.linspace(np.nanmin(x), np.nanmax(x), 500).reshape(-1, 1)
    total_pdf = np.exp(gmm.score_samples(grid))
    ax.plot(grid[:, 0], total_pdf, linewidth=2, label="GMM total")
    means = gmm.means_.reshape(-1)
    vars_ = gmm.covariances_.reshape(gmm.n_components, -1)[:, 0]
    for w, mu, var in zip(gmm.weights_, means, vars_):
        sd = math.sqrt(max(float(var), 1e-12))
        comp_pdf = float(w) * (1.0 / (sd * math.sqrt(2 * math.pi))) * np.exp(-0.5 * ((grid[:, 0] - mu) / sd) ** 2)
        ax.plot(grid[:, 0], comp_pdf, linewidth=1)
    ax.axvline(threshold_t, linestyle="--", linewidth=1.5, label="threshold")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Density")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def phenotype_colors(phenotypes: Iterable[str]) -> Dict[str, object]:
    cmap = plt.get_cmap("tab20")
    return {p: cmap(i % 20) for i, p in enumerate(list(phenotypes))}


def plot_stacked_per_image(proportions: pd.DataFrame, out_png: Path, color_map: Dict[str, object]) -> None:
    if proportions.empty:
        return
    fig, ax = plt.subplots(figsize=(max(10, len(proportions) * 0.45), 7))
    proportions.plot(kind="bar", stacked=True, ax=ax, width=0.8, color=[color_map.get(c) for c in proportions.columns])
    ax.set_xlabel(METADATA_IMAGE_COL)
    ax.set_ylabel("Population / detections in annotation")
    ax.set_title("Phenotype abundance per image, normalized by annotation detections")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=7)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def plot_stacked_by_group(proportions: pd.DataFrame, image_meta: pd.DataFrame, grouping: str, out_dir: Path,
                          color_map: Dict[str, object]) -> None:
    if grouping not in image_meta.columns or proportions.empty:
        return
    tmp = proportions.copy()
    tmp[grouping] = tmp.index.map(image_meta[grouping])
    tmp = tmp.dropna(subset=[grouping])
    if tmp.empty:
        return
    means = tmp.groupby(grouping).mean(numeric_only=True)
    fig, ax = plt.subplots(figsize=(max(8, len(means) * 1.5), 7))
    means.plot(kind="bar", stacked=True, ax=ax, color=[color_map.get(c) for c in means.columns])
    ax.set_ylabel("Mean population / detections in annotation")
    ax.set_title(f"Phenotype abundance by {grouping}")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right")
    fig.tight_layout()
    fig.savefig(out_dir / f"stacked_by_{safe_filename(grouping)}.png", dpi=150)
    plt.close(fig)


def plot_barplot_by_group(proportions: pd.DataFrame, image_meta: pd.DataFrame, grouping: str, out_dir: Path) -> None:
    if grouping not in image_meta.columns or proportions.empty:
        return
    tmp = proportions.copy()
    tmp[grouping] = tmp.index.map(image_meta[grouping])
    tmp = tmp.dropna(subset=[grouping])
    if tmp.empty:
        return
    long = tmp.reset_index().rename(columns={tmp.index.name or "index": METADATA_IMAGE_COL})
    long = long.melt(id_vars=[METADATA_IMAGE_COL, grouping], var_name="phenotype",
                     value_name="normalized_by_annotation_detections")
    groups = sorted(long[grouping].dropna().unique().tolist(), key=lambda x: str(x))
    phenos = sorted(long["phenotype"].dropna().unique().tolist(), key=lambda x: str(x))
    if not groups or not phenos:
        return
    x = np.arange(len(phenos))
    width = 0.8 / max(1, len(groups))
    fig, ax = plt.subplots(figsize=(max(10, len(phenos) * 0.65), 6))
    rng = np.random.default_rng(10)
    for gi, g in enumerate(groups):
        data_g = long[long[grouping] == g]
        means = [data_g.loc[data_g["phenotype"] == p, "normalized_by_annotation_detections"].mean() for p in phenos]
        xpos = x - 0.4 + width / 2 + gi * width
        ax.bar(xpos, means, width=width, label=str(g), alpha=0.75)
        for pi, p in enumerate(phenos):
            vals = data_g.loc[data_g["phenotype"] == p, "normalized_by_annotation_detections"].dropna().to_numpy()
            if len(vals):
                jitter = rng.normal(0, width * 0.08, size=len(vals))
                ax.scatter(np.full(len(vals), xpos[pi]) + jitter, vals, s=14, alpha=0.8, color="black")
    ax.set_xticks(x)
    ax.set_xticklabels(phenos, rotation=45, ha="right")
    ax.set_ylabel("Population / detections in annotation")
    ax.set_title(f"Phenotype abundance by {grouping}")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / f"barplot_by_{safe_filename(grouping)}.png", dpi=150)
    plt.close(fig)


def plot_spatial_by_qupath_image(df: pd.DataFrame, out_dir: Path) -> None:
    if not ENABLE_SPATIAL_PLOTS:
        return
    if X_COL not in df.columns or Y_COL not in df.columns or QUPATH_IMAGE_COL not in df.columns:
        return
    ensure_dir(out_dir)
    colors = phenotype_colors(sorted(df[PHENOTYPE_COL].dropna().unique().tolist()))
    rng = np.random.default_rng(RANDOM_STATE)
    for img, sub in df.groupby(QUPATH_IMAGE_COL, dropna=False, sort=False):
        if SPATIAL_MAX_POINTS_PER_SCENE and len(sub) > SPATIAL_MAX_POINTS_PER_SCENE:
            sub = sub.iloc[rng.choice(len(sub), size=SPATIAL_MAX_POINTS_PER_SCENE, replace=False)]
        fig, ax = plt.subplots(figsize=(7, 7))
        for p, ss in sub.groupby(PHENOTYPE_COL, dropna=False, sort=False):
            ax.scatter(pd.to_numeric(ss[X_COL], errors="coerce"), pd.to_numeric(ss[Y_COL], errors="coerce"),
                       s=3, alpha=0.65, label=str(p), color=colors.get(p))
        ax.invert_yaxis()
        ax.set_xlabel(X_COL)
        ax.set_ylabel(Y_COL)
        ax.set_title(str(img))
        ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=6, markerscale=3)
        fig.tight_layout()
        fig.savefig(out_dir / f"spatial_{safe_filename(img)}.png", dpi=150)
        plt.close(fig)


def phenotype_reports_and_plots(region_df: pd.DataFrame, out_region: Path, plot_groupings: List[str]) -> None:
    """
    Reports y plots normalizados por número de detecciones de la anotación.

    Se llama dentro de una región concreta de Parent, por ejemplo:
    'Annotation (dermis)' o 'Annotation (epidermis)'.

    Por tanto, para cada imagen:
        denominador = total de células/detecciones de esa imagen dentro de esa anotación.

    Para cada población:
        valor_normalizado = n_células_de_la_población / denominador
    """
    reports = ensure_dir(out_region / "reports")
    plots = ensure_dir(out_region / "plots")
    if region_df.empty:
        return

    # -------------------------------------------------------------------------
    # NORMALIZACIÓN POR IMAGEN BASE, usando Image Name.
    # Denominador = total de detecciones/células en esta anotación por imagen.
    # -------------------------------------------------------------------------
    counts_image = pd.crosstab(region_df[METADATA_IMAGE_COL], region_df[PHENOTYPE_COL])
    annotation_detections_image = region_df.groupby(METADATA_IMAGE_COL).size().rename("annotation_detection_count")
    normalized_image = counts_image.div(annotation_detections_image, axis=0).fillna(0)

    report_image = counts_image.add_prefix("count__")
    report_image = report_image.join(annotation_detections_image)
    report_image = report_image.join(normalized_image.add_prefix("normalized_by_annotation_detections__"))

    report_image.to_csv(reports / "phenotype_counts_normalized_by_annotation_detections_per_image.csv")
    # Mantengo el nombre previo para compatibilidad, pero el contenido queda explícitamente normalizado.
    report_image.to_csv(reports / "phenotype_counts_proportions_per_image.csv")

    long_counts_image = counts_image.reset_index().melt(
        id_vars=METADATA_IMAGE_COL,
        var_name=PHENOTYPE_COL,
        value_name="n_cells",
    )
    long_norm_image = normalized_image.reset_index().melt(
        id_vars=METADATA_IMAGE_COL,
        var_name=PHENOTYPE_COL,
        value_name="normalized_by_annotation_detections",
    )
    long_image = long_counts_image.merge(long_norm_image, on=[METADATA_IMAGE_COL, PHENOTYPE_COL], how="left")
    long_image = long_image.merge(annotation_detections_image.reset_index(), on=METADATA_IMAGE_COL, how="left")
    long_image.to_csv(reports / "counts_by_image_phenotype_normalized.csv", index=False)

    # -------------------------------------------------------------------------
    # NORMALIZACIÓN POR SCENE / QUPATH IMAGE para trazabilidad.
    # Denominador = total de detecciones/células en esta anotación por Scene.
    # -------------------------------------------------------------------------
    if QUPATH_IMAGE_COL in region_df.columns:
        counts_scene = pd.crosstab(region_df[QUPATH_IMAGE_COL], region_df[PHENOTYPE_COL])
        annotation_detections_scene = region_df.groupby(QUPATH_IMAGE_COL).size().rename("annotation_detection_count")
        normalized_scene = counts_scene.div(annotation_detections_scene, axis=0).fillna(0)

        report_scene = counts_scene.add_prefix("count__")
        report_scene = report_scene.join(annotation_detections_scene)
        report_scene = report_scene.join(normalized_scene.add_prefix("normalized_by_annotation_detections__"))
        report_scene.to_csv(reports / "phenotype_counts_normalized_by_annotation_detections_per_scene.csv")
        report_scene.to_csv(reports / "phenotype_counts_proportions_per_scene.csv")

        long_counts_scene = counts_scene.reset_index().melt(
            id_vars=QUPATH_IMAGE_COL,
            var_name=PHENOTYPE_COL,
            value_name="n_cells",
        )
        long_norm_scene = normalized_scene.reset_index().melt(
            id_vars=QUPATH_IMAGE_COL,
            var_name=PHENOTYPE_COL,
            value_name="normalized_by_annotation_detections",
        )
        long_scene = long_counts_scene.merge(long_norm_scene, on=[QUPATH_IMAGE_COL, PHENOTYPE_COL], how="left")
        long_scene = long_scene.merge(annotation_detections_scene.reset_index(), on=QUPATH_IMAGE_COL, how="left")
        long_scene.to_csv(reports / "counts_by_scene_phenotype_normalized.csv", index=False)

    if not ENABLE_PHENOTYPE_PLOTS:
        return

    # Este filtro afecta únicamente a los plots. Los reports anteriores
    # mantienen todas las células y sus denominadores originales.
    plot_df = region_df
    if not PLOT_UNCLASSIFIED_CELLS:
        unclassified = {
            str(label).strip().lower()
            for label in UNCLASSIFIED_PHENOTYPE_LABELS
        }
        phenotype_normalized = (
            region_df[PHENOTYPE_COL].fillna("").astype(str).str.strip().str.lower()
        )
        plot_df = region_df.loc[~phenotype_normalized.isin(unclassified)].copy()
        print(
            f"Plots: excluidas {len(region_df) - len(plot_df):,} "
            "células no clasificadas"
        )
    if plot_df.empty:
        print("WARNING: no quedan células clasificadas para generar los plots")
        return

    # Al ocultar las no clasificadas, los plots se renormalizan usando solo
    # las células visibles, de modo que las barras apiladas siguen sumando 1.
    plot_counts_image = pd.crosstab(
        plot_df[METADATA_IMAGE_COL], plot_df[PHENOTYPE_COL]
    )
    plot_denominator_image = plot_df.groupby(METADATA_IMAGE_COL).size()
    plot_normalized_image = (
        plot_counts_image.div(plot_denominator_image, axis=0).fillna(0)
    )
    color_map = phenotype_colors(list(plot_normalized_image.columns))
    plot_stacked_per_image(plot_normalized_image,
                           plots / "stacked_barplot_per_image_normalized_by_annotation_detections.png", color_map)
    # Mantengo también el nombre anterior para compatibilidad.
    plot_stacked_per_image(plot_normalized_image, plots / "stacked_barplot_per_image.png", color_map)
    plot_spatial_by_qupath_image(plot_df, plots / "spatial_by_scene")

    image_meta = plot_df.drop_duplicates(METADATA_IMAGE_COL).set_index(METADATA_IMAGE_COL)
    for grouping in plot_groupings:
        if grouping in image_meta.columns:
            plot_stacked_by_group(plot_normalized_image, image_meta, grouping, plots, color_map)
            plot_barplot_by_group(plot_normalized_image, image_meta, grouping, plots)


# =============================================================================
# CORE
# =============================================================================

def normalize_region_name(value: object) -> str:
    """Normaliza espacios externos del nombre de una región."""
    return str(value).strip()


def hierarchical_region_members(region: str) -> List[str]:
    """Devuelve la región solicitada y sus subregiones configuradas."""
    region = normalize_region_name(region)
    children = HIERARCHICAL_REGIONS.get(region, [])
    members = [region] + [normalize_region_name(x) for x in children]
    return list(dict.fromkeys(members))


def discover_regions_from_loaded(df_all: pd.DataFrame) -> List[str]:
    """
    Detecta regiones presentes en Classification.

    Si aparece una subregión configurada, también añade la región principal
    para generar el resultado combinado, aunque no haya filas etiquetadas
    directamente con el nombre de la región principal.
    """
    observed = (
        df_all[PARENT_COL]
        .dropna()
        .astype(str)
        .map(normalize_region_name)
        .drop_duplicates()
        .tolist()
    )
    observed_set = set(observed)
    regions = list(observed)

    for parent, children in HIERARCHICAL_REGIONS.items():
        parent = normalize_region_name(parent)
        child_set = {normalize_region_name(x) for x in children}
        if parent in observed_set or observed_set.intersection(child_set):
            if parent not in regions:
                regions.append(parent)

    return regions


def select_cells_for_analysis_region(
        df_all: pd.DataFrame,
        region: str,
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Selecciona las células de una región.

    Para una región principal incluye sus propias filas y las filas de todas
    sus subregiones. La columna Classification original se conserva para
    trazabilidad.
    """
    members = hierarchical_region_members(region)
    classification = (
        df_all[PARENT_COL]
        .fillna("")
        .astype(str)
        .map(normalize_region_name)
    )
    region_df = df_all.loc[classification.isin(members)].copy()
    region_df["analysis_region"] = normalize_region_name(region)
    region_df["analysis_region_members"] = " | ".join(members)
    region_df["source_annotation_classification"] = region_df[PARENT_COL]
    return region_df, members


def compute_thresholds_for_region(region_df: pd.DataFrame, region: str, markers: List[str], grouping_col: str,
                                  gmm_plot_dir: Path) -> pd.DataFrame:
    rows = []
    for marker in markers:
        vcol = f"{marker}__value"
        if vcol not in region_df.columns:
            continue
        for group_val, sub in region_df.groupby(grouping_col, dropna=False, sort=False):
            group_val = "UNKNOWN" if pd.isna(group_val) else str(group_val)
            vals = pd.to_numeric(sub[vcol], errors="coerce").to_numpy()
            vals = vals[np.isfinite(vals)]
            res = threshold_one_marker(vals, marker)
            rows.append(make_threshold_row(region, grouping_col, group_val, marker, vals, res))
            save_threshold_outputs(gmm_plot_dir, region, group_val, marker, res)

        vals_all = pd.to_numeric(region_df[vcol], errors="coerce").to_numpy()
        vals_all = vals_all[np.isfinite(vals_all)]
        res_all = threshold_one_marker(vals_all, marker)
        rows.append(make_threshold_row(region, grouping_col, "__ALL__", marker, vals_all, res_all))
        save_threshold_outputs(gmm_plot_dir, region, "__ALL__", marker, res_all)

    return pd.DataFrame(rows)


def make_threshold_row(region, grouping_col, group_val, marker, vals, res):
    method = res.get("threshold_method", threshold_settings_for_marker(marker)["method"])
    if method == "GMM":
        strategy = "mean_of_gmm_component_means_elbow_floor"
    else:
        multiplier = res.get("std_multiplier", np.nan)
        strategy = f"raw_mean_plus_{multiplier:g}_times_sample_std"
    return {
        "analysis_region": region,
        "threshold_grouping_column": grouping_col,
        "threshold_group": group_val,
        # Nombres anteriores conservados para compatibilidad con phenotype_region.
        "gmm_grouping_column": grouping_col,
        "gmm_group": group_val,
        "marker": marker,
        "intensity_compartment": compartment_for_marker(marker),
        "threshold_method": method,
        "status": res.get("status"),
        "n_fit": res.get("n", 0),
        "n_raw_before_sample": res.get("n_raw_before_sample", len(vals)),
        "transform": TRANSFORM if method == "GMM" else "none",
        "raw_mean": res.get("raw_mean", np.nan),
        "raw_std": res.get("raw_std", np.nan),
        "std_multiplier": res.get("std_multiplier", np.nan),
        "k_best": res.get("k_best", np.nan),
        "elbow_k": res.get("elbow_k", np.nan),
        "threshold_transformed": res.get("threshold_transformed", np.nan),
        "threshold_raw": res.get("threshold_raw", np.nan),
        "means_transformed": res.get("means_transformed", ""),
        "weights": res.get("weights", ""),
        "threshold_strategy": strategy,
    }


def save_threshold_outputs(gmm_plot_dir: Path, region: str, group_val: str, marker: str, res: Dict[str, object]) -> None:
    if res.get("status") != "ok":
        return
    prefix = f"{safe_filename(region)}__{safe_filename(group_val)}__{marker}"
    method = res.get("threshold_method", "GMM")
    if method == "Single Distribution":
        save_single_distribution_plot(
            res.get("values_raw"), res.get("raw_mean"), res.get("raw_std"),
            res.get("threshold_raw"), res.get("std_multiplier"),
            gmm_plot_dir / f"{prefix}__hist_single_distribution.png",
            f"{region} | {group_val} | {marker}",
            f"{marker} {compartment_for_marker(marker)} Mean (raw)",
        )
        return
    metrics = res.get("metrics")
    if isinstance(metrics, pd.DataFrame):
        metrics.to_csv(gmm_plot_dir / f"{prefix}__gmm_metrics.csv", index=False)
        save_gmm_evaluation_plot(metrics, gmm_plot_dir / f"{prefix}__gmm_evaluation.png",
                                 f"{region} | {group_val} | {marker}", res.get("k_best"), res.get("elbow_k"))
    save_hist_gmm_plot(res.get("values_transformed"), res.get("gmm"), res.get("threshold_transformed"),
                       gmm_plot_dir / f"{prefix}__hist_gmm.png", f"{region} | {group_val} | {marker}",
                       f"{marker} {compartment_for_marker(marker)} Mean ({TRANSFORM})")


def threshold_lookup(thresholds: pd.DataFrame) -> Dict[Tuple[str, str], float]:
    out = {}
    ok = thresholds[thresholds["status"] == "ok"].copy()
    for _, r in ok.iterrows():
        out[(str(r["gmm_group"]), str(r["marker"]))] = float(r["threshold_raw"])
    return out


def phenotype_region(region_df: pd.DataFrame, thresholds: pd.DataFrame, markers: List[str], populations: List[str],
                     grouping_col: str) -> pd.DataFrame:
    df = region_df.copy()
    tlook = threshold_lookup(thresholds)
    bools: Dict[str, np.ndarray] = {}

    group_values = df[grouping_col].fillna("UNKNOWN").astype(
        str).to_numpy() if grouping_col in df.columns else np.array(["UNKNOWN"] * len(df))

    for marker in markers:
        vcol = f"{marker}__value"
        if vcol not in df.columns:
            vals = np.full(len(df), np.nan)
        else:
            vals = pd.to_numeric(df[vcol], errors="coerce").to_numpy()

        thr_arr = np.array([tlook.get((g, marker), tlook.get(("__ALL__", marker), np.nan)) for g in group_values],
                           dtype=float)
        pos = np.isfinite(vals) & np.isfinite(thr_arr) & (vals > thr_arr)
        bools[marker] = pos
        df[f"{marker}__positive"] = pos
        df[f"{marker}__threshold_raw"] = thr_arr

    df[PHENOTYPE_COL] = assign_population_from_bools(populations, bools)
    df[POPULATION_COL] = df[PHENOTYPE_COL]
    df["gmm_grouping_column"] = grouping_col
    df["gmm_group"] = group_values

    # Denominador de las cuantificaciones de poblaciones.
    # Se agrupa por analysis_region, no por Classification original, porque
    # una región principal puede combinar células de varias subregiones.
    if METADATA_IMAGE_COL in df.columns and "analysis_region" in df.columns:
        df["annotation_detection_count_per_image"] = df.groupby(
            [METADATA_IMAGE_COL, "analysis_region"]
        )["analysis_region"].transform("size")
        df["population_count_per_image_annotation"] = df.groupby(
            [METADATA_IMAGE_COL, "analysis_region", PHENOTYPE_COL]
        )[PHENOTYPE_COL].transform("size")
        df["population_fraction_per_image_annotation"] = (
                df["population_count_per_image_annotation"] /
                df["annotation_detection_count_per_image"].replace(0, np.nan)
        )

    # Denominador por Scene/Image de QuPath para trazabilidad espacial.
    if QUPATH_IMAGE_COL in df.columns and "analysis_region" in df.columns:
        df["annotation_detection_count_per_scene"] = df.groupby(
            [QUPATH_IMAGE_COL, "analysis_region"]
        )["analysis_region"].transform("size")
        df["population_count_per_scene_annotation"] = df.groupby(
            [QUPATH_IMAGE_COL, "analysis_region", PHENOTYPE_COL]
        )[PHENOTYPE_COL].transform("size")
        df["population_fraction_per_scene_annotation"] = (
                df["population_count_per_scene_annotation"] /
                df["annotation_detection_count_per_scene"].replace(0, np.nan)
        )

    return df


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    input_dir = Path(INPUT_DIR)
    output_dir = ensure_dir(Path(OUTPUT_DIR))
    metadata_xlsx = Path(METADATA_XLSX)
    populations_txt = Path(POPULATIONS_TXT)

    if not input_dir.exists():
        raise FileNotFoundError(f"No existe INPUT_DIR: {input_dir}")
    if not metadata_xlsx.exists():
        raise FileNotFoundError(f"No existe METADATA_XLSX: {metadata_xlsx}")
    if not populations_txt.exists():
        raise FileNotFoundError(f"No existe POPULATIONS_TXT: {populations_txt}")

    meta = load_metadata_xlsx(metadata_xlsx, METADATA_SHEET, METADATA_IMAGE_COL)
    populations = read_populations(populations_txt)

    if not populations:
        raise RuntimeError(
            "populations.txt está vacío o no contiene "
            "poblaciones válidas"
        )

    validate_populations(populations)

    markers = build_marker_list_from_populations(
        populations
    )

    if not markers:
        raise RuntimeError(
            "No se detectaron marcadores válidos en populations.txt. "
            "Marcadores permitidos: "
            + ", ".join(ALLOWED_MARKERS)
        )

    print(f"Poblaciones detectadas: {len(populations)}")
    print(f"Marcadores de fenotipado: {markers}")

    files = file_candidates(input_dir, FILE_PATTERN)
    if not files:
        raise RuntimeError(f"No se encontraron archivos con patrón {FILE_PATTERN} en {input_dir}")

    print(f"Archivos detectados: {len(files)}")
    print(f"Marcadores: {markers}")
    print(f"Método de threshold por defecto: {normalize_threshold_method(DEFAULT_THRESHOLD_METHOD)}")
    print("Método de threshold usado por marcador:")
    for marker in markers:
        settings = threshold_settings_for_marker(marker)
        if settings["method"] == "Single Distribution":
            print(
                f"  - {marker}: Single Distribution | "
                f"media + {settings['std_multiplier']:g}*SD"
            )
        else:
            print(f"  - {marker}: GMM")
    print(f"Compartimento por defecto: {normalize_intensity_compartment(DEFAULT_INTENSITY_COMPARTMENT)}")
    print("Compartimentos usados por marcador:")
    for marker in markers:
        print(f"  - {marker}: {compartment_for_marker(marker)}")
    print(f"Incluir no clasificadas en plots: {PLOT_UNCLASSIFIED_CELLS}")
    print("Construyendo metadata por archivo...")
    file_meta = build_file_metadata(files, meta)

    grouping_col = GMM_GROUPING_COLUMN
    if grouping_col not in file_meta.columns:
        found = find_col(file_meta, grouping_col, required=False)
        if found is None:
            raise KeyError(f"La columna GMM_GROUPING_COLUMN='{grouping_col}' no existe en metadata.xlsx")
        grouping_col = found

    file_meta.to_csv(output_dir / "file_metadata_matching.csv", index=False)

    print("Cargando single-cell una sola vez...")
    df_all = load_all_cells_fast(files, file_meta, markers)
    if df_all.empty:
        raise RuntimeError("No se cargaron células")

    if QC_COL and QC_COL in df_all.columns:
        # Filtro QC por archivo y región: > media local
        qc = pd.to_numeric(df_all[QC_COL], errors="coerce")
        means = df_all.assign(__qc=qc).groupby(["source_file", PARENT_COL])["__qc"].transform("mean")
        df_all = df_all[qc > means].copy()

    regions = REGIONS if REGIONS is not None else discover_regions_from_loaded(df_all)
    print(f"Regiones: {regions}")

    all_thresholds = []
    output_paths = []

    for region in regions:
        print(f"\n=== Región: {region} ===")
        out_region = ensure_dir(output_dir / "phenotyping" / safe_filename(region))
        reports = ensure_dir(out_region / "reports")
        gmm_plot_dir = ensure_dir(out_region / "gmm_plots")

        region_df, region_members = select_cells_for_analysis_region(
            df_all,
            region,
        )
        print(f"Classification incluidas: {region_members}")
        print(f"Células en región combinada: {len(region_df):,}")
        if len(region_members) > 1 and not region_df.empty:
            member_counts = (
                region_df[PARENT_COL]
                .fillna("<NA>")
                .astype(str)
                .value_counts(dropna=False)
            )
            print("Desglose por Classification original:")
            for member_name, member_count in member_counts.items():
                print(f"  - {member_name}: {member_count:,}")
        if region_df.empty:
            continue

        print("Calculando thresholds por marcador...")
        thresholds = compute_thresholds_for_region(region_df, region, markers, grouping_col, gmm_plot_dir)
        thresholds.to_csv(reports / "thresholds_by_marker.csv", index=False)
        all_thresholds.append(thresholds)

        print("Asignando phenotypes...")
        phenotyped = phenotype_region(region_df, thresholds, markers, populations, grouping_col)
        single_cell_path = out_region / f"single_cell_phenotyped__{safe_filename(region)}.csv"
        phenotyped.to_csv(single_cell_path, index=False)
        output_paths.append(single_cell_path)
        print(f"Guardado: {single_cell_path}")

        print("Reports y plots...")
        phenotype_reports_and_plots(phenotyped, out_region, PLOT_GROUPINGS)

    if all_thresholds:
        pd.concat(all_thresholds, ignore_index=True).to_csv(output_dir / "thresholds_by_marker_all_regions.csv", index=False)
    pd.DataFrame({"single_cell_phenotyped_csv": [str(p) for p in output_paths]}).to_csv(
        output_dir / "single_cell_phenotyped_outputs_index.csv", index=False)

    print("\nProcesamiento completado")
    print(f"Resultados: {output_dir}")


if __name__ == "__main__":
    main()