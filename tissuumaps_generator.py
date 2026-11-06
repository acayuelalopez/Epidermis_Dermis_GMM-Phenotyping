#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
prepare_tissuumaps_from_single_cell.py

Genera archivos CSV de markers/células para abrir en TissUUmaps a partir de un CSV
single-cell phenotyped que contiene múltiples imágenes/escenas.

Entrada esperada en el CSV:
    - Image
    - Centroid X µm
    - Centroid Y µm
    - phenotype

Cada fila es una célula/detección.
El script genera un CSV de markers por cada imagen/escena, con columnas sencillas:
    - x
    - y
    - phenotype
    - cell_id
    - image
    - parent/opcionales si existen

IMPORTANTE SOBRE COORDENADAS:
TissUUmaps normalmente espera coordenadas en el sistema de coordenadas de la imagen.
Si tu imagen se abre en píxeles y tus centroides vienen en micras, necesitas convertir:
    x_pixels = Centroid X µm / MICRONS_PER_PIXEL_X
    y_pixels = Centroid Y µm / MICRONS_PER_PIXEL_Y
Si quieres conservar exactamente las coordenadas en µm, deja MICRONS_PER_PIXEL_X/Y = 1.0.
"""

from __future__ import annotations

from pathlib import Path
import re
import json
import math
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# =============================================================================
# CONFIGURACIÓN: CAMBIA ESTAS VARIABLES Y EJECUTA EN PYCHARM
# =============================================================================

# Carpeta donde están las imágenes que quieres abrir en TissUUmaps.
# Pueden ser .tif/.tiff/.ome.tif/.ome.tiff/.svs/.ndpi/.png/.jpg, etc.
INPUT_IMAGES_DIR = "/mnt/imgserver/IA/Projects/2026/2026_2_10_lpoluha/images/all"

# Archivo single-cell phenotyped con todas las células de todas las imágenes.
INPUT_SINGLE_CELL_CSV = "/mnt/imgserver/IA/Projects/2026/2026_5_26_ezamorano/output/single_cell_phenotyped_outputs_index.csv"

# Carpeta de salida donde se guardarán los CSV por imagen/escena.
OUTPUT_DIR = "/mnt/imgserver/IA/Projects/2026/2026_5_26_ezamorano/output"

# Columnas del CSV single-cell.
IMAGE_COL = "Image"
X_UM_COL = "Centroid X µm"
Y_UM_COL = "Centroid Y µm"
PHENOTYPE_COL = "phenotype"
PARENT_COL = ""  # opcional, si existe se conserva

# Si quieres agrupar exactamente por cada valor de la columna Image, deja esto en True.
# Es lo recomendado si Image contiene escenas tipo:
#   2025_11_18__0309.czi - Scene #1
#   2025_11_18__0309.czi - Scene #2
GROUP_BY_EXACT_IMAGE_VALUE = False

# Si quieres generar también archivos para imágenes presentes en INPUT_IMAGES_DIR aunque no haya células,
# pon True. Normalmente False.
CREATE_EMPTY_FILES_FOR_IMAGES_WITHOUT_CELLS = False

# Conversión de micras a píxeles.
# Si TissUUmaps necesita coordenadas en píxeles, pon aquí el pixel size de QuPath.
# Ejemplo: si 1 pixel = 0.325 µm, usa 0.325.
# Si quieres dejar x/y en micras tal cual, deja 1.0.
MICRONS_PER_PIXEL_X = 0.3450
MICRONS_PER_PIXEL_Y = 0.3450

# Nombres de columnas que se exportarán como coordenadas para TissUUmaps.
# TissUUmaps puede mapear columnas x/y fácilmente.
OUTPUT_X_COL = "x"
OUTPUT_Y_COL = "y"

# Columnas extra que quieres mantener si existen.
EXTRA_COLUMNS_TO_KEEP = [
    "Image Name",
    "source_file",
    "Parent",
    "analysis_region",
    "population",
    "Ref",
    "DIAGNOSIS",
    "TYPE",
    "DIAGNOSIS_2.0",
    "COMPONENT_TYPE",
    "annotation_detection_count_per_image",
    "population_count_per_image_annotation",
    "population_fraction_per_image_annotation",
]

# Extensiones de imagen que se buscarán en INPUT_IMAGES_DIR.
IMAGE_EXTENSIONS = [
    ".tif", ".tiff", ".ome.tif", ".ome.tiff",
    ".svs", ".ndpi", ".vsi", ".scn", ".qptiff",
    ".png", ".jpg", ".jpeg",
    ".dzi",
    ".czi",  # ojo: TissUUmaps puede requerir conversión si no abre CZI directamente
]

# Separador de salida. Para TissUUmaps, coma suele ir bien.
OUTPUT_SEPARATOR = ","

# Si True, genera un JSON adicional por imagen con metadatos útiles.
WRITE_METADATA_JSON = True


# =============================================================================
# FUNCIONES AUXILIARES
# =============================================================================

def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def detect_sep(path: Path) -> str:
    """Detecta separador tab/coma leyendo los primeros KB."""
    with open(path, "rb") as fh:
        head = fh.read(8192).decode("utf-8", errors="ignore")
    return "\t" if head.count("\t") > head.count(",") else ","


def safe_filename(value: object, max_len: int = 180) -> str:
    """Convierte un nombre de imagen/escena en un nombre seguro para archivo."""
    s = str(value).strip()
    s = s.replace(".czi - Scene #", "_Scene_")
    s = s.replace(" - Scene #", "_Scene_")
    s = s.replace("/", "_").replace("\\", "_")
    s = re.sub(r"[^A-Za-z0-9._\-]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    if len(s) > max_len:
        s = s[:max_len].rstrip("_")
    return s or "image"


def normalize_for_matching(value: object) -> str:
    """Normaliza strings para emparejar nombres de imagen con valores de la columna Image."""
    s = str(value).strip().lower()
    for ext in [".ome.tiff", ".ome.tif", ".tiff", ".tif", ".svs", ".ndpi", ".vsi", ".scn", ".qptiff", ".png", ".jpg", ".jpeg", ".dzi", ".czi"]:
        s = s.replace(ext, "")
    s = s.replace(" - scene #", "_scene_")
    s = s.replace("scene #", "scene_")
    s = re.sub(r"\s+", "", s)
    s = re.sub(r"[^a-z0-9]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s


def list_image_files(input_dir: Path) -> List[Path]:
    files: List[Path] = []
    for p in input_dir.rglob("*"):
        if not p.is_file():
            continue
        name_lower = p.name.lower()
        if any(name_lower.endswith(ext) for ext in IMAGE_EXTENSIONS):
            files.append(p)
    return sorted(files)


def read_single_cell_csv(path: Path) -> pd.DataFrame:
    sep = detect_sep(path)
    df = pd.read_csv(path, sep=sep, engine="python")
    return df


def validate_columns(df: pd.DataFrame) -> None:
    required = [IMAGE_COL, X_UM_COL, Y_UM_COL, PHENOTYPE_COL]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise KeyError(
            "Faltan columnas obligatorias en el CSV single-cell: " + ", ".join(missing) +
            "\nColumnas disponibles:\n" + "\n".join(map(str, df.columns))
        )


def find_matching_image_file(image_value: str, image_files: List[Path]) -> Optional[Path]:
    """
    Intenta encontrar el archivo de imagen correspondiente a un valor de la columna Image.
    Ejemplo de image_value:
        2025_11_18__0309.czi - Scene #1
    """
    if not image_files:
        return None

    target_norm = normalize_for_matching(image_value)

    # 1) match por nombre completo/stem normalizado
    candidates: List[Tuple[int, Path]] = []
    for img_path in image_files:
        name_norm = normalize_for_matching(img_path.name)
        stem_norm = normalize_for_matching(img_path.stem)

        score = 0
        if name_norm == target_norm or stem_norm == target_norm:
            score = 100
        elif name_norm and name_norm in target_norm:
            score = 80 + len(name_norm)
        elif stem_norm and stem_norm in target_norm:
            score = 70 + len(stem_norm)
        elif target_norm and target_norm in name_norm:
            score = 60 + len(target_norm)
        elif target_norm and target_norm in stem_norm:
            score = 50 + len(target_norm)

        if score > 0:
            candidates.append((score, img_path))

    if not candidates:
        return None

    candidates.sort(key=lambda x: x[0], reverse=True)
    return candidates[0][1]


def build_tissuumaps_table(sub: pd.DataFrame) -> pd.DataFrame:
    """Crea la tabla de markers/células para TissUUmaps."""
    out = pd.DataFrame()

    x_um = pd.to_numeric(sub[X_UM_COL], errors="coerce")
    y_um = pd.to_numeric(sub[Y_UM_COL], errors="coerce")

    out[OUTPUT_X_COL] = x_um / float(MICRONS_PER_PIXEL_X)
    out[OUTPUT_Y_COL] = y_um / float(MICRONS_PER_PIXEL_Y)

    # Mantener también coordenadas originales en micras para trazabilidad.
    out["x_um"] = x_um
    out["y_um"] = y_um

    out["phenotype"] = sub[PHENOTYPE_COL].astype(str).fillna("Unknown")
    out["image"] = sub[IMAGE_COL].astype(str)

    # ID simple por célula dentro del archivo exportado.
    out["cell_id"] = np.arange(1, len(out) + 1)

    # Columnas extra si existen.
    for col in EXTRA_COLUMNS_TO_KEEP:
        if col in sub.columns and col not in out.columns:
            out[col] = sub[col].values

    # Quitar filas sin coordenadas válidas.
    out = out[np.isfinite(out[OUTPUT_X_COL]) & np.isfinite(out[OUTPUT_Y_COL])].copy()

    return out


def make_color_table(phenotypes: List[str]) -> pd.DataFrame:
    """
    Tabla auxiliar de colores por phenotype.
    TissUUmaps permite colorear por columna categórica; esta tabla sirve como referencia.
    """
    # Paleta estable y simple. Si hay más fenotipos, se recicla.
    palette = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
        "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173",
        "#3182bd", "#31a354", "#756bb1", "#e6550d", "#969696",
    ]
    rows = []
    for i, ph in enumerate(sorted(map(str, phenotypes))):
        rows.append({"phenotype": ph, "color_hex": palette[i % len(palette)]})
    return pd.DataFrame(rows)


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    input_images_dir = Path(INPUT_IMAGES_DIR)
    input_csv = Path(INPUT_SINGLE_CELL_CSV)
    output_dir = ensure_dir(Path(OUTPUT_DIR))
    markers_dir = ensure_dir(output_dir / "markers_csv")
    metadata_dir = ensure_dir(output_dir / "metadata")

    if not input_csv.exists():
        raise FileNotFoundError(f"No existe INPUT_SINGLE_CELL_CSV: {input_csv}")

    if not input_images_dir.exists():
        print(f"AVISO: INPUT_IMAGES_DIR no existe: {input_images_dir}")
        print("Se generarán CSV por valor de Image, pero image_file quedará vacío en el index.")
        image_files = []
    else:
        image_files = list_image_files(input_images_dir)

    print("Leyendo CSV single-cell...")
    df = read_single_cell_csv(input_csv)
    validate_columns(df)

    print(f"Filas/células totales: {len(df):,}")
    print(f"Imágenes/escenas únicas en columna {IMAGE_COL!r}: {df[IMAGE_COL].nunique():,}")
    print(f"Archivos de imagen encontrados en INPUT_IMAGES_DIR: {len(image_files):,}")

    # Limpiar filas sin Image o phenotype.
    df = df[df[IMAGE_COL].notna()].copy()
    df[PHENOTYPE_COL] = df[PHENOTYPE_COL].fillna("Unknown").astype(str)

    # Si se agrupa por Image exacto, cada valor de Image genera un markers CSV.
    groups = list(df.groupby(IMAGE_COL, sort=True, dropna=False))

    index_rows = []
    all_phenotypes = set()

    for idx, (image_value, sub) in enumerate(groups, start=1):
        image_value = str(image_value)
        print(f"[{idx}/{len(groups)}] Preparando markers para: {image_value}")

        markers = build_tissuumaps_table(sub)
        if markers.empty:
            print("  - Sin coordenadas válidas, se omite.")
            continue

        all_phenotypes.update(markers["phenotype"].dropna().astype(str).unique().tolist())

        image_file = find_matching_image_file(image_value, image_files)
        safe_name = safe_filename(image_value)
        markers_file = markers_dir / f"{safe_name}__tissuumaps_markers.csv"
        markers.to_csv(markers_file, index=False, sep=OUTPUT_SEPARATOR)

        phenotype_counts = markers["phenotype"].value_counts(dropna=False).to_dict()

        index_row = {
            "image_value_in_csv": image_value,
            "matched_image_file": str(image_file) if image_file else "",
            "markers_csv": str(markers_file),
            "n_cells": int(len(markers)),
            "n_phenotypes": int(markers["phenotype"].nunique()),
            "phenotypes": ";".join(sorted(map(str, markers["phenotype"].unique()))),
            "x_column_for_tissuumaps": OUTPUT_X_COL,
            "y_column_for_tissuumaps": OUTPUT_Y_COL,
            "phenotype_column_for_tissuumaps": "phenotype",
            "microns_per_pixel_x_used": MICRONS_PER_PIXEL_X,
            "microns_per_pixel_y_used": MICRONS_PER_PIXEL_Y,
        }
        index_rows.append(index_row)

        if WRITE_METADATA_JSON:
            meta = dict(index_row)
            meta["phenotype_counts"] = {str(k): int(v) for k, v in phenotype_counts.items()}
            json_file = metadata_dir / f"{safe_name}__metadata.json"
            json_file.write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")

    # Tabla de colores por phenotype.
    color_table = make_color_table(sorted(all_phenotypes))
    color_table_path = output_dir / "phenotype_color_reference.tsv"
    color_table.to_csv(color_table_path, index=False, sep="\t")

    # Índice general.
    index_df = pd.DataFrame(index_rows)
    index_path = output_dir / "tissuumaps_files_index.csv"
    index_df.to_csv(index_path, index=False)

    # Instrucciones breves.
    instructions = f"""
Archivos generados para TissUUmaps
===================================

Carpeta markers:
  {markers_dir}

Índice:
  {index_path}

Tabla de colores de referencia:
  {color_table_path}

Cómo abrir en TissUUmaps, por cada imagen:
1. Abre TissUUmaps.
2. Carga la imagen correspondiente de INPUT_IMAGES_DIR.
3. Carga el CSV de markers correspondiente desde la columna 'markers_csv' del índice.
4. En la configuración de markers/CSV, usa:
      X column: {OUTPUT_X_COL}
      Y column: {OUTPUT_Y_COL}
      Color/group column: phenotype
5. Si los puntos no se alinean con la imagen, revisa MICRONS_PER_PIXEL_X/Y en el script.
   Tus columnas originales son '{X_UM_COL}' y '{Y_UM_COL}', en micras.
   Si TissUUmaps espera píxeles, necesitas convertir micras a píxeles usando el pixel size.

Nota:
- Cada CSV generado corresponde a un valor único de la columna '{IMAGE_COL}'.
- Si '{IMAGE_COL}' contiene escenas, por ejemplo '... - Scene #1', se genera un CSV por escena.
""".strip()

    instructions_path = output_dir / "README_how_to_open_in_TissUUmaps.txt"
    instructions_path.write_text(instructions, encoding="utf-8")

    print("\nProceso completado.")
    print(f"CSV markers generados: {len(index_df):,}")
    print(f"Índice: {index_path}")
    print(f"README: {instructions_path}")


if __name__ == "__main__":
    main()