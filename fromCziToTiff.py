#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fromCziToImageJLikeTiff_MosaicSafe_FAST_v2.py

Versión rápida + multicanal + compatible con ImageJ/QuPath.

Corrige dos problemas:
1) En mosaics con una sola scene, read_mosaic NO debe recibir S.
2) ImageJ hyperstack requiere ejes en orden compatible con TZCYXS.
   Por tanto:
      - si conservas Z: se guarda como ZCYX (no CZYX)
      - si colapsas Z: se guarda como CYX

Pensado para QuPath:
- conserva TODOS los canales
- acelera usando Z_MODE='first' por defecto
- permite bajar la resolución del mosaic
- escribe TIFF ImageJ-like con Labels de canal
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import json
import re
import xml.etree.ElementTree as ET
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd
import tifffile
from aicspylibczi import CziFile

# ============================ CONFIGURACIÓN ==================================
INPUT_DIR = "/mnt/imgserver/IA/Projects/2026/2026_2_10_lpoluha/images/sub"
OUTPUT_DIR = "/mnt/imgserver/IA/Projects/2026/2026_2_10_lpoluha/output/images"
RECURSIVE = True
OUTPUT_EXTENSION = ".tif"
SKIP_IF_EXISTS = True
WRITE_METADATA_JSON = True

# Más rápido que max_projection y mantiene multicanal
# Opciones: first | middle | max_projection | keep
Z_MODE = "keep"
T_MODE = "keep"           # first | keep
CHANNEL_MODE = "keep"      # keep | first
MOSAIC_SCALE_FACTOR = 0.5   # prueba 0.25 si sigue lento
TIFF_COMPRESSION = None     # None = más rápido escribiendo
BIGTIFF = True
OUTPUT_MODE = "uint8_rescale"   # original | uint16_rescale | uint8_rescale
RESCALE_LOW_PCT = 0.0
RESCALE_HIGH_PCT = 100.0
N_WORKERS = 2               # sube a 3-4 si hay RAM/disco suficientes

# =============================== UTILS =======================================
def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path

def safe_filename(value: object, max_len: int = 180) -> str:
    s = str(value).strip().replace("/", "_").replace("\\", "_")
    s = re.sub(r"[^A-Za-z0-9._\-]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return (s[:max_len].rstrip("_")) or "scene"

def list_czi_files(input_dir: Path, recursive: bool) -> List[Path]:
    return sorted(input_dir.rglob("*.czi")) if recursive else sorted(input_dir.glob("*.czi"))

def get_dims_shape_safe(czi: CziFile):
    try:
        ds = czi.get_dims_shape()
        if isinstance(ds, list) and ds:
            return ds
        if isinstance(ds, dict):
            return [ds]
    except Exception:
        pass
    return []

def get_scene_indices(czi: CziFile) -> List[int]:
    dims_shapes = get_dims_shape_safe(czi)
    for d in dims_shapes:
        if "S" in d:
            start, size = d["S"]
            return list(range(int(start), int(start) + int(size)))
    try:
        scenes = list(czi.scenes)
        if scenes:
            return list(range(len(scenes)))
    except Exception:
        pass
    return [0]

def get_scene_count(czi: CziFile) -> int:
    dims_shapes = get_dims_shape_safe(czi)
    for d in dims_shapes:
        if "S" in d:
            _, size = d["S"]
            return int(size)
    try:
        scenes = list(czi.scenes)
        if scenes:
            return len(scenes)
    except Exception:
        pass
    return 1

def is_mosaic_safe(czi: CziFile) -> bool:
    try:
        return bool(czi.is_mosaic())
    except Exception:
        return False

def try_get_metadata_xml(czi: CziFile) -> str:
    for attr in ["meta", "metadata", "raw_metadata"]:
        try:
            val = getattr(czi, attr)
            if isinstance(val, str) and val.strip():
                return val
        except Exception:
            pass
    for meth in ["get_metadata", "read_metadata"]:
        try:
            val = getattr(czi, meth)()
            if isinstance(val, str) and val.strip():
                return val
        except Exception:
            pass
    return ""

def extract_channel_names_from_xml(xml_text: str, n_channels: int) -> List[str]:
    if not xml_text:
        return [f"Ch{i+1}" for i in range(n_channels)]
    try:
        root = ET.fromstring(xml_text)
    except Exception:
        return [f"Ch{i+1}" for i in range(n_channels)]
    if root.tag.startswith("{"):
        ns_uri = root.tag.split("}")[0].strip("{")
        ns = {"z": ns_uri}
        elems = root.findall(".//z:Channel", ns)
    else:
        elems = root.findall(".//Channel")
    names = []
    for ch in elems:
        nm = ch.attrib.get("Name") or ch.attrib.get("Fluor") or ch.attrib.get("Id") or ch.attrib.get("ID") or ""
        names.append(str(nm).strip() or f"Ch{len(names)+1}")
    if not names:
        names = [f"Ch{i+1}" for i in range(n_channels)]
    if len(names) < n_channels:
        names.extend([f"Ch{i+1}" for i in range(len(names), n_channels)])
    return names[:n_channels]

def choose_czt_indices(d: Dict[str, Tuple[int, int]]) -> Tuple[List[int], List[int], List[int]]:
    c_start, c_size = d.get("C", (0, 1))
    z_start, z_size = d.get("Z", (0, 1))
    t_start, t_size = d.get("T", (0, 1))
    c_idxs = [int(c_start)] if CHANNEL_MODE == "first" else list(range(int(c_start), int(c_start) + int(c_size)))
    if Z_MODE == "middle" and int(z_size) > 1:
        z_idxs = [int(z_start) + int(z_size) // 2]
    elif Z_MODE == "first" or int(z_size) <= 1:
        z_idxs = [int(z_start)]
    elif Z_MODE in {"max_projection", "keep"}:
        z_idxs = list(range(int(z_start), int(z_start) + int(z_size)))
    else:
        raise ValueError(f"Z_MODE no soportado: {Z_MODE}")
    t_idxs = [int(t_start)] if T_MODE == "first" else list(range(int(t_start), int(t_start) + int(t_size)))
    return c_idxs, z_idxs, t_idxs

def rescale_channel(x: np.ndarray, out_mode: str) -> np.ndarray:
    x = np.asarray(x)
    lo = np.percentile(x, RESCALE_LOW_PCT) if RESCALE_LOW_PCT > 0 else np.min(x)
    hi = np.percentile(x, RESCALE_HIGH_PCT) if RESCALE_HIGH_PCT < 100 else np.max(x)
    if not np.isfinite(lo): lo = 0.0
    if not np.isfinite(hi) or hi <= lo: hi = lo + 1e-6
    y = np.clip((x.astype(np.float32) - float(lo)) / float(hi - lo), 0, 1)
    if out_mode == "uint16_rescale": return (y * 65535.0 + 0.5).astype(np.uint16)
    if out_mode == "uint8_rescale": return (y * 255.0 + 0.5).astype(np.uint8)
    return x

def convert_output_dtype(stack: np.ndarray) -> np.ndarray:
    if OUTPUT_MODE == "original":
        if stack.dtype in (np.uint8, np.uint16, np.float32, np.float64, np.int16, np.int32):
            return stack
        return stack.astype(np.float32)

    # Soporta stack CYX o ZCYX
    if stack.ndim == 3:
        channels = [rescale_channel(stack[i], OUTPUT_MODE) for i in range(stack.shape[0])]
        return np.stack(channels, axis=0)
    elif stack.ndim == 4:
        # ZCYX -> rescale por canal usando todos los Z
        out = []
        for z in range(stack.shape[0]):
            chs = [rescale_channel(stack[z, c], OUTPUT_MODE) for c in range(stack.shape[1])]
            out.append(np.stack(chs, axis=0))
        return np.stack(out, axis=0)
    else:
        return stack

def read_one_plane_regular(czi: CziFile, d: Dict[str, Tuple[int, int]], scene_index: int, c: int, z: int, t: int) -> np.ndarray:
    kwargs = {"S": int(scene_index)}
    if "C" in d: kwargs["C"] = int(c)
    if "Z" in d: kwargs["Z"] = int(z)
    if "T" in d: kwargs["T"] = int(t)
    result = czi.read_image(**kwargs)
    arr = result[0] if isinstance(result, tuple) and len(result) == 2 else result
    return np.squeeze(np.asarray(arr))

def read_one_plane_mosaic(czi: CziFile, d: Dict[str, Tuple[int, int]], scene_index: int, c: int, z: int, t: int, scene_count: int) -> np.ndarray:
    base_kwargs = {"C": int(c), "scale_factor": MOSAIC_SCALE_FACTOR}
    if "Z" in d: base_kwargs["Z"] = int(z)
    if "T" in d: base_kwargs["T"] = int(t)
    if scene_count == 1:
        arr = czi.read_mosaic(**base_kwargs)
        return np.squeeze(np.asarray(arr))
    try:
        arr = czi.read_mosaic(S=int(scene_index), **base_kwargs)
        return np.squeeze(np.asarray(arr))
    except Exception:
        arr = czi.read_mosaic(**base_kwargs)
        return np.squeeze(np.asarray(arr))

def read_scene_multichannel(czi_path: str, scene_index: int) -> Tuple[np.ndarray, List[str], Dict[str, Any]]:
    czi = CziFile(czi_path)
    dims_shapes = get_dims_shape_safe(czi)
    d = dims_shapes[0] if dims_shapes else {}
    c_idxs, z_idxs, t_idxs = choose_czt_indices(d)
    is_mosaic = is_mosaic_safe(czi)
    scene_count = get_scene_count(czi)
    print(scene_count)
    xml_text = try_get_metadata_xml(czi)
    n_channels_total = int(d.get("C", (0, len(c_idxs)))[1])
    channel_names_all = extract_channel_names_from_xml(xml_text, max(len(c_idxs), n_channels_total))
    channel_planes, channel_names = [], []
    for local_c_i, c in enumerate(c_idxs):
        z_planes = []
        for z in z_idxs:
            if is_mosaic:
                arr = read_one_plane_mosaic(czi, d, scene_index, c, z, t_idxs[0], scene_count)
            else:
                arr = read_one_plane_regular(czi, d, scene_index, c, z, t_idxs[0])
            z_planes.append(arr)
        if Z_MODE == "keep" and len(z_planes) > 1:
            ch = np.stack(z_planes, axis=0)   # ZYX
        elif len(z_planes) == 1:
            ch = z_planes[0]                  # YX
        else:
            ch = np.nanmax(np.stack(z_planes, axis=0), axis=0)  # YX
        channel_planes.append(ch)
        ch_name = channel_names_all[c] if c < len(channel_names_all) else f"Ch{local_c_i+1}"
        channel_names.append(ch_name)

    # IMPORTANTE: orden compatible con ImageJ = ZCYX si keep Z, o CYX si no
    if Z_MODE == "keep" and len(channel_planes) > 0 and np.asarray(channel_planes[0]).ndim == 3:
        # channel_planes: lista de ZYX -> apilar en CZYX y reordenar a ZCYX
        czyx = np.stack(channel_planes, axis=0)
        stack = np.moveaxis(czyx, 0, 1)  # ZCYX
        axes = "ZCYX"
    else:
        stack = np.stack(channel_planes, axis=0)   # CYX
        axes = "CYX"

    meta = {
        "read_method": "read_mosaic" if is_mosaic else "read_image",
        "czi_dims": str(getattr(czi, 'dims', '')),
        "dims_shape": str(dims_shapes),
        "is_mosaic": is_mosaic,
        "scene_count": scene_count,
        "channel_names": channel_names,
        "axes": axes,
    }
    return stack, channel_names, meta

def write_imagej_like_tiff(out_path: Path, stack: np.ndarray, channel_labels: Sequence[str], axes: str):
    metadata = {"axes": axes, "mode": "composite", "Labels": list(channel_labels)}
    tifffile.imwrite(str(out_path), stack, imagej=True, metadata=metadata, compression=TIFF_COMPRESSION, bigtiff=BIGTIFF)

def process_one_file(czi_path_str: str, output_dir_str: str, meta_dir_str: str) -> List[Dict[str, Any]]:
    czi_path = Path(czi_path_str)
    output_dir = Path(output_dir_str)
    meta_dir = Path(meta_dir_str)
    rows = []
    try:
        czi = CziFile(str(czi_path))
        scene_indices = get_scene_indices(czi)
    except Exception as e:
        return [{"czi_file": str(czi_path), "status": "error_opening_czi", "error": str(e)}]
    for scene_index in scene_indices:
        scene_label = f"Scene_{scene_index + 1}"
        out_name = f"{safe_filename(czi_path.stem)}__{scene_label}{OUTPUT_EXTENSION}"
        out_path = output_dir / out_name
        if SKIP_IF_EXISTS and out_path.exists():
            rows.append({"czi_file": str(czi_path), "scene_index_S": scene_index, "scene_label": scene_label, "output_tif": str(out_path), "status": "skipped_exists"})
            continue
        try:
            stack, channel_labels, meta = read_scene_multichannel(str(czi_path), int(scene_index))
            original_shape = tuple(np.asarray(stack).shape)
            original_dtype = str(np.asarray(stack).dtype)
            stack_out = convert_output_dtype(np.asarray(stack))
            final_shape = tuple(np.asarray(stack_out).shape)
            final_dtype = str(np.asarray(stack_out).dtype)
            write_imagej_like_tiff(out_path, stack_out, channel_labels, meta["axes"])
            row = {
                "czi_file": str(czi_path), "scene_index_S": scene_index, "scene_label": scene_label,
                "output_tif": str(out_path), "status": "ok", "read_method": meta["read_method"],
                "is_mosaic": meta["is_mosaic"], "scene_count": meta["scene_count"],
                "czi_dims": meta["czi_dims"], "dims_shape": meta["dims_shape"], "axes": meta["axes"],
                "channel_names": ";".join(channel_labels), "n_channels": len(channel_labels),
                "original_shape": str(original_shape), "original_dtype": original_dtype,
                "final_shape": str(final_shape), "final_dtype": final_dtype,
                "z_mode": Z_MODE, "t_mode": T_MODE, "channel_mode": CHANNEL_MODE,
                "mosaic_scale_factor": MOSAIC_SCALE_FACTOR, "output_mode": OUTPUT_MODE, "imagej_like": True,
            }
            rows.append(row)
            if WRITE_METADATA_JSON:
                json_path = meta_dir / f"{safe_filename(czi_path.stem)}__{scene_label}.json"
                json_path.write_text(json.dumps(row, indent=2, ensure_ascii=False), encoding="utf-8")
        except Exception as e:
            rows.append({"czi_file": str(czi_path), "scene_index_S": scene_index, "scene_label": scene_label, "output_tif": str(out_path), "status": "error", "error": str(e)})
    return rows

# ================================= MAIN ======================================
def main() -> None:
    input_dir = Path(INPUT_DIR)
    output_dir = ensure_dir(Path(OUTPUT_DIR))
    meta_dir = ensure_dir(output_dir / "metadata")
    if not input_dir.exists():
        raise FileNotFoundError(f"No existe INPUT_DIR: {input_dir}")
    czi_files = list_czi_files(input_dir, RECURSIVE)
    if not czi_files:
        raise RuntimeError(f"No se encontraron archivos .czi en: {input_dir}")
    print("====================================================")
    print("CZI -> TIFF MULTICANAL IMAGEJ (MOSAIC SAFE, RÁPIDO, AXES FIX)")
    print(f"Input : {input_dir}")
    print(f"Output: {output_dir}")
    print(f"Archivos CZI encontrados: {len(czi_files)}")
    print(f"Z_MODE={Z_MODE} | T_MODE={T_MODE} | CHANNEL_MODE={CHANNEL_MODE} | MOSAIC_SCALE_FACTOR={MOSAIC_SCALE_FACTOR} | N_WORKERS={N_WORKERS}")
    print("====================================================")
    rows = []
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        futures = {ex.submit(process_one_file, str(p), str(output_dir), str(meta_dir)): p for p in czi_files}
        for i, fut in enumerate(as_completed(futures), start=1):
            p = futures[fut]
            try:
                rr = fut.result()
                rows.extend(rr)
                ok = sum(r.get("status") == "ok" for r in rr)
                err = sum(r.get("status") == "error" for r in rr)
                print(f"[{i}/{len(czi_files)}] {p.name}: ok={ok}, error={err}")
            except Exception as e:
                print(f"[{i}/{len(czi_files)}] {p.name}: ERROR general -> {e}")
                rows.append({"czi_file": str(p), "status": "error_future", "error": str(e)})
    pd.DataFrame(rows).to_csv(output_dir / "czi_to_imagej_like_tif_index.csv", index=False)
    print("\n====================================================")
    print("CONVERSIÓN COMPLETADA")
    print(f"Resultados en: {output_dir}")
    print(f"Índice: {output_dir / 'czi_to_imagej_like_tif_index.csv'}")
    print("====================================================")

if __name__ == "__main__":
    main()