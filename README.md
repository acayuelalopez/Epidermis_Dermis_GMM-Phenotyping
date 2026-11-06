# Epidermis_Dermis_GMM-Phenotyping

This workflow converts QuPath single-cell intensity measurements into biologically defined phenotypes, generates QuPath-compatible object classifiers from the calculated thresholds, and prepares per-image marker tables for spatial visualization in TissUUmaps.

The workflow has three consecutive stages:

1. **Single-cell phenotyping** from QuPath intensity exports.
2. **QuPath object-classifier generation** from the calculated marker thresholds.
3. **TissUUmaps export** of phenotyped cells and their spatial coordinates.

## Workflow scripts

```text
.
├── fromQuPathIntensities_phenotyping_hierarchical_regions_threshold.py
├── generate_qupath_gmm_object_classifiers_channels.py
├── tissuumaps_generator.py
├── populations.txt
├── metadata.xlsx
└── README.md
```

## Workflow overview

```text
QuPath single-cell CSV files
          │
          ├── metadata.xlsx
          └── populations.txt
                  │
                  ▼
fromQuPathIntensities_phenotyping_hierarchical_regions_threshold.py
                  │
                  ├── phenotyped single-cell CSV files
                  ├── thresholds_by_marker_all_regions.csv
                  ├── phenotype reports
                  └── threshold and spatial plots
                         │
             ┌───────────┴───────────┐
             ▼                       ▼
generate_qupath_...py       tissuumaps_generator.py
             │                       │
             ▼                       ▼
QuPath classifier JSONs     per-image TissUUmaps marker CSVs
```

---

# 1. Requirements

## Input generated in QuPath

The phenotyping script expects single-cell measurement tables exported from QuPath. By default, it searches for:

```text
*_single_cell.csv
*_single_cell.tsv
*_single_cell.txt
```

The files must contain, or be configured to use, the following columns:

```text
Image
Classification
Centroid X µm
Centroid Y µm
```

They must also contain marker intensity measurements. The current workflow supports:

```text
S100B
K14
MDK
PTN
DAPI
```

The expected QuPath measurement format is:

```text
<MARKER>: <COMPARTMENT>: Mean
```

Examples:

```text
S100B: Cell: Mean
K14: Cell: Mean
MDK: Cell: Mean
PTN: Cell: Mean
```

The phenotyping script can use the following compartments independently for each marker:

```text
Cell
Nucleus
Cytoplasm
Membrane
```

## Python version and dependencies

Recommended version:

```text
Python 3.10 or newer
```

Install the required packages:

```bash
python -m pip install numpy pandas matplotlib scikit-learn openpyxl
```

The TissUUmaps export script uses NumPy and pandas. The phenotyping script additionally requires scikit-learn for Gaussian mixture models and silhouette analysis, and openpyxl for reading the metadata workbook.

---

# 2. Input files

## 2.1 QuPath single-cell files

Place all QuPath single-cell exports in one input directory. Each row must correspond to one cell or detection.

Example:

```text
measurements_cellpose_process/
├── 2025_11_18__0309_scene_1_single_cell.csv
├── 2025_11_18__0309_scene_2_single_cell.csv
├── 2025_11_18__0310_scene_1_single_cell.csv
└── ...
```

The script detects comma-separated and tab-separated tables automatically.

## 2.2 Metadata workbook

The workflow requires an Excel workbook containing one row per biological image or sample. The image identifier is configured with:

```python
METADATA_IMAGE_COL = "Image"
```

The matching logic compares the metadata image identifier with both:

- the single-cell filename;
- the value in the QuPath `Image` column.

Extensions such as `.czi`, `.tif`, `.tiff`, `.ome.tif`, and scene suffixes are normalized during matching.

The current configuration uses metadata fields such as:

```text
Image
Ref
TYPE
DIAGNOSIS
COMPONENT_TYPE
SIMPLIFIED_NAME
```

`Ref` is currently used to calculate thresholds by group:

```python
GMM_GROUPING_COLUMN = "Ref"
```

Every single-cell file should match one unambiguous metadata row. Review the generated `file_metadata_matching.csv` before interpreting downstream results.

## 2.3 Population definition file

Create a plain-text file named `populations.txt`, with one phenotype per line.

Example:

```text
PTN+
S100B+PTN+MDK-
S100B+MDK+PTN+
K14+PTN+MDK-
K14+MDK+PTN+
S100B-K14-PTN+MDK-
S100B-K14-MDK+PTN+
```

Blank lines and lines beginning with `#` are ignored.

Each expression combines marker positivity and negativity:

- `S100B+` means S100B intensity is above its threshold.
- `MDK-` means MDK intensity is at or below its threshold.
- `S100B+MDK+PTN-` requires all three conditions simultaneously.

More specific populations have priority. Populations with the largest number of marker conditions are evaluated first, followed by their original order in `populations.txt`.

DAPI is available as an intensity marker but is excluded from phenotype definitions by default:

```python
EXCLUDED_PHENOTYPE_MARKERS = {"DAPI"}
```

---

# 3. Stage one: phenotyping QuPath single-cell exports

Script:

```text
fromQuPathIntensities_phenotyping_hierarchical_regions_threshold.py
```

## 3.1 Main configuration

Edit the configuration block at the beginning of the script:

```python
INPUT_DIR = "/path/to/qupath/single_cell_files"
OUTPUT_DIR = "/path/to/phenotyping_results"
METADATA_XLSX = "/path/to/metadata.xlsx"
METADATA_SHEET = None
METADATA_IMAGE_COL = "Image"
POPULATIONS_TXT = "/path/to/populations.txt"
```

The script is designed to be edited and run directly from PyCharm or another Python environment.

## 3.2 Region configuration

The region column is currently:

```python
PARENT_COL = "Classification"
```

Regions are detected automatically when:

```python
REGIONS = None
```

To restrict the analysis, provide an explicit list:

```python
REGIONS = [
    "Annotation (Lesion)",
    "Annotation (Lesion_dermis)",
    "Annotation (Lesion_epidermis)",
]
```

### Hierarchical regions

The workflow can combine a main region with its subregions:

```python
HIERARCHICAL_REGIONS = {
    "Annotation (Lesion)": [
        "Annotation (Lesion_dermis)",
        "Annotation (Lesion_epidermis)",
    ],
}
```

With this configuration:

- `Annotation (Lesion)` includes cells directly assigned to the lesion plus cells from lesion dermis and lesion epidermis;
- `Annotation (Lesion_dermis)` still receives its own independent result;
- `Annotation (Lesion_epidermis)` still receives its own independent result;
- the original `Classification` is retained for traceability.

## 3.3 Intensity compartments

The default compartment is:

```python
DEFAULT_INTENSITY_COMPARTMENT = "Cell"
```

Marker-specific overrides can be defined as follows:

```python
MARKER_INTENSITY_COMPARTMENTS = {
    "S100B": "Nucleus",
    "MDK": "Cytoplasm",
    "K14": "Membrane",
    "PTN": "Cell",
}
```

The selected QuPath measurement must exist in the single-cell table. The script prints the column selected for every marker and input file.

## 3.4 Threshold methods

Two threshold methods are supported.

### Single Distribution

The threshold is calculated on the raw intensity scale:

```text
threshold = mean + std_multiplier × sample standard deviation
```

The default configuration is:

```python
DEFAULT_THRESHOLD_METHOD = "Single Distribution"
DEFAULT_SINGLE_DISTRIBUTION_STD_MULTIPLIER = 1.0
```

Therefore, the default threshold is:

```text
mean + 1 × SD
```

### Gaussian mixture model

For markers configured as `GMM`, the workflow:

1. applies `log1p` transformation when configured;
2. evaluates models from `K_MIN = 2` to `K_MAX = 8` components;
3. calculates SSQ, silhouette score, and BIC;
4. identifies an elbow-based lower bound;
5. chooses a model using BIC and silhouette within the allowed elbow margin;
6. defines the positivity threshold as the mean of the fitted component means;
7. converts the threshold back to the raw intensity scale.

The current marker-specific configuration is:

```python
MARKER_THRESHOLD_SETTINGS = {
    "S100B": {"method": "GMM"},
    "MDK": {"method": "GMM"},
    "K14": {"method": "GMM"},
    "PTN": {
        "method": "Single Distribution",
        "std_multiplier": 1.0,
    },
}
```

Markers absent from this dictionary use the default threshold method.

## 3.5 Group-specific thresholds

Thresholds are calculated independently for every value of:

```python
GMM_GROUPING_COLUMN = "Ref"
```

A global `__ALL__` threshold is also calculated for every marker and region. During phenotype assignment, the group-specific threshold is used when available; otherwise, the script falls back to `__ALL__`.

Use a grouping variable that represents a technically meaningful batch, reference, acquisition condition, or staining group. Avoid grouping by a biological outcome if doing so would make case and control positivity definitions systematically different without a justified technical reason.

## 3.6 Running the script

After editing the configuration, run:

```bash
python fromQuPathIntensities_phenotyping_hierarchical_regions_threshold.py
```

## 3.7 Main phenotyping outputs

```text
phenotyping_results/
├── file_metadata_matching.csv
├── thresholds_by_marker_all_regions.csv
├── single_cell_phenotyped_outputs_index.csv
└── phenotyping/
    ├── Annotation_Lesion/
    │   ├── single_cell_phenotyped__Annotation_Lesion.csv
    │   ├── reports/
    │   │   ├── thresholds_by_marker.csv
    │   │   ├── phenotype_counts_normalized_by_annotation_detections_per_image.csv
    │   │   ├── phenotype_counts_normalized_by_annotation_detections_per_scene.csv
    │   │   ├── counts_by_image_phenotype_normalized.csv
    │   │   └── counts_by_scene_phenotype_normalized.csv
    │   ├── gmm_plots/
    │   │   ├── *_gmm_metrics.csv
    │   │   ├── *_gmm_evaluation.png
    │   │   ├── *_hist_gmm.png
    │   │   └── *_hist_single_distribution.png
    │   └── plots/
    │       ├── stacked_barplot_per_image.png
    │       ├── stacked_by_*.png
    │       ├── barplot_by_*.png
    │       └── spatial_by_scene/
    └── ...
```

### Important phenotyped single-cell columns

| Column | Meaning |
|---|---|
| `analysis_region` | Region analyzed, including hierarchical combinations. |
| `source_annotation_classification` | Original QuPath classification. |
| `<MARKER>__value` | Raw marker intensity selected from QuPath. |
| `<MARKER>_positive` | Boolean threshold result. |
| `<MARKER>_threshold_raw` | Threshold applied to that cell. |
| `phenotype` | Assigned population. |
| `population` | Copy of the assigned phenotype. |
| `gmm_group` | Metadata group used for threshold selection. |
| `annotation_detection_count_per_image` | Total detections in the analysis region per image. |
| `population_count_per_image_annotation` | Number of cells in the assigned phenotype. |
| `population_fraction_per_image_annotation` | Population count divided by total region detections. |

The reported population proportions are normalized by the number of detections in the corresponding annotation, not by physical annotation area.

---

# 4. Stage two: generating QuPath object classifiers

Script:

```text
generate_qupath_gmm_object_classifiers_channels.py
```

This stage converts the raw thresholds calculated in stage one into QuPath Object Classifier JSON files.

## 4.1 Configuration

Set the threshold table and output folder:

```python
THRESHOLDS_CSV = Path(
    "/path/to/phenotyping_results/thresholds_by_marker_all_regions.csv"
)

OUTPUT_DIR = Path(
    "/path/to/phenotyping_results/qupath_object_classifiers"
)
```

Optionally restrict the regions or threshold groups:

```python
GROUPS_TO_EXPORT = None
REGIONS_TO_EXPORT = None
```

Examples:

```python
GROUPS_TO_EXPORT = ["__ALL__"]
```

```python
REGIONS_TO_EXPORT = ["Annotation (Lesion)"]
```

## 4.2 Marker-to-channel mapping

This mapping must match the exact measurement names present in QuPath:

```python
MEASUREMENT_BY_MARKER = {
    "S100B": "AF751: Cell: Mean",
    "K14": "AF488: Cell: Mean",
    "MDK": "AF555: Cell: Mean",
    "PTN": "AF647: Cell: Mean",
    "DAPI": "AF405: Cell: Mean",
}
```

This mapping is intentionally separate from the marker names used in the phenotyping script. Verify it against the column names in the original QuPath measurement table before generating the classifiers.

If stage one used a compartment other than `Cell`, update the classifier measurement accordingly. For example:

```python
"S100B": "AF751: Nucleus: Mean"
```

## 4.3 Population definitions

The classifier generator currently uses a Python list:

```python
POPULATIONS = [
    "PTN+",
    "S100B+PTN+MDK-",
    "S100B+MDK+PTN+",
    "K14+PTN+MDK-",
    "K14+MDK+PTN+",
]
```

For consistency, this list should match the populations selected in `populations.txt`. A classifier is skipped if any required marker threshold is missing.

## 4.4 Running the classifier generator

```bash
python generate_qupath_gmm_object_classifiers_channels.py
```

## 4.5 Outputs

```text
qupath_object_classifiers/
├── manifest_qupath_object_classifiers.csv
├── Annotation_Lesion/
│   ├── __ALL__/
│   │   ├── PTN+.json
│   │   ├── S100B+PTN+MDK-.json
│   │   └── ...
│   └── <group>/
│       └── ...
└── ...
```

Each JSON is a QuPath `CompositeClassifier` containing one simple threshold classifier per marker condition.

The classification logic matches the Python phenotyping:

```text
positive: intensity > threshold
negative: intensity <= threshold
```

The manifest records:

- analysis region;
- threshold group;
- population;
- generation status;
- missing markers;
- JSON path.

## 4.6 Using the JSON classifiers in QuPath

Copy or import the generated JSON files into the QuPath project object-classifier location, then open them through QuPath's object-classifier interface. Apply only the classifier corresponding to the correct analysis region and threshold group.

Before batch application, test representative images and confirm that:

- measurement names resolve correctly;
- positive and negative objects match the Python results;
- the classifier is applied to the intended detection objects;
- the selected group-specific threshold is appropriate for the image.

---

# 5. Stage three: preparing TissUUmaps markers

Script:

```text
tissuumaps_generator.py
```

This stage converts a phenotyped single-cell table into one marker CSV per QuPath image or scene.

## 5.1 Important input requirement

`INPUT_SINGLE_CELL_CSV` must point to an actual phenotyped single-cell CSV, for example:

```text
single_cell_phenotyped__Annotation_Lesion.csv
```

Do not point it to:

```text
single_cell_phenotyped_outputs_index.csv
```

The index contains file paths, not one row per cell, and therefore does not contain the required centroid and phenotype columns.

## 5.2 Configuration

```python
INPUT_IMAGES_DIR = "/path/to/images"

INPUT_SINGLE_CELL_CSV = (
    "/path/to/phenotyping/Annotation_Lesion/"
    "single_cell_phenotyped__Annotation_Lesion.csv"
)

OUTPUT_DIR = "/path/to/tissuumaps_output"
```

Required columns are configured as:

```python
IMAGE_COL = "Image"
X_UM_COL = "Centroid X µm"
Y_UM_COL = "Centroid Y µm"
PHENOTYPE_COL = "phenotype"
```

## 5.3 Coordinate conversion

TissUUmaps generally requires marker coordinates in the same coordinate system as the displayed image.

If the centroids are in µm and the image is opened in pixel coordinates:

```text
x_pixels = Centroid X µm / microns_per_pixel_x
y_pixels = Centroid Y µm / microns_per_pixel_y
```

Configure the pixel size obtained from QuPath:

```python
MICRONS_PER_PIXEL_X = 0.3450
MICRONS_PER_PIXEL_Y = 0.3450
```

If the image and marker coordinates already use µm, set both values to:

```python
MICRONS_PER_PIXEL_X = 1.0
MICRONS_PER_PIXEL_Y = 1.0
```

An incorrect pixel size will produce spatially scaled or displaced points even when the phenotyping itself is correct.

## 5.4 Per-scene behavior

The script groups cells by each unique value in the `Image` column. Therefore, values such as:

```text
2025_11_18__0309.czi - Scene #1
2025_11_18__0309.czi - Scene #2
```

produce separate marker files.

The script attempts to match each group with an image in `INPUT_IMAGES_DIR`. CZI files may need to be converted to a format supported by the local TissUUmaps deployment, such as OME-TIFF or a pyramidal web-compatible image.

## 5.5 Running the exporter

```bash
python tissuumaps_generator.py
```

## 5.6 TissUUmaps outputs

```text
tissuumaps_output/
├── markers_csv/
│   ├── <image_or_scene>__tissuumaps_markers.csv
│   └── ...
├── metadata/
│   ├── <image_or_scene>__metadata.json
│   └── ...
├── phenotype_color_reference.tsv
├── tissuumaps_files_index.csv
└── README_how_to_open_in_TissUUmaps.txt
```

Each marker CSV includes:

```text
x
y
x_um
y_um
phenotype
image
cell_id
```

It also preserves configured metadata fields when they are available, including:

```text
Image Name
source_file
Parent
analysis_region
population
Ref
DIAGNOSIS
TYPE
COMPONENT_TYPE
population_fraction_per_image_annotation
```

## 5.7 Opening the results in TissUUmaps

For each image or scene:

1. Open the corresponding image in TissUUmaps.
2. Load the matching CSV listed in `tissuumaps_files_index.csv`.
3. Select `x` as the X coordinate.
4. Select `y` as the Y coordinate.
5. Use `phenotype` as the categorical grouping or color column.
6. Compare the displayed colors with `phenotype_color_reference.tsv` if a stable reference is needed.
7. Confirm alignment using recognizable tissue boundaries before interpreting spatial patterns.

---

# 6. Quality control

## Metadata matching

Inspect:

```text
file_metadata_matching.csv
```

Confirm that every `source_file` is assigned to the correct metadata image and `Ref` group. Incorrect matching changes both threshold estimation and downstream group plots.

## Threshold review

For every marker, region, and group, inspect:

- `thresholds_by_marker.csv`;
- GMM histogram plots;
- GMM evaluation plots;
- single-distribution histograms;
- number of cells used for fitting;
- `status` and `threshold_raw` values.

A threshold with `too_few_values`, `gmm_failed`, or a non-finite value should not be used to classify cells.

## Population validation

Check that:

- all intended populations are present in `populations.txt`;
- the classifier generator uses the same population definitions;
- broader populations do not unintentionally absorb cells before more specific populations;
- the `Negative` class is biologically and technically interpretable.

## Region validation

For hierarchical regions, compare:

```text
analysis_region
source_annotation_classification
analysis_region_members
```

This confirms whether a combined lesion result includes exactly the intended dermal and epidermal subregions.

## Spatial validation

Review the generated spatial plots and the TissUUmaps overlay. Confirm that:

- X and Y columns are correct;
- Y-axis orientation matches the displayed image;
- µm-to-pixel conversion uses the exact image calibration;
- image or scene matching is correct;
- points remain inside the expected tissue and annotation regions.

---

# 7. Common problems

## Marker intensity column not found

The console reports a warning such as:

```text
WARNING: S100B -> Cell: Mean column not found
```

Check the original QuPath export header and update:

```python
DEFAULT_INTENSITY_COMPARTMENT
MARKER_INTENSITY_COMPARTMENTS
```

Also verify marker names and channel names in QuPath.

## Metadata group column does not exist

If `Ref` is unavailable, change:

```python
GMM_GROUPING_COLUMN = "Ref"
```

to a valid metadata column. The selected grouping should correspond to a justified technical thresholding group.

## No cells found for a region

Check the exact strings present in:

```text
Classification
```

Then update `REGIONS` or `HIERARCHICAL_REGIONS`. Region matching is based on the full normalized string.

## Classifier JSON is generated but does not work in QuPath

The most frequent cause is an incorrect measurement name. Confirm that `MEASUREMENT_BY_MARKER` exactly matches the QuPath measurement table, including channel, compartment, and statistic.

## TissUUmaps reports missing columns

Confirm that `INPUT_SINGLE_CELL_CSV` points to a phenotyped single-cell file, not the output index. The input must contain:

```text
Image
Centroid X µm
Centroid Y µm
phenotype
```

## TissUUmaps points are scaled incorrectly

Check:

```python
MICRONS_PER_PIXEL_X
MICRONS_PER_PIXEL_Y
```

These values must match QuPath's pixel calibration for the image used in TissUUmaps.

## TissUUmaps cannot open CZI

Convert the scene or image to a compatible OME-TIFF or pyramidal format while preserving dimensions and pixel calibration. Ensure that the converted image has the same origin and coordinate system as the QuPath measurements.

---

# 8. Reproducibility recommendations

Record the following for each analysis:

- QuPath version.
- Segmentation method and parameters.
- QuPath channel names.
- Intensity compartment used for each marker.
- Pixel size in µm/pixel.
- Metadata workbook version.
- Threshold grouping variable.
- Threshold method for every marker.
- GMM transformation and `k` range.
- Single-distribution SD multiplier.
- Exact `populations.txt` file.
- Hierarchical-region configuration.
- Python and package versions.
- TissUUmaps version or deployment.


