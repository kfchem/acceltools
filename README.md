# ACCeL Tools

**ACCeL Tools** (package name: `acceltools`) is an extension toolkit for [ACCeL](https://github.com/kfchem/accel), offering enhanced functionalities for molecular structure analysis and computational chemistry data handling. It provides utilities for document generation, ECD/UV/NMR analysis, plotting, puckering, and structural scanning, all based on molecular data managed by ACCeL.

---

## 📦 Features

- Export molecular data to TXT / DOCX formats
- Calculate and visualize ECD and UV spectra
- Assign and compare NMR chemical shifts
- Generate energy diagrams and scatter plots
- Analyze Cremer–Pople puckering parameters
- Perform coordinate scanning (bond/angle/dihedral)
- Structured data output via `pandas.DataFrame`

---

## 📄 Installation

```bash
pip install acceltools
```

> ⚠ Requires [ACCeL](https://github.com/kfchem/accel) to be installed and available in your Python environment.

To enable all optional features (e.g., `.docx` export and `pandas`-based data handling), use:

```bash
pip install acceltools[all]
```

---

## 📚 Module Overview

### 🧰 `base.py` — ToolBox Core

- `ToolBox` serves as the base class for all `acceltools` modules.
- Provides unified access to molecular data with filtering via `.get()`.

---

### 📄 `doc.py` — Document Export

- `DocBox` enables customizable formatting via `LayerAbc` subclasses.
- Available layers include `Name`, `Text`, `Data`, `Energy`, `Cartesian`, etc.
- Supports export to `.txt` and `.docx` via `export_txt()` and `export_docx()`.

---

### 🌈 `ecd.py` — ECD & UV Spectrum

- `EcdBox` provides:

  - Averaging across conformers
  - Curve generation using Gaussian broadening
  - Export of spectrum images as PNG
  - Overlay with experimental data

- Spectral data is stored in `System.data`.

---

### 🧪 `nmr.py` — NMR Shift Assignment

- `NmrBox` provides:

  - Import of tensors and experimental shifts from CSV
  - Conversion to chemical shifts via reference values
  - Assignment of calculated shifts to peaks
  - Statistical evaluation (MAE, RMSE, Max Error)

---

### 📊 `plot.py` — Plotting Tools

- `PlotBox` enables:

  - Reaction energy diagrams via `diagram()`
  - Custom scatter plots via `scatter()`
  - Export to `.png` and `.svg`

- Axis values and colors can be mapped via data keys.

---

### 🌀 `pucker.py` — Puckering Analysis

- `PuckerBox` calculates Cremer–Pople puckering parameters (Q, θ, φ).
- Includes conformer type matching (e.g., 4C1, 1C4, B25, etc.).

---

### 🔁 `series.py` — Geometry Series

- `SeriesBox` allows:

  - Generation of coordinate series by modifying bond lengths, angles, or dihedrals
  - Batch analysis of geometric properties
  - Export to `.xyz` trajectory via `write_trjxyz()`

---

### 📋 `table.py` — DataFrame Conversion

- `TableBox` extracts system-level data as `pandas.DataFrame`.
- Use `get_df(["energy", "label", "multiplicity", ...])` to specify columns.

---

## 🗂️ Dependencies

- [ACCeL](https://github.com/kfchem/accel)
- Python ≥ 3.9
- `scipy`, `matplotlib`
- Optional: `python-docx`, `pandas`

---

## 🧪 Testing & Examples

Example workflows and test cases can be found in the [`examples/`](https://github.com/kfchem/acceltools/tree/main/examples) directory.

---

## 📄 License

This project is licensed under the terms of the **MIT License**.

---

## 🗂️ Project Information

- **Repository**: [https://github.com/kfchem/acceltools](https://github.com/kfchem/acceltools)
- **Related Project**: [ACCeL](https://github.com/kfchem/accel)
- **Author**: Keisuke Fukaya
- **License**: MIT
- **Python Compatibility**: Python 3.9+

---

## 📬 Contact

For bug reports or feature requests, please open an issue on the [GitHub repository](https://github.com/kfchem/acceltools/issues).

For academic inquiries or collaboration requests, please contact the author via email: [kfukaya@pu-toyama.ac.jp](mailto:kfukaya@pu-toyama.ac.jp)

Thank you for your interest in ACCeL and acceltools!
