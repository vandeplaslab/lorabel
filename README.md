# LORABEL

**Thermal Background Reduction for Mid-Infrared Imaging by Low-Rank Background and Sparse Point Source Modelling**

![PyPI Version](https://img.shields.io/pypi/v/lorabel) ![License](https://img.shields.io/badge/license-BSD%203--Clause-blue)

LORABEL implements LOw-RAnk Background ELimination, a novel computational method to improve sensitivity in mid-infrared astronomical imaging by modelling and removing varying thermal background noise without classical chopping or nodding.

---

## Table of Contents

1. [Features](#features)
2. [Installation](#installation)
3. [Quick Start](#quick-start)
4. [Algorithm Overview](#algorithm-overview)
5. [Parameters & Tuning](#parameters--tuning)
6. [Examples](#examples)
7. [Performance & Benchmarks](#performance--benchmarks)
8. [Citation](#citation)
9. [License](#license)
10. [Authors & Acknowledgements](#authors--acknowledgements)

---

## Features

* **Low-Rank Background Modelling**: captures quasi-static and slowly varying background patterns
* **Sparse Point Source Extraction**: isolates astrophysical point sources via ℓ₁ regularization
* **Stable Principal Component Pursuit**: solves

  ```text
  minimize ‖B‖_* + θ‖C‖₁
  s.t. ‖A - B - C‖_F ≤ δ
  ```
* **No Nodding Required**: works on chop-only data, reducing operational overhead
* **Flexible**: demonstrated on ground-based VISIR and airborne SOFIA/FORCAST datasets

## Installation

Install the latest release from PyPI:

```bash
pip install lorabel
```

Or install directly from GitHub:

```bash
pip install git+https://github.com/vandeplaslab/lorabel.git
```

---

## Quick Start

```python
import numpy as np
from lorabel import LORABEL

# Load your chop-only time series: a 3D array (frames × height × width)
data = np.load('chop_series.npy')

# Initialize the model (θ and δ can be auto-estimated)
model = LORABEL(data)

# Decompose into low-rank background, sparse sources, and residual noise
background, sources, residual = model.decompose()

# Aggregate the sparse component to detect and photometer point sources
time_avg = sources.mean(axis=0)
# ... proceed with your favorite photometry routine
```

---

## Algorithm Overview

1. **Vectorize Frames**: build matrix $A\in\mathbb{R}^{mn	imes t}$ with each column as a vectorized chop subtraction frame.
2. **Model**: $A = B + C + D$ where:

   * $B$ is low-rank (background)
   * $C$ is sparse (point sources)
   * $D$ is small dense noise
3. **Optimization**: solve Stable Principal Component Pursuit:

   ```text
   minimize ‖B‖_* + θ‖C‖₁
   subject to ‖A - B - C‖_F ≤ δ
   ```
4. **Reconstruct**: reshape $B,C,D$ back to 3D cubes for further analysis.

---

## Parameters & Tuning

* **θ (theta)**: trade-off between sparsity and background rank. Default ≈ 1/√max(m,n).
* **δ (delta)**: noise tolerance, often estimated as √(mn)·σ (σ = noise σ).

Use `model.auto_tune()` or supply custom values:

```python
model = LORABEL(data, theta=1.1/np.sqrt(max(m, n)), delta=0.7*np.linalg.norm(A, 'fro'))
```

---

## Examples

See the [examples/](https://github.com/vandeplaslab/lorabel/tree/main/examples) folder for:

* Synthetic VISIR data with injected sources
* Real SOFIA/FORCAST observations
* Jupyter notebooks demonstrating detection, photometry, and parameter sweeps

---

## Performance & Benchmarks

* **VISIR (SNR < 5)**: up to 5× reduction in photometric error variability vs. chop-nod methods.
* **SOFIA (mid-IR)**: 20–100× decrease in mean background flux while preserving source flux.

Refer to [Fig. 6–9](https://github.com/vandeplaslab/lorabel#results) for full quantitative benchmarks.

---

## Citation

If you use LORABEL in your research, please cite:

> R.A.R. Moens, A.G.M. Pietrow, B. Brandl, & R. Van de Plas (2025), “Thermal Background Reduction for Mid-Infrared Imaging by Low-Rank Background and Sparse Point Source Modelling,” *Astronomy & Astrophysics*, DOI: xx.xxxx/XXXX

---

## License

This project is licensed under the **BSD 3-Clause License**. See [LICENSE](LICENSE) for details.

---

## Authors & Acknowledgements

**Authors**: R.A.R. Moens, A.G.M. Pietrow, B. Brandl, R. Van de Plas

This package leverages open-source tools: NumPy, Matplotlib, Astropy, Photutils.

Supported by: TU Delft Space Institute, DFG, ESO, NASA/SOFIA.

Contributions and issues are welcome on [GitHub](https://github.com/vandeplaslab/lorabel).
