# LORABEL

**Thermal Background Reduction for Mid-Infrared Imaging by Low-Rank Background and Sparse Point Source Modelling**

![License](https://img.shields.io/badge/license-BSD%203--Clause-blue)

LORABEL implements LOw-RAnk Background ELimination, a novel computational method to improve sensitivity in mid-infrared astronomical imaging by modelling and removing varying thermal background noise without classical chopping or nodding.

---

## Table of Contents

1. [Features](#features)
2. [Installation](#installation)
3. [Quick Start](#quick-start)
4. [Algorithm Overview](#algorithm-overview)
5. [Parameters & Tuning](#parameters--tuning)
6. [Citation](#citation)
7. [License](#license)
8. [Authors & Acknowledgements](#authors--acknowledgements)

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

Install directly from GitHub:

```bash
pip install git+https://github.com/vandeplaslab/lorabel.git
```

---

## Quick Start

```python
from astropy.io import fits
from lorabel import LORABEL

# Load your chop subtraction time series: a 3D array (height × width x frames) and reshape to 2D array (pixels x frames)
hdul = fits.open('chop_series.fits')
data = hdul[0].data  # shape (height, width, frames)
hdul.close()
data = data.reshape(data.shape[0], -1)

# Initialize the model (θ and δ can be adjusted)
model = LORABEL()

# Decompose into low-rank background, sparse sources, and residual noise
background, sources, residual = model.decompose(data)

# Aggregate the sparse component to detect and photometer point sources
time_avg = sources.mean(axis=0)
# ... proceed with your favorite photometry routine
```

---

## Algorithm Overview

1. **Vectorize Frames**: build matrix $A\in\mathbb{R}^{mn \times t}$ with each column as a vectorized chop subtraction frame.
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

Supply custom values:

```python
model = LORABEL(data, theta=1.1/np.sqrt(max(m, n)), delta=0.7*np.linalg.norm(A, 'fro'))
```

---

## Citation

If you use LORABEL in your research, please cite (submitted to A&A):

> R.A.R. Moens, A.G.M. Pietrow, B. Brandl, & R. Van de Plas (2025), “Thermal Background Reduction for Mid-Infrared Imaging by Low-Rank Background and Sparse Point Source Modelling,” * *, DOI: xx.xxxx/XXXX

---

## License

This project is licensed under the **BSD 3-Clause License**. See [LICENSE](LICENSE) for details.

---

## Authors & Acknowledgements

**Authors**: R.A.R. Moens, A.G.M. Pietrow, B. Brandl, R. Van de Plas

This package leverages open-source tools: NumPy, Matplotlib, Astropy, Photutils.

Supported by: TU Delft Space Institute, DFG, ESO, NASA/SOFIA.

Contributions and issues are welcome on [GitHub](https://github.com/vandeplaslab/lorabel).
