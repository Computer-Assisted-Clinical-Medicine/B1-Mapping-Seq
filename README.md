# B1-Mapping-Seq — Sodium B1 Mapping (Pulseq)

This repository is part of an abstract submitted to ISMRM 2026. It provides
open-source implementations of two 3D Cartesian B1 mapping approaches for
sodium (²³Na) nuclei, implemented in MATLAB using the
[Pulseq](https://pulseq.github.io/) framework, together with corresponding
reconstruction scripts.

> **Disclaimer:** The accelerated acquisition using block undersampling has
> been tested extensively for the **Bloch-Siegert Shift method only**. TGV
> regularisation parameters (lambda, mu) are provided in `BSS_TGV_vals.mat`
> exclusively for the block undersampling pattern. The Double Flip Angle
> implementation supports undersampling pattern flags but has not been
> validated in the accelerated setting. Use non-block patterns and the DFA
> acceleration with caution.

---

## Overview

Two complementary B1 mapping methods are implemented:

- **Double Flip Angle (DFA)** — two acquisitions at flip angles
  α and 2α; B1 efficiency recovered from the signal magnitude ratio
- **Bloch-Siegert Shift (BSS)** — a 90° excitation followed by an
  off-resonance Fermi pulse at ±Δf; B1 amplitude recovered from the
  phase difference between the two acquisitions via TGV-regularised
  reconstruction (Lesch et al.)

Both sequences support configurable k-space undersampling (block or
Gaussian patterns) and write all acquisition parameters directly into the
`.seq` file for self-contained reconstruction.

---

## Dependencies

| Dependency | Purpose |
|---|---|
| [Pulseq for MATLAB](https://github.com/pulseq/pulseq) | Sequence design and `.seq` file export |
| [mapVBVD](https://github.com/CIC-methods/FID-A) | Raw data loading (Siemens TWIX format) |
| [BSReconFramework](https://github.com/IMTtugraz/BSReconFramework) | TGV-regularised Bloch-Siegert reconstruction (Lesch et al.) |

---

## Usage

### 1. Generate the sequence file

Open `B1map_DFA.m` or `B1map_BS.m` and adapt:

- **System limits** (`MaxGrad`, `MaxSlew`, RF/ADC dead times) to your scanner
- **Sequence parameters** (`fov`, `Nx/Ny/Nz`, `TR`, `alpha`)
- **Undersampling pattern** (`pattern`, `len_x`, `len_y`)
- **Output path** (`savestr`)

Run the script. It will perform a Pulseq timing check and, if passed, write
a `.seq` file to the specified path.

> ⚠️ **Safety notice:** Always verify SAR, gradient amplitude, slew rate, and
> duty cycle limits on your specific scanner before executing any sequence on
> hardware. This code is provided "as is" without warranty of any kind.

### 2. Reconstruct

**Double Flip Angle:**
```matlab
B1map = compB1map_DFA('path/to/B1map_DFA', false);   % relative B1
B1map = compB1map_DFA('path/to/B1map_DFA', true);    % absolute B1 [T]
```

**Bloch-Siegert (TGV):**
```matlab
B1Map = compB1map_BSS('path/to/B1map_BS');
B1Map = compB1map_BSS('path/to/B1map_BS', 64);  % optional: resample to 64³
```

Both functions read all acquisition parameters directly from the `.seq` file
and display three orthogonal central slices of the reconstructed B1 map.

---

## Parameters

### Double Flip Angle

| Parameter | Description | Default |
|---|---|---|
| `alpha1` / `alpha2` | Flip angles [degrees] | `60° / 120°` |
| `fov` | Field of view [m] | `[240e-3, 240e-3, 240e-3]` |
| `Nx / Ny / Nz` | Matrix size | `40 × 40 × 40` |
| `TR` | Repetition time | `250 ms` |
| `pattern` | Undersampling pattern | `'full'` |

### Bloch-Siegert Shift

| Parameter | Description | Default |
|---|---|---|
| `alpha` | Excitation flip angle [degrees] | `90°` |
| `fov` | Field of view [m] | `[240e-3, 240e-3, 240e-3]` |
| `Nx / Ny / Nz` | Matrix size | `40 × 40 × 40` |
| `TR` | Repetition time | `250 ms` |
| `T` | Fermi pulse duration [s] | `2 ms` |
| `bs_offset` | Off-resonance frequency [Hz] | `1000 Hz` |
| `pattern` | Undersampling pattern | `'block'` |

---

## References

1. **Sacolick, L.I. et al. (2010).**
   B1 mapping by Bloch-Siegert shift.
   *Magnetic Resonance in Medicine*, 63(5), 1315–1322.
   https://doi.org/10.1002/mrm.22357

2. **Lesch, A. et al. (2019).**
   Accelerated 3D Bloch-Siegert B1+ mapping using variational spatial
   regularization.
   *Magnetic Resonance in Medicine*, 81(2), 839–851.
   https://doi.org/10.1002/mrm.27434
   — Implementation: https://github.com/IMTtugraz/BSReconFramework


---

## License

Copyright (C) 2025 Valentin Jost

This program is free software: you can redistribute it and/or modify it under
the terms of the **GNU General Public License version 3** as published by the
Free Software Foundation.

This program is distributed in the hope that it will be useful, but **WITHOUT
ANY WARRANTY**; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.en.html)
for more details.

**Note on dependencies:** [Pulseq](https://github.com/pulseq/pulseq) is
licensed under MIT and [BSReconFramework](https://github.com/IMTtugraz/BSReconFramework)
is licensed under MIT, both of which are compatible with GPL-3.0. Any
derivative works that incorporate or modify code from this repository must
also be released under GPL-3.0.
