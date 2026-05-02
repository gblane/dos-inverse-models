# dos-inverse-models

> [!CAUTION]
> **WORK IN PROGRESS**: This repository is currently being reorganized for public release. Documentation and examples are subject to change.

Optical property recovery algorithms and data type conversions for Diffuse Optical Spectroscopy (DOS).

## Overview

This repository focuses on the "Inverse Problem" in DOS: recovering tissue optical properties (absorption $\mu_a$ and reduced scattering $\mu_s'$) from measured data. It includes iterative solvers, spectroscopy tools, and specialized hemodynamic models.

## Repository Structure

### Source Code (`src/`)
- **`solvers/`**: Iterative property recovery algorithms for various measurement types (DSI, DSR, MDR) and non-linear optimizers for layered tissue.
- **`spectroscopy/`**: Tools for multi-wavelength analysis, including extinction coefficient calculation and scattering extrapolation.
- **`conversions/`**: Core logic for frequency-domain data-type processing and differential pathlength/slope factor (DPF/DSF) calculation.
- **`hemodynamics/`**: Legacy scripts and tissue-specific models for calculating blood flow and muscle oxygenation.

### Examples (`examples/`)
- **`broadband_spectroscopy/`**: Supporting code for multi-wavelength tissue analysis (based on Blaney et al., *Appl. Opt.* 2021).

### Shared Data (`data/`)
- Consolidated extinction spectra and multi-subject tissue datasets.

## Citations

If you use this toolkit in your research, please cite the relevant publications:

1.  **Dual-Slope Foundations:** Blaney, G., Sassaroli, A., Pham, T., Fernandez, C., & Fantini, S. (2019). Phase dual-slopes in frequency-domain near-infrared spectroscopy for enhanced sensitivity to brain tissue: First applications to human subjects. *Journal of Biophotonics*, 12(11), e201960018. [https://doi.org/10.1002/jbio.201960018](https://doi.org/10.1002/jbio.201960018)
2.  **Broadband Spectroscopy:** Blaney, G., Curtsmith, P., Sassaroli, A., Fernandez, C., & Fantini, S. (2021). Broadband absorption spectroscopy of heterogeneous biological tissue. *Applied Optics*, 60(25), 7552-7562. [https://doi.org/10.1364/AO.431013](https://doi.org/10.1364/AO.431013)

## Author
Developed by Giles Blaney, Ph.D.

---
*This repository is a reorganized and documented version of a personal codebase, performed by Gemini CLI.*
