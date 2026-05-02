# dos-inverse-models

Optical property recovery algorithms and data type conversions for Diffuse Optical Spectroscopy (DOS).

## Overview

This repository focuses on the "Inverse Problem" in DOS: recovering tissue optical properties (absorption $\mu_a$ and reduced scattering $\mu_s'$) from measured data. It includes iterative solvers, spectroscopy tools, and legacy blood flow models.

## Contents

### Examples (`examples/`)
- **`broadband_spectroscopy/`**: Supporting code for multi-wavelength tissue analysis (based on Blaney et al., *Appl. Opt.* 2021).
  - `exampleLookRawData.m`: Script for visualizing and analyzing broadband spectra.
  - `TissDataAPOPAI2021.mat`: Example dataset from multiple tissue types.

### Iterative Property Recovery (`iterOptPropRecov/`)
Non-linear optimization scripts to recover optical properties:
- **DSI/DSR:** Recovery from Dual-Slope Intensity and Dual-Slope Reflectance measurements (`DSI2mua_iterRecov.m`, `DSR2muamusp_iterRecov.m`).
- **MDR:** Recovery from Multi-Distance Reflectance measurements (`MDR2muamusp_iterRecov.m`).
- Supports both standard and enhanced (EB) measurement types.

### Data Type & Pathlength Utilities (`FDdataTypes/`)
Core conversions for Frequency-Domain (FD) data:
- `calcData_datTyp.m`: Converts raw FD data (AC/DC/Phase) into processed data types.
- `calcPathLen_datTyp.m`: Calculates differential pathlength factors (DPF) for various measurement domains.

### Spectroscopy & Chromophores (`Spectra/`)
Tools for multi-wavelength analysis:
- `makeE.m`: Generates extinction coefficient matrices for standard tissue chromophores (HbO, HbR, Water, Lipid, etc.).
- `extrapMUSP.m`: Extrapolates scattering properties across a wide wavelength range using power-law fitting.

### Legacy Blood Flow & Inverse Models (`abs_multiDist/legacy/` & `Muscle/legacy/`)
Refactored and legacy scripts for calculating hemodynamics:
- **Blood Flow:** Solutions for flow-related metrics (`I2Blood_DCslope.m`, `Ph2Blood.m`).
- **Muscle NIRS:** Legacy implementations for muscular oxygenation monitoring (`I2Blood_Phslope_3pi.m`).

### Advanced Solvers (`TwoLayer/`)
Inverse solvers specifically for layered tissue structures:
- `TwoLayer_InverseMarquardt.m`: Levenberg-Marquardt optimizer for two-layer media.
- `marquardt_2DE_old.m`: Modified legacy solvers for efficient convergence.

## Author
Developed by Giles Blaney, Ph.D.
