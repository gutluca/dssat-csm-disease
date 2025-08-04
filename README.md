# Coupling of an Epidemiological Model for Leaf Fungal Diseases with the CROPGRO-Soybean Crop Simulation Model (DSSAT)

**Authors**: Gustavo de Angelo Luca, Izael Martins Fattori Jr., Fábio Ricardo Marin
**Affiliation**: University of São Paulo – ESALQ/USP  
**Version Date**: August 4, 2025

---

## Abstract

This repository contains a modified version of the [DSSAT](https://dssat.net/) (Decision Support System for Agrotechnology Transfer) source code, integrating a **generic epidemiological model** for foliar diseases into the **CROPGRO-Soybean** module.

The new module, `DISEASE_LEAF`, simulates the impact of leaf diseases—such as **Asian Soybean Rust** (*Phakopsora pachyrhizi*)—on soybean development and yield. It calculates daily diseased leaf area based on environmental variables and pathogen dynamics, dynamically adjusting **photosynthesis** and **leaf senescence** within CROPGRO.

---

## Methodology

### 1. The Epidemiological Model (`DISEASE_LEAF`)

A daily time-step state-transition model simulating pathogen dynamics:

- **Environmental Inputs**: T_min, T_max, and Relative Humidity (RH) from DSSAT weather files.
- **Key Processes**:
  - **Spore Deposition (F_DS)**
  - **Infection Rate (F_IR)** — via a beta response curve based on T and Leaf Wetness Duration (LWD)
  - **Latent Period (F_LR, F_LS)** — temperature-driven
  - **Spore Production (F_PS, F_PPSR)** — using a Verhulst logistic model and lesion aging (F_LAF)

#### Output:
- `DISEASE_LAI`: Daily diseased leaf area (fraction of total LAI)

---

### 2. Fungicide Control Module

Includes a rule-based system to simulate fungicide applications.

- **DVIP Calculation**: Daily value of infection probability (0–3), based on favorability matrices.
- **Application Rule**: Triggered if `SUM7(DVIP) ≥ 6` and `HEALTH_LAI ≥ 1.5`
- **Fungicide Effect**: Reduces infection rate based on:
  - `FUNG_EFFICIENCY`
  - `ResidualDays` (residual protection period)

---

### 3. Coupling with CROPGRO-Soybean

**Integration Points** in `CROPGRO.for`:

- **Declaration**: `DISEASE_LEAF` as `EXTERNAL`
- **Daily Call**: Inserted under `DYNAMIC .EQ. RATE`
- **Data Exchange**:
  - *To* DISEASE_LEAF: `TMIN`, `TMAX`, `RHUM`, `XLAI`
  - *From* DISEASE_LEAF: `DISLA` → used by the `PEST` subroutine

---

## Code Structure

```bash
dssat-csm-os-develop/
├── DISEASE.for               # New disease module (DISEASE_LEAF)
├── CROPGRO.for               # Modified to integrate disease model
├── disease_parameters.txt    # Input file with pathogen-specific parameters
├── DISEASE_DEVELOPMENT.OUT   # Output file with daily disease progress
```

---

## How to Use

### 1. Compilation

- Add `DISEASE.for` to your DSSAT-CSM project
- Replace original `CROPGRO.for` with the modified one
- Recompile DSSAT-CSM

### 2. Input Files

- Create `disease_parameters.txt` in `C:\DSSAT48\Soybean\`
- Ensure format matches the expectations of `READ_DISEASE_PARAMETERS`

### 3. Experiment Activation

In your `.X` experiment file, activate the disease module using the proper switch.

### 4. Output Files

`DISEASE_DEVELOPMENT.OUT` provides a daily log:

| Column         | Description                              |
|----------------|------------------------------------------|
| RUN            | Simulation run ID                        |
| YYDOY          | Year and day of year                     |
| DAS            | Days after sowing                        |
| LAI_HEALTH     | Healthy leaf area                        |
| LA_DISEASE     | Diseased leaf area                       |
| LWD            | Leaf Wetness Duration                    |
| RH             | Relative Humidity                        |
| FAT_TEMP       | Favorable Temperature                    |
| LAI_TOTAL      | Total Leaf Area Index                    |
| SUM7           | 7-day sum of DVIP                        |
| NSPRAYS        | Number of fungicide applications         |
| FUNG_ACT       | Indicator of active fungicide effect     |

---

## Citation

If you use this repository in your research, please cite:

> Luca, G. A.; Fattori Jr., I. M.; Marin, F.R (2025). *Coupling of an Epidemiological Model for Leaf Diseases with the CROPGRO-Soybean Crop Simulation Model (DSSAT)*. GitHub. [https://github.com/mewing3/dssat-csm-disease.git]

---

## References

- **Fungicide Management Logic**:
  - Godoy, C. V., et al. (2024). *Embrapa Technical Circular* (Hypothetical)
  - Beruski, N. D., et al. (2020). *Residual Period of Fungicide Protection* (Hypothetical)

- **Epidemiological Model Basis**:
  - Caubel, J., et al. (2017). *European Journal of Agronomy*, 90, 53–66.

- **DSSAT Reference**:
  - Jones, J. W., et al. (2003). *The DSSAT Cropping System Model*. *European Journal of Agronomy*, 18(3–4), 235–265.

---

