# Universal and Cross-Cultural Variations in Audio-Motor Synchronization Between French and Indian Participants

## Overview

This repository contains the data and code used in the analysis for the article: XXX

## Repository Structure

The repository is organized into two main workflows: **MATLAB preprocessing** (raw data → relative phases) and **R analysis** (statistical and qualitative analysis).

```
.
├── Matlab/
│   ├── DAT/
│   │   ├── XP_French/
│   │   │   ├── FR01/
│   │   │   │   ├── FR01_sync_1.txt
│   │   │   │   ├── FR01_sync_2.txt
│   │   │   │   ├── FR01_sync_3.txt
│   │   │   │   ├── FR01_syncslow_1.txt
│   │   │   │   ├── FR01_syncslow_2.txt
│   │   │   │   ├── FR01_syncslow_3.txt
│   │   │   │   └── FR01_syncslow_4.txt
│   │   │   └── [Additional participants...]
│   │   └── XP_Indian/
│   │       └── [Participant folders with trials...]
│   ├── PRG/
│   └── RES/
│       └── tbl_relphases.txt
├── R/
│   ├── main.Rmd
│   ├── main.pdf
│   ├── DAT/
│   │   └── tbl_relphases.txt
│   ├── PRG/
│   ├── RES/
│   └── Rmd_utils/
```

## Folder Organization

### MATLAB Workflow (Preprocessing)

**Purpose:** Convert raw experimental data to relative phase values.

- **`DAT/`** – Raw data organized by experimental condition
	- `XP_French/` – French participant data (13 participants)
	- `XP_Indian/` – Indian participant data (15 participants)
	- Each participant folder contains trial files (e.g., `FR01_sync_1.txt`, `FR01_syncslow_1.txt`)
- **`PRG/`** – MATLAB functions and preprocessing scripts that handle data extraction, file processing, and relative phase computation
- **`RES/`** – Output results
	- `tbl_relphases.txt` – Processed relative phases table (generated upon execution)

**Execution:** Run `main.m` to execute the complete preprocessing pipeline. The script will call all necessary functions from the `PRG` folder and save results to `RES/`.

**Notes:**
- The `relphase` function was taken from the RelPhase toolbox, authored by Tjeerd Dijkstra
- The `dwell_time` function was initially written by Julien Lagarde and used in Zelic, G., Mottet, D., & Lagarde, J. (2012). Behavioral impact of unisensory and multisensory audio-tactile events: pros and cons for interlimb coordination in juggling. *PLoS One*, *7*(2), e32308.

### R Workflow (Analysis)

**Purpose:** Perform qualitative, quantitative, and statistical analyses on relative phase data.

- **`DAT/`** – Input data
	- `tbl_relphases.txt` – Relative phases table (copied from MATLAB `RES/` folder)
- **`PRG/`** – R analysis functions and utilities
- **`RES/`** – Output results (tables and figures generated during analysis)
- **`Rmd_utils/`** – Utility functions for R Markdown processing

**Execution:** Run `main.Rmd` to perform the analysis. This file contains documentation and executable code. The compiled PDF (`main.pdf`) shows the results but not the code; refer to `main.Rmd` for implementation details.

## Workflow

1. **MATLAB preprocessing:** Execute `Matlab/main.m` to process raw data and generate `tbl_relphases.txt`
2. **Transfer data:** Copy `Matlab/RES/tbl_relphases.txt` to `R/DAT/`
3. **R analysis:** Execute `R/main.Rmd` to perform statistical and qualitative analyses, producing tables and figures in `R/RES/`

## Data Description

### Raw Data Format

The raw `.txt` files contain finger movement (goniometer) and metronome (speakers) data. Both signals were sent to an acquisition card to avoid any delays. The files have the following structure:

| Column | Description                         |
| ------ | ----------------------------------- |
| 1      | Time (goniometer)                   |
| 2      | Goniometer values (finger movement) |
| 3      | Time (speakers)                     |
| 4      | Speakers values (metronome)         |

**Sampling frequency:** 5,000 Hz

**Note:** The two time columns are identical.

### Experimental Conditions

Each trial consists of multiple plateaus at different metronome frequencies, starting from 1.0 Hz and increasing by 0.3 Hz every 15 stimuli.

- **Sync trials** – 1.0 Hz to 6.1 Hz
- **Syncslow trials** – 1.0 Hz to 3.7 Hz

### Processed Data

The `tbl_relphases.txt` file contains relative phase values computed from the raw goniometer and speaker data. It serves as input for statistical and qualitative analysis in R. The file has the following structure:

| Variable | Description |
|----------|-------------|
| `group` | French or Indian |
| `subject` | Participant ID |
| `task` | sync or syncslow |
| `task_number` | Trial number (1 to 4) |
| `frequency` | Metronome frequency for the current plateau (1.0 to 6.1 Hz in 0.3 Hz steps) |
| `tap` | Tap number within the current plateau |
| `DT` | Dwell time for the current plateau |
| `period` | Period between the current tap and the previous one |
| `relphase` | Relative phase value for the current tap |

## Dependencies

### Matlab

The codes were written using Matlab version R2021b

### R

The codes were written using R version 4.3.2 on RStudio version 2024.12.1. The required packages and their version are detailed in the main.Rmd.

## Citation

If you use this code or data, please cite:

XXX
