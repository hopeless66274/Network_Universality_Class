# Temporal Self-Similarity Reveals Percolation Universality Classes in Complex Networks

## Overview

Catastrophic fragmentation and structural transitions are ubiquitous in real-world complex systems. However, identifying their universality classes remains a long-standing challenge due to strong structural heterogeneity and the absence of well-defined critical thresholds.

In this work, we uncover a robust phenomenon of **temporal self-similarity** governing dynamic percolation processes across diverse complex networks. By tracking the full statistics of incremental growth (or fragmentation) events, we show that the dynamics are characterized by two independent Fisher-type critical exponents:

- **τ₁**: governing the statistics of primary gap events  
- **τ₂**: governing higher-order gap statistics  

These two exponents uniquely determine the universality class of the system. All other standard critical exponents can be derived through newly established scaling relations.

After validating the framework on canonical network models, we apply it to extensive empirical datasets. Remarkably, real-world systems (biological, social, infrastructural) systematically exhibit universality classes that deviate from those predicted by idealized models, highlighting the role of higher-order structural correlations.

This repository provides the **simulation codes and data analysis pipelines** used in the paper:

> **“Temporal self-similarity reveals percolation universality classes in complex networks”**

---

## Key Features

- Parameter-free classification of universality classes  
- No need for prior knowledge of critical thresholds  
- Scalable to large real-world networks  
- Applicable to heterogeneous and empirical systems  
- Full pipeline: simulation → data → plotting → scaling analysis  

---

## Repository Structure

## Repository Structure

```text
.
├── prog/           Fortran simulation codes
└── Data_Plot/      Data and plotting scripts for all figures
    ├── Fig1/       Python scripts and data for Fig.1
    ├── Fig2/       Python scripts and data for Fig.2
    ├── Fig3/       Python scripts and data for Fig.3
    ├── Fig4/       Python scripts and data for Fig.4
    ├── Fig5/       Python scripts and data for Fig.5
    ├── FigS1/      Supplementary Fig.S1
    ├── FigS2/      Supplementary Fig.S2
    ├── FigS3/      Supplementary Fig.S3
    ├── FigS4/      Supplementary Fig.S4
    ├── FigS5/      Supplementary Fig.S5
    └── FigS6/      Supplementary Fig.S6

```

---

## Installation & Requirements

### Fortran (simulation)

- Compiler: `gfortran` (recommended)

Compile example:
```bash
gfortran Per2DV0.f90 -o Per2D
gfortran wjob.f -o wjob
./wjob
chmod +x job001
./job001 &

python >= 3.8
numpy
scipy
matplotlib
pandas




