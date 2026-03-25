# ioNERDSS Paper Supplementary Repository

[![View ioNERDSS Repository](https://img.shields.io/badge/GitHub-ioNERDSS_Repo-181717?style=for-the-badge&logo=github)](https://github.com/JohnsonBiophysicsLab/ionerdss)

This repository contains simulation inputs, output data, intermediate processing artifacts, and visualization files used in support of the ioNERDSS paper. It is intended to serve as part of the supplementary information (SI) accompanying the paper and to provide a structured record of the example systems, coarse-grained models, and NERDSS simulation results discussed in the paper.

Each top-level directory corresponds to a specific test case, assembly system, or structure analyzed in the study. Depending on the system, these directories may include:

- `nerdss_files/`: NERDSS-ready input files such as `parms.inp`, molecular definition files (`*.mol`), archived PDB inputs, restart files, and simulation outputs in `DATA/`.
- `visualizations/`: coarse-grained structural representations, PyMOL sessions/scripts, and reference images used to inspect or illustrate the generated models.
- `ode_results/`: deterministic comparison results and associated plots for systems where ODE-based analysis was performed.
- `logs/`: pipeline or processing logs documenting model preparation and execution.
- notebooks or scripts: case-specific analysis notebooks and helper scripts used during model generation, inspection, or figure preparation.

Current example systems in this repository include:

- `8erq/`
- `8y7s/`
- `6bno/`
- `dodecahedron/`
- `proteosome/`

