# Cherenkov Radiation and Scattering of External Dispersive Waves by Two-Color Solitons

This repository contains the source code necessary to reproduce the results published in I. Oreshnikov, A.V. Yulin. Cherenkov Radiation and Scattering of External Dispersive Waves by Two-Color Solitons [arxiv](https://arxiv.org/).

## Layout

The source codes are split between two directories — `code` and `text`.

The first directory, `code`, contains all the simulation and plotting codes, from minor helper functions to publication-ready plotters. There you will find two sub-directories: `scripts` and `common`. `scripts` holds all the concrete scripts that are necessary to produce the paper and some additional exploratory tools. The detailed description of each script can be found either by reading the documentation block in the script itself or by running the script with `--help` argument. `common` is the directory that contains the code shared between the scripts, most importantly the solver (inside `common/solver.py`) and the fiber parameters setup (in `common/fiber.py`).

The second directory, `text`, contains TeX source of the paper and the figures.

## Requirements

The code is written in `python`, version `3.5` is the oldest one that we have tried. To execute the project `Makefile` you need `make` installed on your machine.

## Reproducing the results

We have tried to make the result reproduction as straightforward as possible. All the actions necessary to build the paper from scratch — including the environment setup — are defined in the makefile at the root of the project. The makefile defines the following targets:

1. `venv` — sets up the virtual environment necessary to run the simulations.
2. `npz` — runs all the simulations and produces the raw data files. The simulations themselves take from 30 minutes to 3 hours to run on your average machine and running them sequentially might take _days_. However, if you have a machine powerful enough (with memory size being the most valuable resource), you can run it in parallel with `make -j8 npz` and bring the computation time down to several hours.
3. `fig` — compiles all the plots used in the paper.
4. `pdf` — compiles the pdf file of the paper.
5. `zip` — packages the TeX source and the figures in a zip archive suitable for a journal submission.

As always, you can run all the targets by executing

    make all

## Acknowledgment

Parts of the code are heavily inspired by the work of O. Melchert that he released as a Code Ocean instance https://codeocean.com/capsule/0624342/.
