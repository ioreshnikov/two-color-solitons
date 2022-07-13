# Cherenkov Radiation and Scattering of External Dispersive Waves by Two-Color Solitons

This repository contains the source code necessary to reproduce the results
published in I. Oreshnikov et al. Cherenkov Radiation and Scattering of
External Dispersive Waves by Two-Color Solitons [arxiv](https://arxiv.org/abs/2207.03541).

## Layout

The source codes are split between two directories — `code` and `text`.

The first directory, `code`, contains all the simulation and plotting codes,
from minor helper functions to publication-ready plotters. There you will find
two sub-directories: `scripts` and `common`. `scripts` holds all the concrete
scripts that are necessary to produce the paper and some additional exploratory
tools. The detailed description of each script can be found either by reading
the documentation block in the script itself or by running the script with
`--help` argument. `common` is the directory that contains the code shared
between the scripts, most importantly the solver (inside `common/solver.py`)
and the fiber parameters setup (in `common/fiber.py`).

The second directory, `text`, contains TeX source of the paper and the figures.

## Requirements

The code is written in `python`, version `3.5` is the oldest one that we have
tried. To execute the project `Makefile` you need `make` installed on your
machine.

## Installation

We strongly advise you run an installation in a virtual environment. You can do
this by running

    $ python -m venv venv

from the root of the project. To activate the environment you execute

    $ . venv/bin/activate

Then from inside the environment you can install the package by running

    $ python -m pip install -e .

## Reproducing the results

We have tried to make the result reproduction as straightforward as possible.
All the actions necessary to build the paper from scratch are defined in the
Makefile at the root of the project. To run all the simulations and compile the
draft of the paper you activate the environment from the root of the project

    $ . venv/bin/activate

and then run

    $ make draft

This will produce the draft of the paper and put it inside `text` directory.
