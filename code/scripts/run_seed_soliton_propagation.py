#!/usr/bin/env python3


__doc__ = """
This scripts integrates GNLSE with the initial condition given as a
sum of two fundamental solitons

    u(0, t) = U₁(0, t) + U₂(0, t)

where

    Uₙ(0, t) = Aₙ sech(t/tₙ) exp(-i ωₙ t).

Carrier frequency of the first soliton is passed as parameter `-f1`
of the script and is equal to `1.0` by default. Soliton durations t₁
and t₂ are ad from parameters `-t1` and `-t2` respectively and are
set to 20 fs by default. The amplitudes are calculated based on the
durations and local dispersion coefficients in the fiber.

The computational grid is defined as follows

    t ∈ [-8000, 8000), 2¹⁶ equidistant points
    z ∈ [0, 10000], 1000 equidistant points

In physical units, that corresponds to propagating the initial
conditions for 10 cm of the fiber with a step-size of 100 μm.
While this is good enough if you interested purely in the dynamics,
it is not fitting if you want to analyze the spatial _spectrum_ of
the numerical solutions.
"""


import argparse
import logging

import numpy

from common.fiber import (
    beta, beta1, beta2,
    gamma, kerr_op)
from common.helpers import (
    sech, to_analytic,
    gv_matching_frequencies,
    fundamental_soliton_amplitude)
from common.solver import gnlse
from common.plotter import gpc_setup


logging.basicConfig(
    level=logging.DEBUG,
    format="%(process)s %(asctime)s %(levelname)s %(message)s")


parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "-f1",
    help="carrier frequency of the first soliton",
    type=float,
    default=1.000)
parser.add_argument(
    "-t1",
    help="width of the first soliton in fs",
    type=float,
    default=20.000)
parser.add_argument(
    "-t2",
    help="width of the first soliton in fs",
    type=float,
    default=20.000)
parser.add_argument(
    "output",
    help="path to the output .npz file",
    type=str)
args = parser.parse_args()


# Define the computational grid.
z = numpy.linspace(0, 1E5, int(1E3))
t = numpy.linspace(-8000, +8000, 2**16)

# Pick the carrier frequency of the first pulse at then find a
# frequency for the second pulse by picking a frequency in a region of
# anomalous dispersion so that the group velocity of both pulses match.
f1 = args.f1
f2 = gv_matching_frequencies(f1)[-1]

assert beta2(f1) < 0
assert beta2(f2) < 0


# Define the initial condition as a sum of two fundamental solitons of
# equal width at the selected frequencies.
t1 = args.t1
a1 = fundamental_soliton_amplitude(f1, t1)
u1 = a1 * sech(t/t1) * numpy.exp(-1j * f1 * t)

t2 = args.t2
a2 = fundamental_soliton_amplitude(f2, t2)
u2 = a2 * sech(t/t2) * numpy.exp(-1j * f2 * t)

u0 = to_analytic((u1 + u2).real)


# Construct a frequency filtered and group velocity compensated
# dispersive profile.
def filtered_beta(f):
    b = beta(f) - beta(f1) - beta1(f1) * (f - f1)
    b[(f <= 0) | (f >= 0.75 * f.max())] = 0
    return b


# Frequency filter the nonlinear coefficient.
def filtered_gamma(f):
    g = gamma(f)
    g[(f <= 0) | (f >= 0.75 * f.max())] = 0
    return g


# Setup interactive plotting
gpc = gpc_setup(z, t, title=r"$\omega_{{1}} = {:.3f}$".format(f1))


# Integrate the initial condition.
result = gnlse(
    z, t, u0,
    filtered_beta, filtered_gamma,
    kerr_op, dt=10, gpc=gpc)

if not result.successful:
    logging.error(
        "integration failed with error code %d",
        result.error_code)
    logging.error(
        "intermediate integration result is written in the output file.")


# Save the result.
numpy.savez(
    args.output,
    z=z, t=t, u=result.u, v=result.v,          # grids and fields
    a1=a1, a2=a2, t1=t1, t2=t2, f1=f1, f2=f2)  # soliton parameters
