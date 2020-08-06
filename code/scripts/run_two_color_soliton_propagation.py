#!/usr/bin/env python3


from functools import partial
import argparse
import logging

import numpy

from common.fiber import (
    beta, beta1, beta2,
    gamma, kerr_op, linear_absorption_op,
    gv_matching_frequencies,
    fundamental_soliton_amplitude)
from common.helpers import sech, to_analytic
from common.solver import gnlse_propagate


parser = argparse.ArgumentParser()
parser.add_argument(
    "output",
    help="path to the output .npz file",
    type=str)
args = parser.parse_args()


logging.basicConfig(level=logging.INFO)


# Defined the computational grid.
z = numpy.linspace(0, 1E5, int(1E3))
t = numpy.linspace(-2000, +2000, 2**13)

# Pick the carrier frequency of the first pulse at then find a
# frequency for the second pulse by picking a frequency in a region of
# anomalous dispersion so that the group velocity of both pulses match.
f1 = 3.0
f2 = gv_matching_frequencies(f1)[0]

assert beta2(f1) < 0
assert beta2(f2) < 0


# Define the initial condition as a sum of two fundamental solitons of
# equal width at the selected frequencies.
t1 = 25
a1 = fundamental_soliton_amplitude(f1, t1)
u1 = a1 * sech(t / t1) * numpy.exp(-1j * f1 * t)

t2 = 25
a2 = fundamental_soliton_amplitude(f2, t2)
u2 = a2 * sech(t / t2) * numpy.exp(-1j * f2 * t)

u0 = to_analytic((u1 + u2).real)


# Construct a frequency filtered and group velocity compensated
# dispersive profile.
def filtered_beta(f):
    b = beta(f) - beta1(f1) * f
    b[(f <= 0) | (f >= 0.75 * f.max())] = 0
    return b


# Frequency filter the nonlinear coefficient.
def filtered_gamma(f):
    g = gamma(f)
    g[(f <= 0) | (f >= 0.75 * f.max())] = 0
    return g


# Add radiation absorber close to the edges of the computational
# domain.
absorber = (
    0.001 * sech((t - t.min()) / 20) +
    0.001 * sech((t - t.max()) / 20))


# Integrate the initial condition.
result = gnlse_propagate(
    z, t, u0,
    filtered_beta, filtered_gamma,
    kerr_op, partial(linear_absorption_op, profile=absorber), dt=10)

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
