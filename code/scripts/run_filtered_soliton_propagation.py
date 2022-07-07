#!/usr/bin/env python3


__doc__ = """
This script takes an output of `run_seed_soliton_propagation`, filters
the radiation tails from the soliton, shifts the center back to zero
and then integrates it further for another 10 cm with the same step of
10 μm.

This script is here mostly so that we can check if the stationary
state produced by the seed soliton is stationary enough. It is also
very useful to be used as a template for scattering simulations, since
one can reuse the filtering and shifting code as is.
"""


import argparse
import logging

import numpy

from common.fiber import beta, beta1, gamma, kerr_op
from common.helpers import frame_of_reference
from common.solver import gnlse


logging.basicConfig(
    level=logging.DEBUG,
    format="%(process)s %(asctime)s %(levelname)s %(message)s")


parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "input",
    help="path to the input .npz file",
    type=str)
parser.add_argument(
    "output",
    help="path to the output .npz file",
    type=str,
    default=1.000)
args = parser.parse_args()


npz = numpy.load(args.input)
t = npz["t"]
u = npz["u"]


# Those parameters are not used in the simulation, but it's nice to
# have them in the output file.
a1 = npz["a1"]
a2 = npz["a2"]
t1 = npz["t1"]
t2 = npz["t2"]
f1 = npz["f1"]
f2 = npz["f2"]


# We run the simulation for 10 cm more. Should be enough to really see
# the structure of the radiation in a clean setting.
z = numpy.linspace(0, 1E5, int(1E3))


# We construct a super-Gaussian window around the soliton peak.
# Multiplying the soliton by that window we suppress the radiation
# tails.
u0 = u[-1, :]
idx_max = abs(u0).argmax()
t0 = t[idx_max]
w0 = 100

window = numpy.exp(- ((t - t0)/w0)**8)
u0 = window * u0


# We also perform a shift of the solution so that the center of the
# soliton is at the origin t=0.
idx_zero = abs(t).argmin()
u0 = numpy.roll(u0, idx_zero - idx_max)


# Construct a frequency filtered and group velocity compensated
# dispersive profile.
def filtered_beta(f):
    b = beta(f) - frame_of_reference(f, f1)
    b[(f <= 0.5) | (f >= 0.75 * f.max())] = 0
    return b


# Frequency filter the nonlinear coefficient.
def filtered_gamma(f):
    g = gamma(f)
    g[(f <= 0.5) | (f >= 0.75 * f.max())] = 0
    return g


# Integrate the initial condition.
result = gnlse(
    z, t, u0,
    filtered_beta, filtered_gamma,
    kerr_op, dt=10)

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
