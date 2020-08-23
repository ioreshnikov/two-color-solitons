#!/usr/bin/env python3


__doc__ = """
This script continues GNLSE integration on a very fine grid
with an initial condition that is read from another simulation
output. The main goal is to construct a numeric solution on a
grid that allows really fine spatial spectral resolution.

The computational grid is defined so that the total wavenumber
span is equal π and resolution is equal to 0.001 rad/μm. This
should be enough to produce a publication-ready 2-way Fourier
image of the solution.
"""


import argparse
import logging

import numpy

from common.fiber import beta, beta1, gamma, kerr_op
from common.solver import gnlse


logging.basicConfig(level=logging.INFO)


parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "input",
    help="path to the input .npz file",
    type=str,
    default=1.000)
parser.add_argument(
    "output",
    help="path to the output .npz file",
    type=str,
    default=1.000)
parser.add_argument(
    "-z0",
    help="use this z coordinate as the starting point",
    type=float)
args = parser.parse_args()


npz = numpy.load(args.input)
z = npz["z"]
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


# Some sanity checks.
if args.z0 is not None and not z[0] <= args.z0 <= z[-1]:
    logging.fatal("-z0 is outside of [{}, {}]".format(z[0], z[-1]))
    exit()


# Find the index that corresponds to the initial condition.
if not args.z0:
    idx = -1
else:
    idx = abs(z - args.z0).argmin()
z0 = z[idx]


# Read the initial condition. It is assumed to be analytic at this
# point.
u0 = u[idx, :]


# Pick a new z grid so that it has a desired spectral resolution of
# 0.001 rad/μm and a total spectral span of π.
dk = 0.001
z = numpy.arange(z0, z0 + int(2*numpy.pi/dk), 2)


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


# Integrate the initial condition.
result = gnlse(
    z, t, u0,
    filtered_beta, filtered_gamma, kerr_op)

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
