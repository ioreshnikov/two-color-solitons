#!/usr/bin/env python3


__doc__ = """
This script simulates scattering of a weak dispersive wave on a soliton.
The soliton is taken from an output of `run_seed_soliton_propagation`,
filtered and shifted back to zero. To that initial condition added
a weak dispersive wave in the shape of a weak Gaussian pulse. The result
is then propagated for some distance that is enough for scattering to
happen and for the reflected and the transmitted components to be fully
separated from the soliton.

The time grid is taken from the soliton's .npz file. The distance grid
is automatically computed so that at the output of the simulation the
transmitted component of the dispersive wave is put at a position symmetric
w.r.t the soliton.
"""


import argparse
import logging

import numpy

from common.fiber import beta, beta1, gamma, kerr_op
from common.solver import gnlse
from common.helpers import filter_tails, to_analytic


logging.basicConfig(level=logging.INFO)


parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "-fi",
    help="frequency of the incident dispersive wave",
    type=float,
    default=2.000)
parser.add_argument(
    "input",
    help="path to the input .npz file",
    type=str)
parser.add_argument(
    "output",
    help="path to the output .npz file",
    type=str)
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


# We suppress the radiation tails.
u0 = u[-1, :]
u0 = filter_tails(t, u0)


# We also perform a shift of the solution so that the center of the
# soliton is at the origin t=0.
idx_zero = abs(t).argmin()
idx_max = abs(u0).argmax()
u0 = numpy.roll(u0, idx_zero - idx_max)


# We're done with the soliton. Let us proceed with the dispersive wave.
# first we pick the amplitude.
ai = 0.01 * abs(u0).max()
fi = args.fi
ti =  300
t0 = 1000
logging.info("DW amplitude is {:.3f}".format(ai))


# Let's calculate the parameters of the scattering process.
vg = beta1(fi) - beta1(f1)

if t0 * vg > 0:
    logging.info(
        "DW initial delay is {:.3f} and group velocity is {:.3f}"
        .format(t0, vg))
    logging.info("swapping the sign of t0")
    t0 = -t0
    logging.info("new soliton position t0 = {}".format(t0))

dw = ai * numpy.exp(-(t - t0)**2 / ti**2) * numpy.exp(1j * fi * t)
u0 = to_analytic(u0.real + dw.real)


# Find the coordinate that should correspond to the middle of the
# collision process.
zc = abs(t0) / abs(vg)
logging.info("z = {:.3f} cm should be the middle of collision".format(zc / 10000))


# We create the distance grid symmetric to the collision point.
zmax = min(2 * zc, 2E5)
logging.info("picking z = {:.3f} cm as simulation endpoint".format(zmax / 10000))

z = numpy.linspace(0, zmax, int(1E3))


# Construct a frequency filtered and group velocity compensated
# dispersive profile.
def filtered_beta(f):
    b = beta(f) - beta(f1) - beta1(f1) * f
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
    filtered_beta, filtered_gamma, kerr_op, dt=10)

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
    a1=a1, a2=a2, t1=t1, t2=t2, f1=f1, f2=f2,  # soliton parameters
    ai=ai, fi=fi)                              # DW parameters
