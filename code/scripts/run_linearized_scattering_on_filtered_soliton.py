#!/usr/bin/env python3


__doc__ = """
This script
"""


import argparse
import logging

import numpy

from common.fiber import beta, beta1, beta2, gamma
from common.solver import linsfeq
from common.helpers import (
    estimate_soliton_parameters, frame_of_reference,
    fundamental_soliton_dispersive_relation, freqs, sech)
from common.plotter import gpc_setup


logging.basicConfig(
    level=logging.DEBUG,
    format="%(process)s %(asctime)s %(levelname)s %(message)s")


parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "-fi",
    help="frequency of the incident dispersive wave",
    type=float,
    default=None)
parser.add_argument(
    "-t", "--terms",
    nargs="*",
    help="resonances terms to include",
    type=str,
    choices=(
        "pot.1", "pot.2", "pot",
          "7.1",   "7.2",   "7",
          "8.1",   "8.2",   "8",
         "12.1",  "12.2",  "12",
         "13.1",  "13.2",  "13",
         "14.1",  "14.2",  "14",
           "15"),
    default=())
parser.add_argument(
    "--title",
    help="title to be used for interactive plotter",
    type=str)
parser.add_argument(
    "input",
    help="path to the input .npz file",
    type=str)
parser.add_argument(
    "output",
    help="path to the output .npz file",
    type=str)
# parser.add_argument(
#     "canon",
#     help="path to the canonical output .npz file",
#     type=str)
args = parser.parse_args()


npz = numpy.load(args.input)
t = npz["t"]
u = npz["u"]


# We start by estimating the soliton parameters. The values saved in the seed
# npz are used as an initial approximation.
a1 = npz["a1"]
a2 = npz["a2"]
t1 = npz["t1"]
t2 = npz["t2"]
f1 = npz["f1"]
f2 = npz["f2"]

u0 = u[-1, :]
a1, a2, t1, t2, f1, f2 = estimate_soliton_parameters(
    t, u0, a1, a2, t1, t2, f1, f2)


# We downsize the grid. The original one comes from the seed simulations where
# it was necessary so that the Cherenkov radiation does now travel around.
tw = (t > -3000) & (t <= 3000)
t = t[tw]
f = freqs(t)
fmin = f.min()
fmax = f.max()

u10 = a1 * sech(t / t1) * numpy.exp(-1j * f1 * t)
u20 = a2 * sech(t / t2) * numpy.exp(-1j * f2 * t)

k1 = (
    fundamental_soliton_dispersive_relation(f1, t1, f1)
    - frame_of_reference(f1, f1)
    + gamma(f1) * a2**2)
k2 = (
    fundamental_soliton_dispersive_relation(f2, t2, f2)
    - frame_of_reference(f1, f2)
    + gamma(f2) * a1**2)


if args.fi:
    # Let us proceed with the dispersive wave first we pick the amplitude.
    ai = 0.01 * abs(u0).max()
    fi = args.fi
    ti =  300
    t0 = 1000
    t0 = 100
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

    # Find the coordinate that should correspond to the middle of the
    # collision process.
    zc = abs(t0) / abs(vg)
    logging.info("z = {:.3f} cm should be the middle of collision".format(zc / 10000))

    # We create the distance grid symmetric to the collision point.
    zmax = min(2 * zc, 2E5)
    logging.info("picking z = {:.3f} cm as simulation endpoint".format(zmax / 10000))

    z = numpy.linspace(0, zmax, int(1E3))

    # Finally we construct the incident radiation
    pi0 = ai * numpy.exp(-(t - t0)**2 / ti**2) * numpy.exp(- 1j * fi * t)
else:
    z = npz["z"]
    pi0 = numpy.zeros_like(t)


# Construct a frequency filtered and group velocity compensated
# dispersive profile

def filtered_beta(f):
    b = beta(f) - frame_of_reference(f, f1)
    b[(f <= 0.5) | (f >= 0.75 * fmax)] = 0
    return b


def filtered_beta1(f):
    b = (
        beta(f1) +
        beta1(f1) * (f - f1) +
        0.5 * beta2(f1) * (f - f1)**2 -
        frame_of_reference(f, f1))
    b[(f <= 0.5) | (f >= 0.75 * fmax)] = 0
    return b


def filtered_beta2(f):
    b = (
        beta(f2) +
        beta1(f2) * (f - f2) +
        0.5 * beta2(f2) * (f - f2)**2 -
        frame_of_reference(f, f1))
    b[(f <= 0.5) | (f >= 0.75 * fmax)] = 0
    return b


def filtered_gamma(f):
    g = gamma(f)
    g[(f <= 0.5) | (f >= 0.75 * fmax)] = 0
    return g


# Setup interactive plotting
if args.title:
    title = args.title
else:
    if args.fi:
        title = r"$\omega_{{1}} = {:.3f}$, $\omega_{{i}} = {:.3f}$".format(f1, fi)
    else:
        title = r"$\omega_{{1}} = {:.3f}$".format(f1)

gpc = gpc_setup(z, t, title=title, fmin=0.5, fmax=4, tmin=-1, tmax=+1)

# Integrate the initial condition.
result = linsfeq(
    z, t, u10, u20, k1, k2, pi0,
    filtered_beta,
    filtered_beta1,
    filtered_beta2,
    filtered_gamma,
    terms=args.terms,
    dt=10, gpc=gpc)

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
