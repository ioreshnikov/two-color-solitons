#!/usr/bin/env python3


__doc__ = """
This script extracts parameters of individual soliton components from the
integration result produced by `run_...` scripts. The parameters are
estimated by approximating the simulated spectrum as a sum of two
hyperbolic secants. For performance reason, in this script we do not
require the group velocities of the individual components to be equal.
"""


import argparse

from matplotlib import pyplot as plot
import numpy

from common.helpers import (
    _envelope_ansatz, freqs, estimate_soliton_parameters)


parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "-i", "--interactive",
    help="Plot estimated spectrum as we go",
    action="store_true")
parser.add_argument(
    "input",
    help="path to an .npz file with the simulation results",
    type=str)
parser.add_argument(
    "output",
    help="path to the output image",
    type=str)
args = parser.parse_args()


npz = numpy.load(args.input)
z = npz["z"]
t = npz["t"]
u = npz["u"]
v = npz["v"]

f = freqs(t, shift=True)


# From the full output spectrum we estimate the final soliton
# parameters using the original values as an initial guess
a1 = npz["a1"]; a2 = npz["a2"]
t1 = npz["t1"]; t2 = npz["t2"]
f1 = npz["f1"]; f2 = npz["f2"]

a1s = []; a2s = []
t1s = []; t2s = []
f1s = []; f2s = []


if args.interactive:
    plot.ion()

    plot.figure(1)
    up, *_ = plot.plot(f, 0 * f, label=r"$\tilde u$")
    ep, *_ = plot.plot(f, 0 * f, label=r"est.")

    n0 = abs(v[0, :]).max()
    plot.xlim(0.5, 4.0)
    plot.ylim(0.0, 1.25 * n0)

    plot.xlabel(r"Frequency $\omega$, rad/fs")
    plot.legend()

    plot.show(block=False)


for idx in range(len(z)):
    print("{}/{}".format(idx, len(z)))

    u0 = u[idx, :]
    a1, a2, t1, t2, f1, f2 = estimate_soliton_parameters(
        t, u0, a1, a2, t1, t2, f1, f2, match_f=False)

    a1s.append(a1); a2s.append(a2)
    t1s.append(t1); t2s.append(t2)
    f1s.append(f1); f2s.append(f2)

    if args.interactive:
        v0 = v[idx, :]
        up.set_data(f, abs(v0))
        ep.set_data(f, _envelope_ansatz(t, f, a1, a2, t1, t2, f1, f2))
        plot.pause(1E-3)


a1s = numpy.array(a1s); a2s = numpy.array(a2s)
t1s = numpy.array(t1s); t2s = numpy.array(t2s)
f1s = numpy.array(f1s); f2s = numpy.array(f2s)


numpy.savez(
    args.output,
    z=z,
    a1=npz["a1"], a2=npz["a2"],
    t1=npz["t1"], t2=npz["t2"],
    f1=npz["f1"], f2=npz["f2"],
    a1s=a1s, a2s=a2s,
    t1s=t1s, t2s=t2s,
    f1s=f1s, f2s=f2s)
