#!/usr/bin/env python3


from argparse import ArgumentParser

import numpy
from matplotlib import pyplot as plot

from common.fiber import beta, beta1, gamma
from common.helpers import fundamental_soliton_dispersive_relation
from common.plotter import pr_setup, pr_publish


parser = ArgumentParser()
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

a1 = npz["a1"]; a1s = npz["a1s"]
a2 = npz["a2"]; a2s = npz["a2s"]

t1 = npz["t1"]; t1s = npz["t1s"]
t2 = npz["t2"]; t2s = npz["t2s"]

f1 = npz["f1"]; f1s = npz["f1s"]
f2 = npz["f2"]; f2s = npz["f2s"]


# Some parameters of the soliton.

# Group velocity of individual frequency components.
gv1 = beta1(f1s)
gv2 = beta1(f2s)

# Wavenumber of individual solitons interacting in a non-coherent
# manner.
f0 = 1.0
frame = beta(f1s) + beta1(f1s) * (f0 - f1s)
k1 = fundamental_soliton_dispersive_relation(f1s, t1s, f0) - frame + gamma(f1s) * a2s**2
k2 = fundamental_soliton_dispersive_relation(f2s, t2s, f0) - frame + gamma(f2s) * a1s**2


z /= 1E4  # switch to centimeters

pr_setup()

plot.figure(figsize=(3.2, 1.4 * 5))

plot.subplot(5, 1, 1)
plot.plot(z, a1s, color="black", label=r"$A_{1}$")
plot.plot(z, a2s, color="grey",  label=r"$A_{2}$")
plot.xlabel(r"Distance $z$, cm")
plot.ylabel(r"Amplitude, a.u")
plot.xlim(z.min(), z.max())
plot.legend(ncol=2, loc="upper center")

plot.subplot(5, 1, 2)
plot.plot(z, t1s, color="black", label=r"$T_{1}$")
plot.plot(z, t2s, color="grey",  label=r"$T_{2}$")
plot.ylabel(r"Duration $t$, fs")
plot.xlabel(r"Distance $z$, cm")
plot.xlim(z.min(), z.max())
plot.legend(ncol=2, loc="upper center")

plot.subplot(5, 1, 3)
plot.plot(z, f1s - f1, color="black", label=r"$\Delta \omega_{1}$")
plot.plot(z, f2s - f2, color="grey",  label=r"$\Delta \omega_{2}$")
plot.xlabel(r"Distance $z$, cm")
plot.ylabel(r"Fr. offset $\Delta \omega$, rad/fs")
plot.xlim(z.min(), z.max())
plot.legend(ncol=2, loc="upper center")

plot.subplot(5, 1, 4)
plot.plot(z, gv1, color="black", label=r"$GV_{1}$")
plot.plot(z, gv2, color="grey",  label=r"$GV_{2}$")
plot.xlabel(r"Distance $z$, cm")
plot.ylabel(r"Group velocity")
plot.xlim(z.min(), z.max())
plot.legend(ncol=2, loc="upper center")

plot.subplot(5, 1, 5)
plot.plot(z, k1, color="black", label=r"$k_{1}$")
plot.plot(z, k2, color="grey",  label=r"$k_{2}$")
plot.xlabel(r"Distance $z$, cm")
plot.ylabel(r"Wavenumber $k$, rad/$\mu$m")
plot.xlim(z.min(), z.max())
plot.legend(ncol=2, loc="upper center")

pr_publish(args.output)
