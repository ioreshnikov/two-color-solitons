#!/usr/bin/env python3


from argparse import ArgumentParser

from matplotlib import pyplot as plot
from matplotlib import colors
from matplotlib.ticker import MultipleLocator
from numpy import fft
import numpy

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
t = npz["t"]

f1 = npz["f1"]
f2 = npz["f2"]
t1 = npz["t1"]
t2 = npz["t2"]

f = fft.fftfreq(len(t), t[1] - t[0])
f = 2 * numpy.pi * fft.fftshift(f)


# We don't need the full time-domain plot anymore and we can filter
# and downsample it.
u = npz["u"]

tw = (t > -500) & (t < +500)
t = t[tw]
u = u[:, tw]

ssx = int(len(t) / 1000) or 1
t = t[::ssx]
u = u[:, ::ssx]

u = abs(u) / abs(u).max()


# From the spectrum we need to take the input and the output slices.
v = npz["v"]
v0 = v[ 0, :]
v1 = v[-1, :]

del v

vn = abs(v0).max()
v0 = abs(v0) / vn
v1 = abs(v1) / vn


# To makes the plots more readable, we rescale the time to picoseconds
# and distance to centimeters.
t = t / 1000
z = z / 10000


pr_setup()

plot.figure(figsize=(3.4, 2.8))

plot.subplot(2, 1, 1)
plot.pcolormesh(
    z, t, u.T**2,
    cmap="magma",
    norm=colors.LogNorm(vmin=1E-6),
    rasterized=True,
    shading="auto")
plot.xlim(z.min(), z.max())
plot.ylim(-0.5, +0.5)
plot.yticks([-0.5, 0.0, +0.5])
plot.xlabel(r"Distance $z$, cm")
plot.ylabel(r"Delay $t$, ps")
plot.colorbar()

plot.annotate(
    "(a)",
    xy=(1.0, 0.0),
    xytext=(-18.0, +10.0),
    xycoords="axes fraction",
    textcoords="offset points",
    color="white")

ax = plot.subplot(2, 1, 2)
plot.semilogy(f, v0, color="gray", label="in", zorder=10)
plot.semilogy(f, v1, color="black", label="out")
plot.legend(ncol=2, loc="upper center")
plot.xlim(0.5, 4.0)
plot.ylim(1E-4, 5.0)
plot.xlabel(r"Frequency $\omega$, rad/fs")
plot.ylabel(r"$\left| \tilde u(\omega) \right|$, a.\,u.")
ax.xaxis.set_major_locator(MultipleLocator(0.5))

plot.annotate(
    "(b)",
    xy=(1.0, 1.0),
    xytext=(-18.0, -12.0),
    xycoords="axes fraction",
    textcoords="offset points")

plot.tight_layout(pad=1.0)
pr_publish(args.output)
