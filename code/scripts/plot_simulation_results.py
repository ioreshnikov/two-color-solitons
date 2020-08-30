#!/usr/bin/env python3


from argparse import ArgumentParser

from matplotlib import colors
from matplotlib import pyplot as plot
from numpy import fft
import numpy


parser = ArgumentParser()
parser.add_argument(
    "input",
    help="path to an .npz file with the simulation results",
    type=str)
args = parser.parse_args()


npz = numpy.load(args.input)
z = npz["z"]
t = npz["t"]
u = npz["u"]
v = npz["v"]

f = fft.fftfreq(len(t), t[1] - t[0])
f = 2 * numpy.pi * fft.fftshift(f)


# We do some aggressive windowing and sub-sampling to reduce the image
# size. Realistically speaking, we don't need images with the
# resolution larger than 1000×1000 pixels.

tw = (t > -500) & (t < +1000)
fw = (f > 0) & (f < 5)

t = t[tw]
f = f[fw]
u = u[:, tw]
v = v[:, fw]

nx = len(t)
ny = len(z)
ssx = int(nx / 1000) or 1
ssy = int(ny / 1000) or 1

t = t[::ssx]
f = f[::ssx]
z = z[::ssy]

u = u[::ssy, ::ssx]
v = v[::ssy, ::ssx]

# Next step is to normalize everything to the maximum value.
u = u / abs(u).max()
v = v / abs(v).max()


# To makes the plots more readable, we rescale the time to picoseconds
# and distance to centimeters.
t = t / 1000
z = z / 10000


plot.figure()

plot.subplot(1, 2, 1)
plot.pcolormesh(
    t, z, abs(u)**2,
    cmap="magma",
    norm=colors.LogNorm(vmin=1E-6),
    shading="auto")
plot.xlabel(r"$t$, ps")
plot.ylabel(r"$z$, cm")

plot.subplot(1, 2, 2)
plot.pcolormesh(
    f, z, abs(v)**2,
    cmap="jet",
    norm=colors.LogNorm(vmin=1E-6),
    shading="auto")
plot.xlabel(r"$\omega$, rad/fs")
plot.ylabel(r"$z$, cm")

plot.show()
