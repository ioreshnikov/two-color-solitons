#!/usr/bin/env python3


from argparse import ArgumentParser

from matplotlib import pyplot as plot
from matplotlib import colors
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
u = npz["u"]
v = npz["v"]

f = fft.fftfreq(len(t), t[1] - t[0])
f = 2 * numpy.pi * fft.fftshift(f)


# We do some aggressive windowing and sub-sampling to reduce the image
# size. Realistically speaking, we don't need images with the
# resolution larger than 1000×1000 pixels.

tw = (t > -500) & (t < +1000)
fw = (f > 0.5) & (f < 4)

t = t[tw]
f = f[fw]
u = u[:, tw]
v = v[:, fw]

nx = len(t)
ny = len(z)
ssx = int(nx / 1000)
ssy = int(ny / 1000)

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


# Finally, we plot everything
pr_setup()

plot.figure(figsize=(4, 3))

# Time domain subplot.
plot.subplot(2, 1, 1)
plot.pcolormesh(
    z, t, abs(u.T)**2,
    cmap="magma",
    norm=colors.LogNorm(vmin=1E-6),
    rasterized=True,
    shading="auto")
plot.xlim(0, z.max())
plot.ylim(-0.5, +0.5)
plot.yticks([-0.5, 0.0, +0.5])
plot.xlabel(r"Distance $z$, cm")
plot.ylabel(r"Delay $t$, ps")

plot.colorbar()

# Frequency domain subplot.
plot.subplot(2, 1, 2)
plot.pcolormesh(
    z, f, abs(v.T)**2,
    cmap="jet",
    norm=colors.LogNorm(vmin=1E-6),
    rasterized=True,
    shading="auto")
plot.xlim(0, z.max())
plot.ylim(0.5, 4)
plot.yticks(list(range(1, 5)))
plot.xlabel(r"Distance $z$, cm")
plot.ylabel(r"Freq. $\omega$, rad/fs")

plot.colorbar()

plot.tight_layout(pad=1.0)
pr_publish(args.output)
