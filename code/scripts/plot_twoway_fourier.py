#!/usr/bin/env python3


from argparse import ArgumentParser

from matplotlib import colors
from matplotlib import pyplot as plot
from numpy import fft
import numpy

from common.fiber import beta, beta1, fundamental_soliton_dispersive_relation
from common.plotter import pr_setup, pr_publish


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

f1 = npz["f1"]
f2 = npz["f2"]
t1 = npz["t1"]
t2 = npz["t2"]

f = fft.fftfreq(len(t), t[1] - t[0])
f = 2 * numpy.pi * fft.fftshift(f)

k = fft.fftfreq(len(z), z[1] - z[0])
k = 2 * numpy.pi * fft.fftshift(k)


# We take a Fourier transform over z coordinate.
w = fft.fft(v, axis=0)
w = fft.fftshift(w, axes=0)


# We do some aggressive windowing and sub-sampling to reduce the image
# size. Realistically speaking, we don't need images with the
# resolution larger than 1000×1000 pixels.
tw = (t > -500) & (t < 1000)
fw = (f > 0) & (f < 5)

t = t[tw]
f = f[fw]
u = u[:, tw]
v = v[:, fw]
w = w[:, fw]

nx = len(f)
ssx = int(nx / 1000)

f = f[::ssx]
v = v[:, ::ssx]
w = w[:, ::ssx]

nx = len(t)
ssx = int(nx / 1000)

t = t[::ssx]
u = u[:, ::ssx]

# Next step is to normalize everything to the maximum value.
u = u / abs(u).max()
v = v / abs(v).max()
w = w / abs(w).max()


z = z / 1000


plot.figure()

plot.subplot(1, 3, 1)
plot.pcolormesh(
    t, z, abs(u)**2,
    cmap="magma",
    norm=colors.LogNorm(vmin=1E-4),
    rasterized=True,
    shading="auto")
plot.xlabel(r"$t$, fs")
plot.ylabel(r"$z$, mm")

plot.subplot(1, 3, 2)
plot.pcolormesh(
    f, z, abs(v)**2,
    cmap="jet",
    norm=colors.LogNorm(vmin=1E-4),
    rasterized=True,
    shading="auto")
plot.xlabel(r"$\omega$, rad/fs")
plot.ylabel(r"$z$, mm")

# The last panel we plot in a wavenumber and group-velocity
# compensated frame of reference.
frame = beta(f1) + beta1(f1) * (f - f1)

plot.subplot(1, 3, 3)
plot.pcolormesh(
    f, k, abs(w)**2,
    cmap="jet",
    norm=colors.LogNorm(vmin=1E-6),
    rasterized=True,
    shading="auto")
plot.plot(
    f, fundamental_soliton_dispersive_relation(f1, t1, f) - frame,
    color="magenta",
    linewidth=0.5)
plot.plot(
    f, fundamental_soliton_dispersive_relation(f2, t2, f) - frame,
    color="cyan",
    linewidth=0.5)
plot.plot(
    f, beta(f) - frame,
    color="white",
    linewidth=0.5)
plot.ylim(k.min(), k.max())
plot.xlabel(r"$\omega$, rad/fs")
plot.ylabel(r"$k$, rad/$\mu$m")

plot.show()
