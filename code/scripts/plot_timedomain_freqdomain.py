#!/usr/bin/env python3


from argparse import ArgumentParser

from matplotlib import pyplot as plot
from matplotlib import colors

import numpy
from numpy import fft

from common.plotter import pr_setup, pr_publish


parser = ArgumentParser()
parser.add_argument(
    "--vmin",
    help="minimun value for intensity plots",
    type=float,
    default=1E-6)
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

z = z / 1E4
t = t / 1E3
u = abs(u) / abs(u).max()
v = abs(v) / abs(v).max()


pr_setup()
plot.figure(figsize=(6, 8))

_, nt = u.shape
sst = int(nt / 1024)

_, nf = v[:, (f >= 0) & (f <= 4)].shape
ssf = int(nf / 1024)

u = u[:, ::sst]
v = v[:, ::ssf]

plot.subplot(2, 1, 1)
plot.pcolormesh(
    z, t[::sst], u.T**2,
    cmap="jet",
    norm=colors.LogNorm(vmin=args.vmin),
    rasterized=True,
    shading="auto")
plot.tight_layout(pad=4.0)
plot.xlabel(r"Distance $z$,~cm")
plot.ylabel(r"Delay $t$,~ps")

plot.subplot(2, 1, 2)
plot.pcolormesh(
    z, f[::ssf], v.T**2,
    cmap="jet",
    norm=colors.LogNorm(vmin=args.vmin),
    rasterized=True,
    shading="auto")
plot.tight_layout(pad=4.0)
plot.xlabel(r"Distance $z$,~cm")
plot.ylabel(r"Frequency $\omega$,~rad/fs")
plot.ylim(0, 4)

pr_publish(args.output)
