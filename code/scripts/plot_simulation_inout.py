#!/usr/bin/env python3


from argparse import ArgumentParser

from matplotlib import pyplot as plot
from scipy import fft
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


plot.figure()

plot.subplot(1, 2, 1)
plot.plot(t, u[ 0, :].real / abs(u[0, :].real).max())
plot.plot(t, u[-1, :].real / abs(u[0, :].real).max())
plot.xlim(-500, +1000)
plot.xlabel("Delay, ps")

plot.subplot(1, 2, 2)
plot.plot(f, abs(v[ 0, :]) / abs(v[0, :]).max())
plot.plot(f, abs(v[-1, :]) / abs(v[0, :]).max())
plot.xlim(0, 4)
plot.xlabel("Frequency, rad/fs")

plot.show()
