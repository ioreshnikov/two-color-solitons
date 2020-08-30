#!/usr/bin/env python3


from argparse import ArgumentParser

from matplotlib import colors
from matplotlib import pyplot as plot
from matplotlib.ticker import MultipleLocator
from numpy import fft
from scipy.optimize import curve_fit
import numpy

from common.fiber import (
    beta, beta1, gamma,
    fundamental_soliton_dispersive_relation,
    fundamental_soliton_width)
from common.helpers import filter_tails, sech, zeros_on_a_grid
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
t = npz["t"]
z = npz["z"]

f0 = npz["f1"]
fi = npz.get("fi")

# ^^^ the frequency above is the reference frequency which coincides
# with the initial frequency of the first soliton. However, during the
# soliton evolution it sheds Cherenkov radiation and accumulates
# significant frequency shift in both of the components. We have to
# deduce the current soliton parameters from the simulation data.

# It is a bit involved process. To determine the frequencies we cannot
# simply use the raw spectrum, due to the spectral interference with
# Cherenkov radiation. However, the radiation at this point is well
# separated from the soliton, therefore we can filter it in the time
# domain (see `filter_tails` utility). Then we can transform back into
# frequency domain, separate the solitons in there and then find the
# spectral maxima. Then we transform each of the separated solitons
# back in time domain to determine each of the soliton's amplitude,
# which in turn is used to estimate the soliton width.

f = fft.fftfreq(len(t), t[1] - t[0])
f = 2 * numpy.pi * fft.fftshift(f)

u = npz["u"]
u0 = u[0, :]

del u  # delete the time-domain matrix to free some memory

u0 = filter_tails(t, u0)
v0 = fft.fftshift(fft.ifft(u0))

# We take a weighted average in left and right hand's sides of the
# spectrum using square of the spectrum as a weight.
fl = f[f < 2]
fr = f[f > 2]
vl = (abs(v0)**2)[f < 2]
vr = (abs(v0)**2)[f > 2]

# The distributions are a bit lopsided, so we suppress the tails and
# weight using only the dome.
vl[vl < 0.5 * vl.max()] = 0
vr[vr < 0.5 * vr.max()] = 0

f1 = sum(fl * vl) / sum(vl)
f2 = sum(fr * vr) / sum(vr)

u1 = fft.fftshift(fft.fft(v0[f < 2]))
u2 = fft.fftshift(fft.fft(v0[f > 2]))

a1 = abs(u1).max()
a2 = abs(u2).max()

t1 = fundamental_soliton_width(a1, f1)
t2 = fundamental_soliton_width(a2, f2)


# Now that we're done with the spectrum of the cleaned signal we can
# proceed to the real messy data :)
v = npz["v"]


# We do some aggressive windowing and sub-sampling to reduce the image
# size. Realistically speaking, we don't need images with the
# resolution larger than 1000×1000 pixels.

# We do it in two steps. First we narrow the frequency domain.
fw = (f > 0.5) & (f < 4.0)

f = f[fw]
v = v[:, fw]

nf = len(f)
ssx = int(nf / 1000) or 1

f = f[::ssx]
v = v[:, ::ssx]


# `v` now a much smaller matrix and we should be able to take a
# Fourier transform over it on your average machine.
w = fft.fft(v, axis=0)
w = fft.fftshift(w, axes=0)

del v  # delete the frequency-domain matrix to free some memory

k = fft.fftfreq(len(z), z[1] - z[0])
k = 2 * numpy.pi * fft.fftshift(k)


# Now we can narrow in wavenumber domain. We calculate the bounds so
# that the dispersive curve fits into the picture somewhat nicely.
frame = beta(f0) + beta1(f0) * (f - f0)

b = beta(f) - frame
k1 = fundamental_soliton_dispersive_relation(f1, t1, f) - frame + gamma(f1) * a2**2
k2 = fundamental_soliton_dispersive_relation(f2, t2, f) - frame + gamma(f2) * a1**2

kmin = min(b.min(), k1.min(), k2.min())
kmax = max(b.max(), k1.max(), k2.max())
margin = 0.25 * (kmax - kmin)
kmax = kmax + margin

kw = (k > kmin) & (k < kmax)

k = k[kw]
w = w[kw, :]

nk = len(k)
ssy = int(nk / 1000) or 1

w = w[::ssy, :]
w = w / abs(w).max()


pr_setup()

plot.figure(figsize=(3.4, 2.8))

plot.pcolormesh(
    f, k, abs(w)**2,
    cmap="jet",
    norm=colors.LogNorm(vmin=1E-8),
    rasterized=True,
    shading="auto")

plot.plot(
    f, b,
    color="white",
    linewidth=0.5)
plot.annotate(
    r"$\beta(\omega)$",
    xy=(f[0], b[0]),
    xytext=(5.0, 0.0),
    xycoords="data",
    textcoords="offset points",
    color="white")

if not fi:
    plot.plot(
        f, k1,
        color="white",
        linewidth=0.5)
    plot.annotate(
        r"$k_{1}(\omega)$",
        xy=(f[-1], k1[-1]),
        xytext=(-28.0, -12.0),
        xycoords="data",
        textcoords="offset points",
        color="white")

if not fi:
    plot.plot(
        f, k2,
        color="white",
        linewidth=0.5)
    plot.annotate(
        r"$k_{2}(\omega)$",
        xy=(f[-1], k2[-1]),
        xytext=(-28.0, -12.0),
        xycoords="data",
        textcoords="offset points",
        color="white")

if not fi:
    roots = zeros_on_a_grid(f, 2*k1 - k2 - b)
    if len(roots):
        plot.plot(
            f, 2*k1 - k2,
            color="white",
            linewidth=0.5,
            linestyle="dashed")

if not fi:
    roots = zeros_on_a_grid(f, 2*k2 - k1 - b)
    if len(roots):
        plot.plot(
            f, 2*k2 - k1,
            color="white",
            linewidth=0.5,
            linestyle="dashed")

bi = None
if fi:
    bi = beta(fi) - beta(f0) - beta1(f0) * (fi - f0)
    plot.plot(
        f, bi * numpy.ones_like(f),
        color="white",
        linewidth=0.5)

if fi:
    roots = zeros_on_a_grid(f, 2 * k1 - bi - b)
    if len(roots):
        plot.plot(
            f, 2*k1 - bi,
            color="white",
            linewidth=0.5,
            linestyle="dashed")

if fi:
    roots = zeros_on_a_grid(f, 2 * k2 - bi - b)
    if len(roots):
        plot.plot(
            f, 2*k2 - bi,
            color="white",
            linewidth=0.5,
            linestyle="dashed")

if fi:
    roots = zeros_on_a_grid(f, k1 + k2 - bi - b)
    if len(roots):
        plot.plot(
            f, k1 + k2 - bi,
            color="white",
            linewidth=0.5,
            linestyle="dashed")

plot.ylim(k.min(), k.max())
plot.xlabel(r"Frequency $\omega$, rad/fs")
plot.ylabel(r"Wavenbr. $k$, rad/$\mu$m")
plot.gca().xaxis.set_major_locator(MultipleLocator(0.5))

plot.annotate(
    "(c)",
    xy=(1.0, 1.0),
    xytext=(-18.0, -12.0),
    xycoords="axes fraction",
    textcoords="offset points",
    color="white")

plot.tight_layout(pad=1.0)
pr_publish(args.output)
