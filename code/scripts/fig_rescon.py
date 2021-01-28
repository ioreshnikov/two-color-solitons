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
    estimate_soliton_parameters,
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
u0 = u[-1, :]

del u  # delete the time-domain matrix to free some memory


# From the full output spectrum we estimate the final soliton
# parameters using the original values as an estimate
a1 = npz["a1"]; a2 = npz["a2"]
t1 = npz["t1"]; t2 = npz["t2"]
f1 = npz["f1"]; f2 = npz["f2"]
a1, a2, t1, t2, f1, f2 = estimate_soliton_parameters(
    t, u0, a1, a2, t1, t2, f1, f2)


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

k1 = fundamental_soliton_dispersive_relation(f1, t1, f) - frame + gamma(f1) * a2**2
k2 = fundamental_soliton_dispersive_relation(f2, t2, f) - frame + gamma(f2) * a1**2
b0 = beta(f) - frame

kmin = min(
    k1.min(), k2.min(),
    (2*k1 - k2).min(), (2*k2 - k1).min(),
    b0[(f > f1) & (f < f2)].min())
kmax = max(k1.max(), k2.max(), b0.max())
margin = 0.25 * (kmax - kmin)

kmin = kmin - margin
kmax = kmax + margin


# The narrowing itself
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

plot.plot(f, b0, color="white")
plot.plot(f, k1, color="white", linestyle="solid", label="$k_{1}$")
plot.plot(f, k2, color="white", linestyle="dashed", label="$k_{2}$")

nlabels = 2
if any((2*k2 - k1) < b0):
    nlabels += 1
    plot.plot(
        f, 2*k2 - k1,
        color="white",
        linestyle="dotted",
        label="$2k_{2} - k_{1}$")
if any((2*k1 - k2) < b0):
    nlabels += 1
    plot.plot(
        f, 2*k1 - k2,
        color="yellow",
        linestyle="dotted",
        label="$2k_{1} - k_{2}$")


plot.ylim(k.min(), k.max())
plot.xlabel(r"Frequency $\omega$, rad/fs")
plot.ylabel(r"Wavenbr. $k$, rad/$\mu$m")
plot.gca().xaxis.set_major_locator(MultipleLocator(0.5))
legend = plot.legend(ncol=nlabels, loc="upper center")
for text in legend.get_texts():
    text.set_color("white")

plot.tight_layout(pad=1.0)
pr_publish(args.output)
