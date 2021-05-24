#!/usr/bin/env python3


__doc__ = """
This script produces a five-panel plot of the simulation results for
the specific case of high-intensity scattering with ω₁=1.010 and ωᵢ = 1.100.
The first panel is the time domain view, the second and the third one are
in a single row with soliton amplitude and soliton frequency plots
side-by-side, the fourth one is the input and the output spectra and the
last one is the resonance condition plot.

It is no way generalizable and we hard-code all the extra parameters in here.
"""


import argparse

from matplotlib import pyplot as plot
from matplotlib import colors
from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator
from matplotlib.legend_handler import HandlerBase

from numpy import fft
import numpy

from common.fiber import beta, beta1, gamma
from common.helpers import (
    estimate_soliton_parameters,
    fundamental_soliton_dispersive_relation, zeros_on_a_grid)
from common.plotter import XSMALL_FONT_SIZE, pr_setup, pr_publish


parser = argparse.ArgumentParser()
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
    "parameters",
    help="path to an .npz file with extracted soliton parameters",
    type=str)
parser.add_argument(
    "output",
    help="path to the output image",
    type=str)
args = parser.parse_args()


npz = numpy.load(args.input)
z = npz["z"]
t = npz["t"]

f = fft.fftfreq(len(t), t[1] - t[0])
f = 2 * numpy.pi * fft.fftshift(f)

u = npz["u"]

# From the full output spectrum we estimate the final soliton
# parameters using the original values as an estimate
u0 = u[0, :]
a1 = npz["a1"]; a2 = npz["a2"]
t1 = npz["t1"]; t2 = npz["t2"]
f1 = npz["f1"]; f2 = npz["f2"]
a1, a2, t1, t2, f1, f2 = estimate_soliton_parameters(
    t, u0, a1, a2, t1, t2, f1, f2)

# Read incident frequency
fi = npz["fi"]


# Now that we have estimated everything, we don't need the full output
# at the original resolution, so we can start windowing and
# subsampling the matrices.
tw = (t > -1000) & (t < +1000)
t = t[tw]
u = u[:, tw]

ssx = int(len(t) / 1000) or 1
t = t[::ssx]
u = u[:, ::ssx]

u = abs(u) / abs(u).max()


# From the spectrum we need to take the input and the output slices.
# To make our life easier when calculating the resonance condition, we
# also narrow the frequency range in here.
fw = (f > 0.5) & (f < 4.0)

f = f[fw]
v = npz["v"]
v0 = v[ 0, fw]
v1 = v[-1, fw]

del v

vn = abs(v0).max()
v0 = abs(v0) / vn
v1 = abs(v1) / vn


# To makes the plots more readable, we rescale the time to picoseconds
# and distance to centimeters.
t = t / 1000
z = z / 10000

# Frame of reference for resonance condition plotting and dispersive
# relations for the solitons and the medium.
frame = beta(f1) + beta1(f1) * (f - f1)

k1 = fundamental_soliton_dispersive_relation(f1, t1, f) - frame + gamma(f1) * a2**2
k2 = fundamental_soliton_dispersive_relation(f2, t2, f) - frame + gamma(f2) * a1**2
b0 = beta(f) - frame

# In this specific scattering process we have an oscillating soliton. This means
# that when we see β(ωᵢ) in the resonance condition, the corresponding line will
# always be split in (generally speaking) infinite series of lines, each
# separated by a wavenumber 2π / z₀, where in this specific case z₀ = 2 mm.
z0 = 0.2 * 10000  # period in μm
k0 = 2 * numpy.pi / z0

harmonics = range(-3, 4)


rfs = []

framei = beta(f1) + beta1(f1) * (fi - f1)
bi = beta(fi) - framei

# 1. Resonance due to β(ωₛ) = β(ωᵢ) (but exclude ωᵢ)
srfi = []
for n in harmonics:
    srfi.append([
        rf for rf in zeros_on_a_grid(f, b0 - bi - n*k0)
        if abs(rf - fi) > 0.01
    ])

# 2. Resonances due to β(ω) = k₂(ω) - k₁(ω) + β(ωᵢ)
srf21 = []
for n in harmonics:
    srf21.append(list(zeros_on_a_grid(f, b0 - k2 + k1 - bi - n*k0)))

rfs.extend(srfi)
rfs.extend(srf21)


# Load soliton parameters
npz = numpy.load(args.parameters)
a1s = npz["a1s"]
a2s = npz["a2s"]
f1s = npz["f1s"]
f2s = npz["f2s"]


# Finally, we plot everything.
pr_setup()

height_ratios = [2/3, 3/4, 2/3, 1]
plot.figure(
    figsize=(3.2, sum(1.4 * h for h in height_ratios)))
gs = gridspec.GridSpec(4, 2, height_ratios=height_ratios)


# First panel (optional): time-domain plot
plot.subplot(gs[0, :])
plot.pcolormesh(
    z, t, u.T**2,
    cmap="jet",
    norm=colors.LogNorm(vmin=args.vmin),
    rasterized=True,
    shading="auto")
plot.xlim(z.min(), z.max())
plot.ylim(-1.0, +1.0)
plot.yticks([-1.0, 0.0, +1.0])
plot.xlabel(r"Distance $z$, cm")
plot.ylabel(r"Delay $t$, ps")
plot.colorbar()

plot.annotate(
    "(a)",
    xy=(0.0, 0.0),
    xytext=(+6.0, +6.0),
    xycoords="axes fraction",
    textcoords="offset points",
    color="white")


# Second and third panel: soliton parameters
plot.subplot(gs[1, 0])

da1 = a1s - a1s[0]
da2 = a2s - a2s[0]
plot.plot(z, da1, color="black", label="1")
plot.plot(z, da2, color="gray",  label="2")
plot.legend(ncol=2, loc="upper center")

plot.xlim(z.min(), z.max())

ymin = min(da1.min(), da2.min())
ymax = max(da1.max(), da2.max())
margin = 0.25 * (ymax - ymin)
plot.ylim(ymin - margin, ymax + 2 * margin)
plot.xlabel(r"Distance $z$, cm")
plot.ylabel(r"$\Delta A_n$, a.u")

plot.annotate(
    "(b.1)",
    xy=(0.0, 0.0),
    xytext=(+6.0, +6.0),
    xycoords="axes fraction",
    textcoords="offset points")

plot.subplot(gs[1, 1])

df1 = f1s - f1s[0]
df2 = f2s - f2s[0]

plot.plot(z, df1, color="black", label="1")
plot.plot(z, df2, color="gray",  label="2")
plot.legend(ncol=2, loc="upper center")

plot.xlim(z.min(), z.max())

ymin = min(df1.min(), df2.min())
ymax = max(df1.max(), df2.max())
margin = 0.25 * (ymax - ymin)

plot.ylim(ymin - margin, ymax + 2 * margin)
plot.xlabel(r"Distance $z$, cm")
plot.ylabel(r"$\Delta \omega_{n}$, rad/fs")

plot.annotate(
    "(b.2)",
    xy=(0.0, 0.0),
    xytext=(+6.0, +6.0),
    xycoords="axes fraction",
    textcoords="offset points")


# Fourth panel: input and output spectra
ax = plot.subplot(gs[2, :])
plot.plot(f, v0**0.5, color="black", linewidth=0.5, label="in",  zorder=10)
plot.plot(f, v1**0.5, color="gray",  linewidth=1.0, label="out", alpha=0.75)

plot.legend(ncol=2, loc="upper center")
plot.xlim(0.5, 4.0)
plot.ylim(0.0, 1.7)
plot.xlabel(r"Frequency $\omega$, rad/fs")
plot.ylabel(r"$\left| \tilde u(\omega) \right|^{1/2}$, a.\,u.")
ax.xaxis.set_major_locator(MultipleLocator(0.5))

plot.annotate(
    "(c)",
    xy=(1.0, 1.0),
    xytext=(-16.0, -12.0),
    xycoords="axes fraction",
    textcoords="offset points",
    bbox=dict(boxstyle="circle", pad=0.1, fc="white", ec="white"))


for rf in sum(rfs, []):
    plot.axvline(
        x=rf, color="gray",
        linewidth=0.25,
        linestyle="dotted",
        zorder=-10)

peaks = [0.90, 1.75, 2.70]
level = 0.50

for n, pf in enumerate(peaks, start=1):
    plot.annotate(
        str(n),
        xy=(pf, level),
        xytext=(-2.0, +12.0),
        textcoords="offset points",
        bbox=dict(boxstyle="circle", pad=0.15, fc="white", ec="black", linewidth=0.5),
        fontsize=XSMALL_FONT_SIZE,
        zorder=10)

plot.annotate(
    r"$i$",
    xy=(fi + 0.10, level),
    xytext=(-2.0, +12.0),
    textcoords="offset points",
    bbox=dict(boxstyle="circle", pad=0.15, fc="white", ec="black", linewidth=0.5),
    fontsize=XSMALL_FONT_SIZE,
    zorder=10)


# Fifth panel: Resonance condition
ax = plot.subplot(gs[3, :])

# Limits in k (better not to rely on the automatic ones)
kmins = [b0[(f > f1) & (f < f2)].min()]
kmaxs = [b0.max()]

# plot the scattering resonances.
ncolumns = 2

for n, rf in zip(harmonics, srfi):
    if not rf:
        continue
    plot.plot(
        f, bi * numpy.ones_like(f) + n * k0,
        color="black",
        linestyle="solid",
        label=(r"$\beta(\omega_{i})$" if n == 0 else None),
        linewidth=(1.0 if n == 0 else 0.5))

for n, rf in zip(harmonics, srf21):
    if not rf:
        continue
    plot.plot(
        f, k2 - k1 + bi + n * k0,
        color="gray",
        linestyle="dashed",
        label=(r"$k_{2} - k_{1} + \beta(\omega_{i})$" if n == 0 else None),
        linewidth=(1.0 if n == 0 else 0.5))
    kmins.append((k2 - k1 + bi + n * k0).min())
    kmaxs.append((k2 - k1 + bi + n * k0).max())

kmin = min(kmins)
kmax = max(kmaxs)
margin = 0.25 * (kmax - kmin)

pb0 = plot.plot(f, b0, color="black")
plot.xlim(0.5, 4.0)
plot.ylim(kmin - margin, kmax + 2*margin)
plot.ylabel(r"$k$, rad/$\mu$m")
plot.xlabel(r"Frequency $\omega$, rad/fs")
ax.xaxis.set_major_locator(MultipleLocator(0.5))

plot.legend(ncol=ncolumns, loc="upper center")

plot.annotate(
    "(d)",
    xy=(0.0, 0.0),
    xytext=(+16.0, +6.0),
    xycoords="axes fraction",
    textcoords="offset points",
    bbox=dict(boxstyle="circle", pad=0.1, fc="white", ec="white"))

for rf in sum(rfs, []):
    plot.axvline(
        x=rf, color="gray",
        linewidth=0.25,
        linestyle="dotted",
        zorder=-10)

level = -0.025
for n, pf in enumerate(peaks, start=1):
    idx = abs(f - pf).argmin()
    plot.annotate(
        str(n),
        xy=(pf, level),
        xytext=(-1.0, -1.0),
        textcoords="offset points",
        bbox=dict(boxstyle="circle", pad=0.15, fc="white", ec="black", linewidth=0.5),
        fontsize=XSMALL_FONT_SIZE,
        zorder=10)

plot.annotate(
    r"$i$",
    xy=(fi + 0.1, level),
    xytext=(-1.0, -1.0),
    textcoords="offset points",
    bbox=dict(boxstyle="circle", pad=0.15, fc="white", ec="black", linewidth=0.5),
    fontsize=XSMALL_FONT_SIZE,
    zorder=10)


plot.tight_layout(pad=0.25)
pr_publish(args.output)
