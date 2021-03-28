#!/usr/bin/env python3


from argparse import ArgumentParser

from matplotlib import pyplot as plot
from matplotlib import colors
from matplotlib.ticker import MultipleLocator
from numpy import fft
import numpy

from common.fiber import (
    beta, beta1, gamma,
    estimate_soliton_parameters,
    fundamental_soliton_amplitude,
    fundamental_soliton_dispersive_relation)
from common.helpers import zeros_on_a_grid
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

f = fft.fftfreq(len(t), t[1] - t[0])
f = 2 * numpy.pi * fft.fftshift(f)

u = npz["u"]

# From the full output spectrum we estimate the final soliton
# parameters using the original values as an estimate
u0 = u[-1, :]
a1 = npz["a1"]; a2 = npz["a2"]
t1 = npz["t1"]; t2 = npz["t2"]
f1 = npz["f1"]; f2 = npz["f2"]
a1, a2, t1, t2, f1, f2 = estimate_soliton_parameters(
    t, u0, a1, a2, t1, t2, f1, f2)


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

a1 = fundamental_soliton_amplitude(f1, t1)
a2 = fundamental_soliton_amplitude(f2, t2)

k1 = fundamental_soliton_dispersive_relation(f1, t1, f) - frame + gamma(f1) * a2**2
k2 = fundamental_soliton_dispersive_relation(f2, t2, f) - frame + gamma(f2) * a1**2
b0 = beta(f) - frame

k12_too_close = abs(k1 - k2).max() < 0.01
# ^^^ for a almost perfect soliton all those resonances will fall
# very close to the first kind of resonance, so there is no point
# in plotting them.

rfs = []
if "fi" not in npz:
    # Resonance frequencies for Cherenkov radiation as predicted by
    # the theory:
    # 1. Cherenkov radiation due to first soliton component.
    crf1 = zeros_on_a_grid(f, b0 - k1)
    # 2. Cherenkov radiation due to second soliton component.
    crf2 = zeros_on_a_grid(f, b0 - k2)
    # 3. Cherenkov radiation due to FWVM processes between the components.
    crf12 = zeros_on_a_grid(f, b0 - 2*k1 + k2)
    crf21 = zeros_on_a_grid(f, b0 - 2*k2 + k1)

    rfs.extend(crf1)
    rfs.extend(crf2)
    rfs.extend(crf12)
    rfs.extend(crf21)
else:
    # If it's a scattering problem, then we are interested only in the
    # scattering resonances.
    fi = npz["fi"]
    framei = beta(f1) + beta1(f1) * (fi - f1)
    bi = beta(fi) - framei

    # 1. Resonance due to β(ωₛ) = β(ωᵢ)
    srfi = zeros_on_a_grid(f, b0 - bi)

    # 2. Resonances due to β(ω) = k₂(ω) - k₁(ω) + β(ωᵢ)
    srf21 = zeros_on_a_grid(f, b0 - k2 + k1 - bi)

    rfs.extend(srfi)
    if not k12_too_close:
        rfs.extend(srf21)


# Finally, we plot everything.
pr_setup()

plot.figure(figsize=(3.4, 4.2))

plot.subplot(3, 1, 1)
plot.pcolormesh(
    z, t, u.T**2,
    cmap="magma",
    norm=colors.LogNorm(vmin=1E-6),
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
    xy=(1.0, 0.0),
    xytext=(-18.0, +10.0),
    xycoords="axes fraction",
    textcoords="offset points",
    color="white")

ax = plot.subplot(3, 1, 2)
plot.plot(f, v0**0.5, color="gray",  linewidth=1.5,  label="in", zorder=-10)
plot.plot(f, v1**0.5, color="black", linewidth=0.5, label="out")

for rf in rfs:
    plot.axvline(
        x=rf, color="gray",
        linewidth=0.25,
        linestyle="dotted",
        zorder=-10)

plot.legend(ncol=2, loc="upper center")
plot.xlim(0.5, 4.0)
plot.ylim(0.0, 1.1)
plot.xlabel(r"Frequency $\omega$, rad/fs")
plot.ylabel(r"$\left| \tilde u(\omega) \right|^{1/2}$, a.\,u.")
ax.xaxis.set_major_locator(MultipleLocator(0.5))

plot.annotate(
    "(b)",
    xy=(1.0, 1.0),
    xytext=(-18.0, -12.0),
    xycoords="axes fraction",
    textcoords="offset points")

ax = plot.subplot(3, 1, 3)

for rf in rfs:
    plot.axvline(
        x=rf, color="gray",
        linewidth=0.25,
        linestyle="dotted",
        zorder=-10)

# Limits in k (better not to rely on the automatic ones)
kmins = [k1.min(), k2.min(), b0[(f > f1) & (f < f2)].min()]
kmaxs = [k1.max(), k2.max(), b0.max()]

if "fi" not in npz:
    # Again, if it's not a scattering problem, we deal only with
    # Cherenkov radiation.
    plot.plot(f, k1, color="black", linestyle="solid", label="$k_{1}$")
    plot.plot(f, k2, color="black", linestyle="dashed", label="$k_{2}$")

    ncolumns = 2
    if any((2*k2 - k1) < b0):
        ncolumns += 1
        plot.plot(
            f, 2*k2 - k1,
            color="black",
            linestyle="dotted",
            label="$2k_{2} - k_{1}$")
    if any((2*k1 - k2) < b0):
        ncolumns += 1
        plot.plot(
            f, 2*k1 - k2,
            color="gray",
            linestyle="dotted",
            label="$2k_{1} - k_{2}$")

    kmins.extend([
        (2*k1 - k2).min(),
        (2*k2 - k1).min()
    ])
else:
    # In case of a scattering problem we plot only the scattering
    # resonances.
    ncolumns = 0

    if len(srfi) > 1:
        ncolumns += 1
        plot.plot(
            f, bi * numpy.ones_like(f),
            color="black",
            linestyle="solid",
            label=r"$\beta(\omega_{i})$")

    # if len(srf12) > 0 and not k12_too_close:
    #     ncolumns += 1
    #     plot.plot(
    #         f, k1 - k2 + bi,
    #         color="black",
    #         linestyle="dashed",
    #         label=r"$k_{1} - k_{2} + \beta(\omega_{i})$")

    if len(srf21) > 0 and not k12_too_close:
        ncolumns += 1
        plot.plot(
            f, k2 - k1 + bi,
            color="gray",
            linestyle="dashed",
            label=r"$k_{2} - k_{1} + \beta(\omega_{i})$")

    kmins.extend([
        # (k1 - k2 + bi).min(),
        (k2 - k1 + bi).min(),
    ])

kmin = min(kmins)
kmax = max(kmaxs)
margin = 0.25 * (kmax - kmin)

plot.plot(f, b0, color="black")
plot.xlim(0.5, 4.0)
plot.ylim(kmin - margin, kmax + 2*margin)
plot.xlabel(r"Frequency $\omega$, rad/fs")
ax.xaxis.set_major_locator(MultipleLocator(0.5))
plot.legend(ncol=ncolumns, loc="upper center")

plot.annotate(
    "(c)",
    xy=(1.0, 1.0),
    xytext=(-18.0, -12.0),
    xycoords="axes fraction",
    textcoords="offset points")

plot.tight_layout(pad=1.0)
pr_publish(args.output)
