#!/usr/bin/env python3


__doc__ = """
This script produces a three-panel plot of the simulation results. The
first panel is the time domain view, the second one is the input
and the output spectrum of the simulation and the third one is the
resonance condition diagram.
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
    fundamental_soliton_dispersive_relation,
    peaks_close_to, zeros_on_a_grid)
from common.plotter import XSMALL_FONT_SIZE, pr_setup, pr_publish


parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawTextHelpFormatter)
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

# Read the incident frequency (if present)
fi = None
if "fi" in npz:
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

k12_too_close = abs(k1 - k2).max() < 0.01
# ^^^ for a almost perfect soliton all those resonances will fall
# very close to the first kind of resonance, so there is no point
# in plotting them.

rfs = []
if not fi:
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
    framei = beta(f1) + beta1(f1) * (fi - f1)
    bi = beta(fi) - framei

    # 1. Resonance due to β(ωₛ) = β(ωᵢ) (but exclude ωᵢ)
    srfi = [
        rf for rf in zeros_on_a_grid(f, b0 - bi)
        if abs(rf - fi) > 0.01
    ]

    # 2. Resonances due to β(ω) = k₂(ω) - k₁(ω) + β(ωᵢ)
    srf21 = zeros_on_a_grid(f, b0 - k2 + k1 - bi)

    # 3. Resonances due to β(ω) = k₁(ω) - k₂(ω) + β(ωᵢ)
    srf12 = zeros_on_a_grid(f, b0 - k1 + k2 - bi)

    rfs.extend(srfi)
    if not k12_too_close:
        rfs.extend(srf21)
        rfs.extend(srf12)


# Filter the resonance frequencies that correspond to a peak in the
# output spectrum.
peaks = peaks_close_to(f, v1, rfs)


# Finally, we plot everything.
pr_setup()

height_ratios = [2/3, 2/3, 1]
plot.figure(
    figsize=(3.2, sum(1.4 * h for h in height_ratios)))
gs = gridspec.GridSpec(3, 1, height_ratios=height_ratios)


# First panel (optional): time-domain plot
plot.subplot(gs[0])
plot.pcolormesh(
    z, t, u.T**2,
    # cmap="inferno",
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


# Second panel: input and output spectra
ax = plot.subplot(gs[1])
plot.plot(f, v0**0.5, color="black", linewidth=0.5, label="in",  zorder=10)
plot.plot(f, v1**0.5, color="gray",  linewidth=1.0, label="out", alpha=0.75)

plot.legend(ncol=2, loc="upper center")
plot.xlim(0.5, 4.0)
plot.ylim(0.0, 1.4)
plot.xlabel(r"Frequency $\omega$, rad/fs")
plot.ylabel(r"$\left| \tilde u(\omega) \right|^{1/2}$, a.\,u.")
ax.xaxis.set_major_locator(MultipleLocator(0.5))

plot.annotate(
    "(b)",
    xy=(1.0, 1.0),
    xytext=(-16.0, -12.0),
    xycoords="axes fraction",
    textcoords="offset points",
    bbox=dict(boxstyle="circle", pad=0.1, fc="white", ec="white"))


for rf in rfs:
    plot.axvline(
        x=rf, color="gray",
        linewidth=0.25,
        linestyle="dotted",
        zorder=-10)

level = max(pv for _, pv in peaks)
if fi:
    idx = abs(f - fi).argmin()
    level = max(level, v1[idx])

for n, (pf, _) in enumerate(peaks, start=1):
    plot.annotate(
        str(n),
        xy=(pf, level),
        xytext=(-2.0, +12.0),
        textcoords="offset points",
        bbox=dict(boxstyle="circle", pad=0.15, fc="white", ec="black", linewidth=0.5),
        fontsize=XSMALL_FONT_SIZE)

if fi:
    plot.annotate(
        r"$i$",
        xy=(fi, level),
        xytext=(-2.0, +12.0),
        textcoords="offset points",
        bbox=dict(boxstyle="circle", pad=0.15, fc="white", ec="black", linewidth=0.5),
        fontsize=XSMALL_FONT_SIZE)


# Third panel (optional): Resonance condition
ax = plot.subplot(gs[2])

# Limits in k (better not to rely on the automatic ones)
kmins = [b0[(f > f1) & (f < f2)].min()]
kmaxs = [b0.max()]

if not fi:
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
        k1.min(),
        k2.min(),
        (2*k1 - k2).min(),
        (2*k2 - k1).min()
    ])
    kmaxs.extend([
        k1.max(),
        k2.max(),
        (2*k1 - k2).max(),
        (2*k2 - k1).max()
    ])
else:
    # In case of a scattering problem we plot only the scattering
    # resonances.
    ncolumns = 0

    if srfi:
        ncolumns += 1
        plot.plot(
            f, bi * numpy.ones_like(f),
            color="black",
            linestyle="solid",
            label=r"$\beta(\omega_{i})$")

    if len(srf21) and not k12_too_close:
        ncolumns += 1
        plot.plot(
            f, k2 - k1 + bi,
            color="gray",
            linestyle="dashed",
            label=r"$k_{2} - k_{1} + \beta(\omega_{i})$")
        kmins.append((k2 - k1 + bi).min())
        kmaxs.append((k2 - k1 + bi).max())

    if len(srf12) and not k12_too_close:
        ncolumns += 1
        plot.plot(
            f, k1 - k2 + bi,
            color="gray",
            linestyle="dotted",
            label=r"$k_{1} - k_{2} + \beta(\omega_{i})$")
        kmins.append((k1 - k2 + bi).min())
        kmaxs.append((k1 - k2 + bi).max())

kmin = min(kmins)
kmax = max(kmaxs)
margin = 0.25 * (kmax - kmin)

pb0 = plot.plot(f, b0, color="black")
plot.xlim(0.5, 4.0)
plot.ylim(kmin - margin, kmax + 2*margin)
plot.ylabel(r"$k$, rad/$\mu$m")
plot.xlabel(r"Frequency $\omega$, rad/fs")
ax.xaxis.set_major_locator(MultipleLocator(0.5))

if fi and len(srf21) and len(srf12) and not k12_too_close:
    # Ok, this is a bit more complicated :) If we try to draw the
    # legend as usual then the labels will not fit. What we want to do
    # is to combine the labels for FWVM resonances into one (with ±
    # and ∓ signs) and stack the line markers according to the signs.
    # This requires a separate legend handler. I am following this answer:
    #
    #     https://stackoverflow.com/a/41765095

    class Stub:
        pass

    class FWVMHanlder(HandlerBase):
        def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
            line1 = plot.Line2D(
                [x0, x0 + width],
                [0.7 * height, 0.7 * height],
                color="gray",
                linestyle="dashed")
            line2 = plot.Line2D(
                [x0, x0 + width],
                [0.3 * height, 0.3 * height],
                color="gray",
                linestyle="dotted")
            return [line1, line2]

    plot.legend(
        [pb0[0], Stub],
        [r"$\beta(\omega_{i})$", r"$\pm k_{2} \mp k_{1} + \beta(\omega_{i})$"],
        ncol=ncolumns,
        loc="upper center",
        handler_map={Stub: FWVMHanlder()})
else:
    plot.legend(ncol=ncolumns, loc="upper center")

plot.annotate(
    "(c)",
    xy=(0.0, 0.0),
    xytext=(+16.0, +6.0),
    xycoords="axes fraction",
    textcoords="offset points",
    bbox=dict(boxstyle="circle", pad=0.1, fc="white", ec="white"))

for rf in rfs:
    plot.axvline(
        x=rf, color="gray",
        linewidth=0.25,
        linestyle="dotted",
        zorder=-10)

for n, (pf, _) in enumerate(peaks, start=1):
    idx = abs(f - pf).argmin()
    level = b0[idx]
    plot.annotate(
        str(n),
        xy=(pf, level),
        xytext=(-1.0, -1.0),
        textcoords="offset points",
        bbox=dict(boxstyle="circle", pad=0.15, fc="white", ec="black", linewidth=0.5),
        fontsize=XSMALL_FONT_SIZE)

if fi:
    idx = abs(f - fi).argmin()
    level = b0[idx]
    plot.annotate(
        r"$i$",
        xy=(fi, level),
        xytext=(-1.0, -1.0),
        textcoords="offset points",
        bbox=dict(boxstyle="circle", pad=0.15, fc="white", ec="black", linewidth=0.5),
        fontsize=XSMALL_FONT_SIZE)


plot.tight_layout(pad=0.25)
pr_publish(args.output)
