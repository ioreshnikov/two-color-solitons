#!/usr/bin/env python3


__doc__ = """
This script produces a resonance curve plot.
"""


import argparse
import os

from matplotlib import pyplot as plot
import numpy
from scipy.ndimage import uniform_filter1d
from scipy.optimize import curve_fit

from common.fiber import beta, beta1, beta2, gamma
from common.helpers import frame_of_reference, zeros_on_a_grid
from common.plotter import COLOR_BLACK, COLOR_GREY, pr_setup, pr_publish


parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "--interactive",
    help="show the interactive plots of the parameter averaging for every data point",
    action="store_true")
parser.add_argument(
    "seed",
    help="path to the seed soliton simulation",
    type=str)
parser.add_argument(
    "input",
    help="paths to the soliton parameters .npz files",
    type=str,
    nargs="+")
parser.add_argument(
    "output",
    help="path to the output image",
    type=str)
args = parser.parse_args()


npz = numpy.load(args.seed)

u0 = npz["u"][-1, :]
a1 = npz["a1"]; a2 = npz["a2"]  # noqa: E702
t1 = npz["t1"]; t2 = npz["t2"]  # noqa: E702
f1 = npz["f1"]; f2 = npz["f2"]  # noqa: E702


fis = []
dfs = []
das = []
for filename in args.input:
    # It's a bit hacky, but we attempt to extract the incident frequency from
    # the filename.
    basename = os.path.basename(filename)
    chunks = basename.split("_")
    _, fi, *_ = chunks
    fi = float(fi)

    npz = numpy.load(filename)

    z = npz["z"]
    f1s = npz["f1s"]
    f2s = npz["f2s"]
    a1s = npz["a1s"]
    a2s = npz["a2s"]

    # Compute frame-adjusted wavenumber of the incident radiation.
    f1 = f1s[0]
    ki = beta(fi) - beta(f1) - beta1(f1) * (fi - f1)

    nz = int(0.25 * len(z))

    f0s = uniform_filter1d(f2s, nz)
    a0s = uniform_filter1d(a2s, nz)

    df = f2s - f0s
    da = a2s - a0s

    fis.append(fi)
    dfs.append((df.max() - df.min()) / 2)
    das.append((da.max() - da.min()) / 2)

    if not args.interactive:
        continue

    plot.figure()

    plot.subplot(2, 2, 1)
    plot.plot(z / 1E4, f2s)
    plot.plot(z / 1E4, f0s)
    plot.xlabel("Distance z, cm")

    plot.subplot(2, 2, 2)
    plot.plot(z / 1E4, df)
    plot.xlabel("Distance z, cm")

    plot.subplot(2, 2, 3)
    plot.plot(z / 1E4, a2s)
    plot.plot(z / 1E4, a0s)
    plot.xlabel("Distance z, cm")

    plot.subplot(2, 2, 4)
    plot.plot(z / 1E4, da)
    plot.xlabel("Distance z, cm")

    plot.suptitle(r"$\omega_{{i}} = {:.3f}$".format(fi))
    plot.tight_layout()
    plot.show()


fis = numpy.array(fis)
dfs = numpy.array(dfs)

dfs = dfs[fis < 1.4]
fis = fis[fis < 1.4]


f = numpy.linspace(fis.min(), fis.max(), 1000)
frame = frame_of_reference(f, f1)
b = beta(f) - frame
k0sq = (
    -8 * gamma(f1) * gamma(f2) /
    (t1**2 + t2**2) ** (3/2) * (
        beta2(f1) * t1 * a1**2 + beta2(f2) * t2 * a2**2
    ))
k0 = numpy.sqrt(k0sq) * numpy.ones_like(f)
f0 = zeros_on_a_grid(f, b + k0)


pr_setup()
plot.figure(figsize=(3.2, 2.6))


plot.subplot(2, 1, 1)


# Let's try to find a Lorentzian in the left part of the plot
def lorentzian(x, *params):  # noqa: D401
    """A Lorentzian distribution."""
    a, x0, dx = params
    return a / (1 + (x - x0)**2 / dx**2)

params, _ = curve_fit(lorentzian, fis, dfs, p0=[0.1, 1.0, 0.5])
af, xf, dxf = params


plot.scatter(fis, dfs, marker="x", s=10, label="Simulation", color=COLOR_BLACK)
plot.plot(
    f, lorentzian(f, af, xf, dxf),
    label="Lorentzian fit",
    color=COLOR_GREY)

plot.axvline(
    f0,
    color=COLOR_GREY,
    linewidth=0.33,
    linestyle="dotted",
    zorder=-10)
plot.xlim(f.min(), f.max())
plot.ylim(0.00, 0.03)
plot.xlabel(r"Incident frequency $\omega_{i}$, rad/fs")
plot.ylabel(r"max. $\Delta \omega_{2}$, rad/fs")
plot.legend()

plot.annotate(
    "(a)",
    xy=(1.0, 0.0),
    xytext=(-16.0, +12.0),
    xycoords="axes fraction",
    textcoords="offset points",
    bbox=dict(boxstyle="circle", pad=0.1, fc="white", ec="white"))

plot.subplot(2, 1, 2)

plot.plot(f, b, label=r"$\beta(\omega)$", color=COLOR_BLACK)
plot.plot(f, +k0, label=r"$\pm K_0$", color=COLOR_GREY, linestyle="dotted")
plot.plot(f, -k0, color=COLOR_GREY, linestyle="dotted")
plot.axvline(
    f0,
    color=COLOR_GREY,
    linewidth=0.33,
    linestyle="dotted",
    zorder=-10)
plot.xlim(f.min(), f.max())
plot.ylim(-0.01, +0.01)
plot.xlabel(r"Frequency $\omega$, rad/fs")
plot.ylabel(r"$k$, rad/$\mu$m")
plot.legend(ncol=2)

plot.annotate(
    "(b)",
    xy=(1.0, 0.0),
    xytext=(-16.0, +6.0),
    xycoords="axes fraction",
    textcoords="offset points",
    bbox=dict(boxstyle="circle", pad=0.1, fc="white", ec="white"))


plot.tight_layout()
pr_publish(args.output)
