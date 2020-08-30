#!/usr/bin/env python3


import argparse

from matplotlib import pyplot as plot
from matplotlib.widgets import TextBox
import numpy


from common.fiber import (
    beta, beta1, gv_matching_frequencies,
    fundamental_soliton_dispersive_relation)
from common.helpers import zeros_on_a_grid


parser = argparse.ArgumentParser()
parser.add_argument(
    "-f1",
    help="carrier frequency of the first soliton",
    type=float,
    default=1.000)
parser.add_argument(
    "-t0",
    help="width of the fundamental soliton",
    type=float,
    default=20.0)
parser.add_argument(
    "-fi",
    help="carrier frequency of the incident dispersive wave",
    type=float)
args = parser.parse_args()


f = numpy.linspace(0.5, 4, 1000)


# Set up the plots and the text boxes
fig, ax = plot.subplots()
plot.subplots_adjust(bottom=0.25)

f1vline = plot.axvline(0, color="gray", linewidth=1.0)
f2vline = plot.axvline(0, color="gray", linewidth=1.0)
gvmvline = plot.axvline(0, color="gray", linewidth=1.0, linestyle="dotted")

bplot, *_ = plot.plot(
    f, numpy.zeros_like(f),
    label=r"$\beta(\omega)$",
    color="#9400d3")
n1plot, *_ = plot.plot(
    f, numpy.zeros_like(f),
    label=r"$\eta_{1}(\omega)$",
    color="#e51e10")
n2plot, *_ = plot.plot(
    f, numpy.zeros_like(f),
    label=r"$\eta_{2}(\omega)$",
    color="#0072b2")

cf1scatter = plot.scatter([], [], s=50, color="#e51e10", zorder=10)
cf2scatter = plot.scatter([], [], s=50, color="#0072b2", zorder=10)

cf12plot, *_ = plot.plot(
    f, numpy.zeros_like(f),
    label=r"$2 \eta_{1}(\omega) - \eta_{2}(\omega)$",
    color="#e51e10",
    linestyle="dotted")
cf12scatter = plot.scatter([], [], s=10, color="#e51e10", zorder=10)

cf21plot, *_ = plot.plot(
    f, numpy.zeros_like(f),
    label=r"$2 \eta_{2}(\omega) - \eta_{1}(\omega)$",
    color="#0072b2", linestyle="dotted")
cf21scatter = plot.scatter([], [], s=10, color="#0072b2", zorder=10)

biplot, *_ = plot.plot(
    f, numpy.zeros_like(f),
    color="#009e73", label=r"$\beta(\omega_{inc})$")
fivline = plot.axvline(0, color="#009e73", linestyle="dotted")

rf0scatter = plot.scatter([], [], s=50, color="#009e73", zorder=10)

rf1plot, *_ = plot.plot(
    f, numpy.zeros_like(f),
    color="#e69f00",
    label=r"$\beta(\omega) = 2 \eta_{1}(\omega) - \beta(\omega_{inc})$")
rf1scatter = plot.scatter([], [], s=50, color="#e69f00", zorder=10)
rf2plot, *_ = plot.plot(
    f, numpy.zeros_like(f),
    color="#f0e442",
    label=r"$\beta(\omega) = 2 \eta_{2}(\omega) - \beta(\omega_{inc})$")
rf2scatter = plot.scatter([], [], s=50, color="#f0e442", zorder=10)
rf12plot, *_ = plot.plot(
    f, numpy.zeros_like(f),
    color="#56b4e9",
    label=r"$\beta(\omega) = \eta_{1}(\omega) + \eta_{2}(\omega) - \beta(\omega_{inc})$")
rf12scatter = plot.scatter([], [], s=50, color="#56b4e9", zorder=10)

plot.xlim(f.min(), f.max())
plot.ylim(-0.2, 0.1)
plot.xlabel(r"$\omega$, rad/fs")
plot.legend()


def update(*_, **__):
    f1 = f1box.text
    fi = fibox.text

    f1 = float(f1)
    fi = float(fi) if fi else None

    # This dispersive relation defines a frame of reference co-moving
    # with the first soliton.
    frame = beta(f1) + beta1(f1) * (f - f1)

    fs = gv_matching_frequencies(f1)
    gvm, f2 = fs

    b = beta(f) - frame
    n1 = fundamental_soliton_dispersive_relation(f1, args.t0, f) - frame
    n2 = fundamental_soliton_dispersive_relation(f2, args.t0, f) - frame

    f1vline.set_xdata(f1)
    f2vline.set_xdata(f2)
    gvmvline.set_xdata(gvm)

    bplot.set_ydata(b)
    n1plot.set_ydata(n1)
    n2plot.set_ydata(n2)

    # Find the Cherenkov resonances due to the individual solitons.
    cf1 = zeros_on_a_grid(f, b - n1)
    cf2 = zeros_on_a_grid(f, b - n2)

    offsets = numpy.zeros((len(cf1), 2))
    offsets[:, 0] = cf1
    offsets[:, 1] = n1[-1]
    cf1scatter.set_offsets(offsets)

    offsets = numpy.zeros((len(cf2), 2))
    offsets[:, 0] = cf2
    offsets[:, 1] = n2[-1]
    cf2scatter.set_offsets(offsets)

    # Find the Cherenkov resonances due to FWVM process.
    cf12 = zeros_on_a_grid(f, 2 * n1 - n2 - b)
    cf21 = zeros_on_a_grid(f, 2 * n2 - n1 - b)

    offsets = numpy.zeros((len(cf12), 2))
    offsets[:, 0] = cf12
    offsets[:, 1] = 2 * n1[-1] - n2[-1]
    cf12scatter.set_offsets(offsets)

    offsets = numpy.zeros((len(cf21), 2))
    offsets[:, 0] = cf21
    offsets[:, 1] = 2 * n2[-1] - n1[-1]
    cf21scatter.set_offsets(offsets)

    if len(cf12):
        cf12plot.set_visible(True)
        cf12plot.set_ydata(2 * n1 - n2)
    else:
        cf12plot.set_visible(False)

    if len(cf21):
        cf21plot.set_visible(True)
        cf21plot.set_ydata(2 * n2 - n1)
    else:
        cf21plot.set_visible(False)

    if not fi:
        biplot.set_visible(False)
        rf1plot.set_visible(False)
        rf2plot.set_visible(False)
        rf12plot.set_visible(False)
        return

    fivline.set_xdata(fi)

    bi = beta(fi) - beta(f1) - beta1(f1) * (fi - f1)
    rf0 = numpy.array([
        f0 for f0 in zeros_on_a_grid(f, b - bi)
        if 0.1 < abs(f0 - fi) < 1.0
    ])
    if len(rf0):
        biplot.set_visible(True)
        biplot.set_ydata(bi)
    else:
        biplot.set_visible(False)

    offsets = numpy.zeros((len(rf0), 2))
    offsets[:, 0] = rf0
    offsets[:, 1] = bi
    rf0scatter.set_offsets(offsets)

    rf1 = zeros_on_a_grid(f, 2 * n1 - bi - b)
    if len(rf1):
        rf1plot.set_visible(True)
        rf1plot.set_ydata(2 * n1 - bi)
    else:
        rf1plot.set_visible(False)

    offsets = numpy.zeros((len(rf1), 2))
    offsets[:, 0] = rf1
    offsets[:, 1] = 2 * n1[-1] - bi
    rf1scatter.set_offsets(offsets)

    rf2 = zeros_on_a_grid(f, 2 * n2 - bi - b)
    if len(rf2):
        rf2plot.set_visible(True)
        rf2plot.set_ydata(2 * n2 - bi)
    else:
        rf2plot.set_visible(False)

    offsets = numpy.zeros((len(rf2), 2))
    offsets[:, 0] = rf2
    offsets[:, 1] = 2 * n2[-1] - bi
    rf2scatter.set_offsets(offsets)

    rf12 = zeros_on_a_grid(f, n1 + n2 - bi - b)
    if len(rf12):
        rf12plot.set_visible(True)
        rf12plot.set_ydata(n1 + n2 - bi)
    else:
        rf12plot.set_visible(False)

    offsets = numpy.zeros((len(rf12), 2))
    offsets[:, 0] = rf12
    offsets[:, 1] = 2 * n1[-1] + n2[-1] - bi
    rf12scatter.set_offsets(offsets)


f1ax = plot.axes([0.12, 0.05, 0.35, 0.075])
f1box = TextBox(f1ax, r"$\omega_{1}$", initial=args.f1)
f1box.on_submit(update)

fiax = plot.axes([0.55, 0.05, 0.35, 0.075])
fibox = TextBox(fiax, r"$\omega_{inc}$", initial=args.fi)
fibox.on_submit(update)

update()

plot.show()
