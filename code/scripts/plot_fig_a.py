#!/usr/bin/env python3


__doc__ = """
This script plots dispersive operator β(ω) together with it's direvatives.
"""


import argparse

import numpy
from matplotlib import pyplot as plot
from matplotlib import ticker

from common.fiber import beta, beta1, beta2
from common.helpers import zeros_on_a_grid
from common.plotter import (
    COLOR_BLACK, COLOR_GREY,
    XSMALL_FONT_SIZE,
    pr_setup, pr_publish)


parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "--interactive",
    help="show the interactive plot without producing an image",
    action="store_true")
parser.add_argument(
    "output",
    help="path to the output image",
    type=str)
args = parser.parse_args()


f = numpy.linspace(0.5, 3.5, 2**10)
b1 = beta1(f)
b2 = beta2(f)

fa1, fa2 = zeros_on_a_grid(f, b2)
f1, _, f2 = zeros_on_a_grid(f, beta1(f) - beta1(1.010))

if not args.interactive:
    pr_setup()

plot.figure(figsize=(3.6, 3.6))

ax = plot.subplot(2, 1, 1)

plot.plot(f, 1/b1 * 10**3, color=COLOR_BLACK)
plot.plot(
    f, 1/beta1(f1) * 10**3 * numpy.ones_like(f),
    linewidth=0.5,
    color=COLOR_GREY,
    zorder=-10)
plot.scatter(
    [f1, f2],
    [1/beta1(f1) * 10**3, 1/beta1(f2) * 10**3],
    s=70, marker='.', color=COLOR_BLACK)

top_axis = ax.secondary_xaxis("top")
top_axis.xaxis.set_major_locator(ticker.FixedLocator((f1, f2)))
top_axis.xaxis.set_major_formatter(
    ticker.FixedFormatter((r"$\omega_1$", r"$\omega_2$")))

plot.axvline(
    x=f1, color=COLOR_GREY,
    linewidth=0.50,
    zorder=-10)
plot.axvline(
    x=f2, color=COLOR_GREY,
    linewidth=0.50,
    zorder=-10)
plot.axvline(
    x=fa1, color=COLOR_GREY,
    linewidth=0.33,
    linestyle="dotted",
    zorder=-10)
plot.axvline(
    x=fa2, color=COLOR_GREY,
    linewidth=0.33,
    linestyle="dotted",
    zorder=-10)
plot.axvspan(fa1, fa2, alpha=0.25, color=COLOR_GREY)

plot.annotate(
    r"$A_{1}$",
    xy=(1.0, 76.0),
    xytext=(0.0, -20.0),
    textcoords="offset points",
    bbox=dict(boxstyle="circle", pad=0.15, fc="white", ec="black", linewidth=0.5),
    fontsize=XSMALL_FONT_SIZE)
plot.annotate(
    r"$N$",
    xy=(2.0, 76.0),
    xytext=(0.0, -20.0),
    textcoords="offset points",
    bbox=dict(boxstyle="circle", pad=0.30, fc="white", ec="black", linewidth=0.5),
    fontsize=XSMALL_FONT_SIZE)
plot.annotate(
    r"$A_{2}$",
    xy=(3.0, 76.0),
    xytext=(0.0, -20.0),
    textcoords="offset points",
    bbox=dict(boxstyle="circle", pad=0.15, fc="white", ec="black", linewidth=0.5),
    fontsize=XSMALL_FONT_SIZE)

plot.xlim(0.5, 3.5)
plot.ylim(73.0, 76.0)
plot.xlabel(r"Frequency $\omega$ (rad/fs)")
plot.ylabel(r"$v_g$ (nm/fs)")

ax = plot.subplot(2, 1, 2)
plot.plot(f, b2, color=COLOR_BLACK)
plot.plot(
    f, numpy.zeros_like(f),
    linewidth=0.5,
    color=COLOR_GREY,
    zorder=-10)

plot.axvline(
    x=fa1, color=COLOR_GREY,
    linewidth=0.33,
    linestyle="dotted",
    zorder=-10)
plot.axvline(
    x=fa2, color=COLOR_GREY,
    linewidth=0.33,
    linestyle="dotted",
    zorder=-10)
plot.axvspan(fa1, fa2, alpha=0.25, color=COLOR_GREY)

plot.annotate(
    r"$A_{1}$",
    xy=(1.0, 0.5),
    xytext=(0.0, -20.0),
    textcoords="offset points",
    bbox=dict(boxstyle="circle", pad=0.15, fc="white", ec="black", linewidth=0.5),
    fontsize=XSMALL_FONT_SIZE)
plot.annotate(
    r"$N$",
    xy=(2.0, 0.5),
    xytext=(0.0, -20.0),
    textcoords="offset points",
    bbox=dict(boxstyle="circle", pad=0.30, fc="white", ec="black", linewidth=0.5),
    fontsize=XSMALL_FONT_SIZE)
plot.annotate(
    r"$A_{2}$",
    xy=(3.0, 0.5),
    xytext=(0.0, -20.0),
    textcoords="offset points",
    bbox=dict(boxstyle="circle", pad=0.15, fc="white", ec="black", linewidth=0.5),
    fontsize=XSMALL_FONT_SIZE)

plot.xlim(0.5, 3.5)
plot.ylim(-0.5, 0.5)
plot.xlabel(r"Frequency $\omega$ (rad/fs)")
plot.ylabel(r"$\beta''(\omega)$, (fs$^2$/$\mu$m)")

plot.tight_layout()

if args.interactive:
    plot.show()
else:
    pr_publish(args.output)
