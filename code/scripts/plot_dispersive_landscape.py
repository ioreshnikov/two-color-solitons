#!/usr/bin/env python3


import argparse

from matplotlib import pyplot as plot
import numpy

from common.fiber import beta, beta1, beta2, gv_matching_frequencies


parser = argparse.ArgumentParser()
parser.add_argument(
    "-f0",
    help="central frequency for slope compensation",
    type=float)
args = parser.parse_args()


def compensated_beta(f):
    if args.f0:
        return beta(f) - beta1(args.f0) * (f - args.f0)
    else:
        return beta(f)


f = numpy.linspace(0.5, 6, 100)
if args.f0:
    fs = gv_matching_frequencies(args.f0)


plot.figure(figsize=(3, 6))

plot.subplot(3, 1, 1)
plot.plot(f, compensated_beta(f))
if args.f0:
    plot.gca().axvline(args.f0, color="gray", linewidth=0.5)
    plot.scatter(args.f0, compensated_beta(args.f0))
    for f1 in fs:
        plot.gca().axvline(f1, color="gray", linewidth=0.5)
        plot.scatter(f1, compensated_beta(f1))
    plot.xticks(list(plot.xticks()[0]) + [args.f0] + list(fs))
plot.xlim(f.min(), f.max())
plot.xlabel(r"Angular frequency $\omega$, rad/fs")
plot.ylabel(r"$\beta(\omega)$")

plot.subplot(3, 1, 2)
plot.plot(f, beta1(f))
if args.f0:
    plot.plot(f, beta1(args.f0) * numpy.ones_like(f), color="gray", linewidth=0.5)
    plot.gca().axvline(args.f0, color="gray", linewidth=0.5)
    plot.scatter(args.f0, beta1(args.f0))
    for f1 in fs:
        plot.gca().axvline(f1, color="gray", linewidth=0.5)
        plot.scatter(f1, beta1(f1))
    plot.xticks(list(plot.xticks()[0]) + [args.f0] + list(fs))
plot.xlim(f.min(), f.max())
plot.xlabel(r"Angular frequency $\omega$, rad/fs")
plot.ylabel(r"$\beta_{1}(\omega)$")

plot.subplot(3, 1, 3)
plot.plot(f, numpy.zeros_like(f), color="gray", linewidth=0.5)
plot.plot(f, beta2(f))
if args.f0:
    plot.gca().axvline(args.f0, color="gray", linewidth=0.5)
    plot.scatter(args.f0, beta2(args.f0))
    for f1 in fs:
        plot.gca().axvline(f1, color="gray", linewidth=0.5)
        plot.scatter(f1, beta2(f1))
    plot.xticks(list(plot.xticks()[0]) + [args.f0] + list(fs))
plot.xlim(f.min(), f.max())
plot.xlabel(r"Angular frequency $\omega$, rad/fs")
plot.ylabel(r"$\beta_{2}(\omega)$")

plot.tight_layout()
plot.show()
