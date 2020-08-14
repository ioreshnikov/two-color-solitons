#!/usr/bin/env python3


import argparse

from matplotlib import pyplot as plot
import numpy

from common.fiber import (
    beta, gv_matching_frequencies,
    fundamental_soliton_dispersive_relation)


parser = argparse.ArgumentParser()
parser.add_argument(
    "-f1",
    help="first soliton frequency",
    type=float,
    required=True)
parser.add_argument(
    "-t1",
    help="first soliton width",
    type=float,
    required=True)
parser.add_argument(
    "-f2",
    help="second soliton frequency",
    type=float)
parser.add_argument(
    "-t2",
    help="second soliton width",
    type=float)
args = parser.parse_args()


f = numpy.linspace(0.5, 6, 100)

if not args.f2:
    fs = gv_matching_frequencies(args.f1)
    if args.f1 < 2.0:
        args.f2 = fs[1]
    else:
        args.f2 = fs[0]

if not args.t2:
    args.t2 = args.t1

b = beta(f)
k1 = fundamental_soliton_dispersive_relation(args.f1, args.t1, f)
k2 = fundamental_soliton_dispersive_relation(args.f2, args.t2, f)

plot.plot(f, k1 - b)
plot.plot(f, k2 - b)
plot.plot(f, numpy.zeros_like(f), color="gray", linewidth=0.5)
plot.gca().axvline(args.f1, color="gray", linewidth=0.5)
plot.gca().axvline(args.f2, color="gray", linewidth=0.5)
plot.xlabel(r"$\omega$, rad/fs")
plot.show()
