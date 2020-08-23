#!/usr/bin/env python3


import argparse

from matplotlib import pyplot as plot
from matplotlib.widgets import Slider
import numpy

from common.fiber import (
    beta, gv_matching_frequencies,
    fundamental_soliton_dispersive_relation)
from common.helpers import zeros_on_a_grid


parser = argparse.ArgumentParser()
args = parser.parse_args()


f1 = 1.0
t1 = 20


fig, axs = plot.subplots()
plot.subplots_adjust(bottom=0.25)

f_ = numpy.arange(0.5, 5, 0.01)

f2 = gv_matching_frequencies(f1)[-1]

b0 = beta(f_)
n1 = fundamental_soliton_dispersive_relation(f1, t1, f_)
n2 = fundamental_soliton_dispersive_relation(f2, t1, f_)

p1, *_ = plot.plot(f_, b0 - n1, color="#d62728")
p2, *_ = plot.plot(f_, b0 - n2, color="#1f77b4")

res1 = zeros_on_a_grid(f_, b0 - n1)
res2 = zeros_on_a_grid(f_, b0 - n2)

s1 = plot.scatter(res1, 0 * res1, color="#d62728")
s2 = plot.scatter(res2, 0 * res2, color="#1f77b2")

plot.plot(f_, numpy.zeros_like(f_), color="gray", linewidth=0.5)
v1 = plot.axvline(f1, color="gray", linewidth=0.5)
v2 = plot.axvline(f2, color="gray", linewidth=0.5)

plot.xlim(f_.min(), f_.max())
plot.xlabel(r"$\omega$, rad/fs")
plot.ylabel(r"$k$, rad/$\mu$m")
plot.ylim(-0.5, +0.5)


axfreq = plot.axes([0.10, 0.10, 0.80, 0.02])
axwidth = plot.axes([0.10, 0.15, 0.80, 0.02])

sfreq = Slider(
    axfreq, r"$\omega_{1}$, rad/fs",
    1.0, 2.5, valinit=f1, valstep=0.01)
swidth = Slider(
    axwidth, r"$t_{1}$, fs",
    10, 50, valinit=t1, valstep=1)


def update(*_, **__):
    f1 = sfreq.val
    t1 = swidth.val

    f2 = gv_matching_frequencies(f1)[-1]

    b0 = beta(f_)
    n1 = fundamental_soliton_dispersive_relation(f1, t1, f_)
    n2 = fundamental_soliton_dispersive_relation(f2, t1, f_)

    res1 = zeros_on_a_grid(f_, b0 - n1)
    res2 = zeros_on_a_grid(f_, b0 - n2)

    p1.set_ydata(b0 - n1)
    p2.set_ydata(b0 - n2)

    v1.set_xdata(f1)
    v2.set_xdata(f2)

    d1 = numpy.zeros((len(res1), 2))
    d2 = numpy.zeros((len(res2), 2))
    d1[:, 0] = res1
    d2[:, 0] = res2

    s1.set_offsets(d1)
    s2.set_offsets(d2)


sfreq.on_changed(update)
swidth.on_changed(update)

plot.show()
