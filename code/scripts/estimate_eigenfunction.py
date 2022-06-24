#!/usr/bin/env python3


from argparse import ArgumentParser

import matplotlib.pyplot as plot
import numpy
from scipy import fft
from scipy import interpolate

from common.fiber import beta, beta1, gamma
from common.helpers import freqs
from common.solver import lopmodes


parser = ArgumentParser()
parser.add_argument(
    "input",
    help="path to the input .npz file")
args = parser.parse_args()


npz = numpy.load(args.input)
f1 = npz["f1"]
fi = npz["fi"]

t_ = npz["t"]
f_ = freqs(t_)
u0_ = npz["u"][0, :]

nt = 2001
t = numpy.linspace(-200, 200, nt)
f = freqs(t)

reu0 = interpolate.interp1d(t_, u0_.real)
imu0 = interpolate.interp1d(t_, u0_.imag)
u0 = reu0(t) + 1j * imu0(t)

v0 = fft.ifft(u0, norm="ortho")

v1 = v0.copy()
v1[f > 2.0] = 0.0

v2 = v0.copy()
v2[f < 2.0] = 0.0

u1 = fft.fft(v1, norm="ortho")
u2 = fft.fft(v2, norm="ortho")

import matplotlib.pyplot as plot
plot.plot(t, abs(u1)**2 + abs(u2)**2)
plot.plot(t, abs(u0)**2)
plot.show()

# plot.subplot(2, 1, 1)
# plot.plot(t, abs(u0))
# plot.xlim(t.min(), t.max())
# plot.xlabel(r"$t$, fs")
# plot.ylabel(r"$|u_0|$, a.u.")

# plot.subplot(2, 1, 2)
# plot.plot(f, abs(v0))
# plot.xlim(0.5, 4.0)
# plot.xlabel(r"$\omega$, rad/fs")
# plot.ylabel(r"$|v_0|$, a.u")

# plot.show()


# Construct a frequency filtered and group velocity compensated
# dispersive profile.
def filtered_beta(f):
    b = beta(f) - beta(f1) - beta1(f1) * (f - f1)
    # NOTE: We don't suppress negative frequencies here. Otherwise the solver
    # can put ANY trash in the negative frequency range and it still will be a
    # valid solution.
    # b[(f <= 0)] = 0
    return b


# Frequency filter the nonlinear coefficient.
def filtered_gamma(f):
    g = gamma(f)
    g[(f <= 0)] = 0
    # NOTE: Here we have some problematic division, so we'd rather filter the
    # negative frequency region.
    return g

ki = filtered_beta(fi)
states = lopmodes(
    t, filtered_beta, filtered_gamma,
    2 * abs(u0)**2, u0**2, ki, nmodes=32)


rhs = filtered_gamma(f) * fft.ifft(
    2 * (abs(u1)**2 + abs(u2)**2) * numpy.exp(-1j * fi * t),
    norm="ortho")


spec = 0
for k, pk in states:
    coef = sum(rhs * pk.conjugate() + rhs.conjugate() * pk)
    norm = sum(2 * pk * pk.conjugate())
    spec += coef / norm * pk


import matplotlib.pyplot as plot
plot.plot(f, abs(spec))
plot.xlim(0.5, 4.0)
plot.show()
