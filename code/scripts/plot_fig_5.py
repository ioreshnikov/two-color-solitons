#!/usr/bin/env python3


from argparse import ArgumentParser

from matplotlib import pyplot as plot
from matplotlib import colors
from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator

import numpy

from common.helpers import freqs
from common.plotter import pr_setup, pr_publish


parser = ArgumentParser()
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
u = npz["u"]
v = npz["v"]

f = freqs(t, shift=True)


# Now that we have estimated everything, we don't need the full output
# at the original resolution, so we can start windowing and
# subsampling the matrices.
tw = (t > -1000) & (t < +2000)
t = t[tw]
u = u[:, tw]

ssx = len(t) // 1000 or 1
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

vn = abs(v0).max()
v0 = abs(v0) / vn
v1 = abs(v1) / vn

del v

# To makes the plots more readable, we rescale the time to picoseconds
# and distance to centimeters.
t = t / 1000
z = z / 10000


# Load soliton parameters
npz = numpy.load(args.parameters)
a1s = npz["a1s"]
a2s = npz["a2s"]
f1s = npz["f1s"]
f2s = npz["f2s"]


pr_setup()

height_ratios = [2/3, 3/4, 2/3]
plot.figure(
    figsize=(3.2, sum(1.4 * h for h in height_ratios)))
gs = gridspec.GridSpec(3, 2, height_ratios=height_ratios)


plot.subplot(gs[0, :])
plot.pcolormesh(
    z, t, u.T**2,
    cmap="jet",
    norm=colors.LogNorm(vmin=args.vmin),
    rasterized=True,
    shading="auto")
plot.xlim(z.min(), z.max())
plot.ylim(-0.25, +1.75)
plot.yticks([0.0, +0.5, +1.0, +1.5])
plot.xlabel(r"Distance $z$, cm")
plot.ylabel(r"Delay $t$, ps")
plot.colorbar()

plot.annotate(
    "(a)",
    xy=(1.0, 0.0),
    xytext=(-16.0, +6.0),
    xycoords="axes fraction",
    textcoords="offset points",
    color="white")


# Second and third panel: soliton parameters
plot.subplot(gs[1, 0])

iz = len(z) // 10

da1 = a1s - a1s[iz]
da2 = a2s - a2s[iz]
plot.plot(z, da1, color="black", label="1")
plot.plot(z, da2, color="gray",  label="2")
plot.legend(ncol=2, loc="upper center")

plot.xlim(z.min(), z.max())

ymin = min(da1.min(), da2.min())
ymax = max(da1.max(), da2.max())
margin = 0.25 * (ymax - ymin)
plot.ylim(-0.01, +0.01)
plot.xlabel(r"Distance $z$, cm")
plot.ylabel(r"$\Delta A_n$, a.u")

plot.annotate(
    "(b.1)",
    xy=(0.0, 0.0),
    xytext=(+6.0, +6.0),
    xycoords="axes fraction",
    textcoords="offset points")

plot.subplot(gs[1, 1])

df1 = f1s - f1s[iz]
df2 = f2s - f2s[iz]

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


plot.tight_layout(pad=0.25)
pr_publish(args.output)
