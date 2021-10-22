import logging

import matplotlib
from matplotlib import colors
from matplotlib import pyplot as plot
import numpy

from .helpers import freqs


logging.getLogger("matplotlib").setLevel(logging.WARNING)


COLOR_BLUE1 = "#0072b2"
COLOR_GREEN = "#009e73"
COLOR_VIOLET = "#9400d3"
COLOR_ORANGE = "#e69f00"
COLOR_YELLOW = "#f0e442"
COLOR_BLUE2 = "#56b4e9"
COLOR_RED = "#e51e10"
COLOR_BLACK = "#000000"
COLOR_GREY = "#7f7f7f"

XSMALL_FONT_SIZE = 6
SMALL_FONT_SIZE = 7
REGULAR_FONT_SIZE = 8

PGF_PARAMS = {
    "figure.figsize": (3.1, 2.5),

    "text.usetex": True,
    "pgf.texsystem": "pdflatex",
    "pgf.preamble": r"""
        \usepackage[utf8]{inputenc}
        \usepackage[T2A]{fontenc}
        \renewcommand{\familydefault}{\sfdefault}
    """,
    "axes.unicode_minus": False,

    "font.family": "serif",
    "font.size": REGULAR_FONT_SIZE,

    "font.monospace": [],
    "font.sans-serif": [],
    "font.serif": [],

    "legend.fontsize": SMALL_FONT_SIZE,
    "legend.frameon": False,

    "xtick.labelsize": XSMALL_FONT_SIZE,
    "ytick.labelsize": XSMALL_FONT_SIZE,
    "grid.linewidth": 0.5,
    "axes.labelsize": REGULAR_FONT_SIZE,
    "axes.titlesize": REGULAR_FONT_SIZE,

    "lines.linewidth": 0.8,
    "axes.linewidth": 0.5,

    "axes.prop_cycle": plot.cycler(
        "color", [
            COLOR_BLUE1, COLOR_GREEN, COLOR_VIOLET,
            COLOR_ORANGE, COLOR_YELLOW, COLOR_BLUE2,
            COLOR_RED, COLOR_BLACK, COLOR_GREY
        ]),
}


def pr_setup():
    """
    Setup publication-ready plotting frontend the way I like it.
    """
    plot.switch_backend("pgf")
    matplotlib.rcParams.update(PGF_PARAMS)


def pr_publish(filename, dpi=None):
    """
    Save the plot as a publication-ready image.

    Parameters
    ----------
    filename : path
        path to the output document
    dpi : int
        DPI of the output image. If nothing is passed, 400 is used for
        pdf documents and 150 otherwise.
    """
    if dpi is None:
        pdf = filename.endswith(".pdf")
        dpi = 400 if pdf else 150

    plot.tight_layout(pad=0.5)
    plot.gca().tick_params(which="both", direction="out")
    plot.savefig(filename, dpi=dpi)


def gpc_setup(
    z, t,
    tmin=None, tmax=None, fmin=None, fmax=None,
    vmin=1E-6, title=None, cols=1):
    """
    Setup gridpoint callback to be passed to the solver.

    This function creates an interactive plotter window and returns a
    callback that updates the time- and frequency domain subplots in
    there interactively as we get new data from the solver.

    Returns
    -------
    gpc : callable with the signature of gpc(z, t, f, u, v)
        the plot updating function. This can be passed directly to
        the solver.
    """

    # Calculate the frequency grid. Usually not computed by the caller
    # at this point, but will come from the solver.
    nz = len(z)
    nt = len(t)

    f = freqs(t, shift=True)

    if fmin is None:
        fmin = f.min()

    if fmax is None:
        fmax = f.max()

    if tmin is None:
        tmin = t.min()

    if tmax is None:
        tmax = t.max()

    # The images we get are excessively large. We are going to
    # sub-sample everything for faster rendering.
    if nt > 1024:
        sst = nt // 1024
    else:
        sst = 1
    nt = nt // sst

    # Prepare the stand-in for the future data matrices
    _ = numpy.empty((nz, nt))
    _[:] = numpy.nan

    plot.figure(figsize=(8, 8))
    if title:
        plot.suptitle(title)

    # Scale the axis to the approriate units
    z = z / 1E4
    t = t / 1E3

    # Prepare the panels with the stand-in matrix
    meshes = []
    for col in range(cols):
        plot.subplot(2, cols, 2 * col + 1)
        mesh1 = plot.pcolormesh(
            z, t[::sst], _.T,
            cmap="jet",
            norm=colors.LogNorm(vmin=vmin, vmax=1))
        plot.ylim(tmin, tmax)
        plot.xlabel(r"Distance $z$, cm")
        plot.ylabel(r"Delay $t$, ps")

        plot.subplot(2, cols, 2 * col + 2)
        mesh2 = plot.pcolormesh(
            z, f[::sst], _.T,
            cmap="jet",
            norm=colors.LogNorm(vmin=vmin, vmax=1))
        plot.ylim(fmin, fmax)
        plot.xlabel(r"Distance $z$, cm")
        plot.ylabel(r"Frequency $\omega$, rad/fs")

        meshes.append([mesh1, mesh2])

    plot.ion()
    plot.show(block=False)
    plot.pause(1E-3)

    # Construct the callback and return it
    def gpc(z, t, f, u, v):
        for col in range(cols):
            if cols == 1:
                u_ = u[:, ::sst]
                v_ = v[:, ::sst]
            else:
                u_ = u[:, ::sst, col]
                v_ = v[:, ::sst, col]

            # Scale to [0, 1] interval
            u_ = abs(u_) / abs(u).max()
            v_ = abs(v_) / abs(v).max()

            # Rotate it appropriately and then switch to instantaneous
            # intensity
            u_ = u_.T**2
            v_ = v_.T**2

            meshes[col][0].set_array(u_[:-1, :-1].ravel())
            meshes[col][1].set_array(v_[:-1, :-1].ravel())
            #                          ^^^ yes, super weird :)

        plot.pause(1E-3)

    return gpc

