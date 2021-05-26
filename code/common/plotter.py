import matplotlib
from matplotlib import colors
from matplotlib import pyplot as plot
import numpy
from numpy import fft


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


def gpc_setup(z, t, vmin=1E-6, title=None):
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
    dt = t[1] - t[0]

    f = 2 * numpy.pi * fft.fftshift(fft.fftfreq(nt, dt))

    # The images we get are excessively large. We are going to
    # sub-sample everything for faster rendering.
    sst = int(nt / 1024)
    nt = int(nt / sst)

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
    plot.subplot(2, 1, 1)
    mesh1 = plot.pcolormesh(
        z, t[::sst], _.T,
        cmap="jet",
        norm=colors.LogNorm(vmin=vmin, vmax=1))
    plot.xlabel(r"Distance $z$, cm")
    plot.ylabel(r"Delay $t$, ps")

    plot.subplot(2, 1, 2)
    mesh2 = plot.pcolormesh(
        z, f[::sst], _.T,
        cmap="jet",
        norm=colors.LogNorm(vmin=vmin, vmax=1))
    plot.ylim(0.0, 4.0)
    plot.xlabel(r"Distance $z$, cm")
    plot.ylabel(r"Frequency $\omega$, rad/fs")

    plot.ion()
    plot.show(block=False)
    plot.pause(1E-3)

    # Construct the callback and return it
    def gpc(z, t, f, u, v):
        u = u[:, ::sst]
        v = v[:, ::sst]

        # Scale to [0, 1] interval
        u = abs(u) / abs(u).max()
        v = abs(v) / abs(v).max()

        # Rotate it appropriately and then switch to instantaneous
        # intensity
        u = u.T**2
        v = v.T**2

        mesh1.set_array(u[:-1, :-1].ravel())
        mesh2.set_array(v[:-1, :-1].ravel())
        #                 ^^^ yes, super weird :)

        plot.pause(1E-3)

    return gpc
