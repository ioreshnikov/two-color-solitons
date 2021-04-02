import matplotlib
from matplotlib import pyplot as plot

COLOR_BLUE1 = "#0072b2"
COLOR_GREEN = "#009e73"
COLOR_VIOLET = "#9400d3"
COLOR_ORANGE = "#e69f00"
COLOR_YELLOW = "#f0e442"
COLOR_BLUE2 = "#56b4e9"
COLOR_RED = "#e51e10"
COLOR_BLACK = "#000000"
COLOR_GREY = "#7f7f7f"

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
    "font.size": 10,

    "font.monospace": [],
    "font.sans-serif": [],
    "font.serif": [],

    "legend.fontsize": 8,
    "legend.frameon": False,

    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "grid.linewidth": 0.5,
    "axes.labelsize": 10,
    "axes.titlesize": 10,

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
