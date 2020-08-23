import matplotlib
from matplotlib import pyplot as plot


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
            "#0072b2", "#009e73", "#9400d3",
            "#e69f00", "#f0e442", "#56b4e9",
            "#e51e10", "#000000", "#7f7f7f"
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
