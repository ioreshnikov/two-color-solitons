"""
This module defines the physical parameters for the fiber used in
Melchert et al. Soliton Molecules with Two Frequencies. Phys. Rev.
Lett. 123, 243905. The original paper uses micrometers and
femtoseconds as fundamental unit of distance and time. Since we are
copying much of the constants from the paper directly, it is safer to
use the same unit system throughout the code.
"""


from scipy.constants import c
from sympy import Symbol
from sympy.utilities.lambdify import lambdify
import numpy


cumfs = 1E6 * c / 1E15  # speed of light in μm/fs


# We define the propagation constant symbolically. This allows us to
# automatically calculate the symbolic derivatives which we then wrap
# as numpy functions.
_omega = Symbol(r"\omega")
_numer = 9.654 - 39.739 * _omega + 16.885 * _omega**2 - 2.746 * _omega**3
_denom = 1.000 -  9.496 * _omega +  4.221 * _omega**2 - 0.703 * _omega**3  # noqa

_beta = _omega / cumfs * _numer / _denom


def beta(f):
    """
    This function defines the propagation constant in the model fiber
    as a function of frequency.

    Parameters
    ----------
    f : array_like
        an array of frequencies

    Returns
    -------
    beta(f) : array_like
        the propagation constant at given frequencies
    """
    return lambdify(_omega, _beta, 'numpy')(f)


def beta1(f):
    """
    First derivative of the propagation constant in the model fiber.

    Parameters
    ----------
    f : array_like
        an array of frequencies

    Returns
    -------
    beta'(f) : array_like
        the propagation constant derivative at given frequencies
    """
    return lambdify(_omega, _beta.diff(_omega), 'numpy')(f)


def beta2(f):
    """
    Second derviative of the propagation constant in the model fiber.

    Parameters
    ----------

    f : array_like
        an array of frequencies

    Returns
    -------
    beta''(f) : array_like
        the second derivative of the propagation constant at given frequencies
    """
    return lambdify(_omega, _beta.diff(_omega, 2), 'numpy')(f)


def gamma(f):
    """
    This function defines the frequency-dependent Kerr nonlinearity
    coefficient.

    Parameters
    ----------
    f : array_like
        an array of frequencies

    Returns
    -------
    gamma(f) : array_like
        the coefficient for given frequencies
    """
    return 3/8 * f**2 / cumfs**2 / beta(f)


def kerr_op(z, t, u):
    """
    Kerr nonlinearity operator with a signature suitable for GLNSE
    solver.

    Parameters
    ----------
    z : float, ignored
        current fiber position, has no effect
    t : array_like, ignored
        coordinate grid, has no effect
    u : array_like
        instantaneous field in time-domain representation
    """
    return abs(u)**2 * u


def linear_absorption_op(z, t, u, profile):
    """
    Linear absorption operator with the given absorption profile.

    Parameters
    ----------
    z : float, ignored
        current fiber position, has no effect
    t : array_like, ignored
        coordinate grid, has no effect
    u : array_like
        instantaneous field in time-domain representation
    absorber : array_like
        absorption coefficient profile
    """
    return 1j * profile * u


def gv_matching_frequencies(f1, beta1=beta1):
    """
    Given a frequency find another ones where the group velocity is the
    same.

    Parameters
    ----------
    f1 : float
        the first frequency

    Returns
    -------
    f2 : list(float)
        a list of frequencies f such that β₁(f1) = β₁(f2)`

    Note
    ----
    Root finding is very crude -- we evaluate beta1() on a very fine grid
    in the interval [f1 + ε, 5] or [0.5, f1 - ε] (depending on f1) and
    then return the frequencies where the difference β₁(f) - β₁(f1)
    changes sign.
    """

    d = 1E-6
    if f1 < 2:
        f = numpy.arange(f1 + d, 5.0, d)
    else:
        f = numpy.arange(0.5, f1 - d, d)
    y = beta1(f) - beta1(f1)

    yp = y[:-1]
    yn = y[1:]

    mask = (yp * yn) < 0

    fp = f[:-1][mask]
    fn = f[1:][mask]

    return (fp + fn) / 2


def fundamental_soliton_amplitude(f0, t0=1.0):
    """
    Calculate fundamental soliton amplitude for the given soliton
    central frequency and the given width.

    Parameters
    ----------
    f0 : float
        soliton carrier frequency

    t0 : float, optional
        soliton width, 1.0 by default

    Returns
    -------
    amp : float
        amplitude of the fundamental soliton in the model fiber
    """
    return numpy.sqrt(abs(beta2(f0)) / gamma(f0)) / t0


def fundamental_soliton_dispersive_relation(f0, t0, f):
    """
    Dispersive relation for a fundamental soliton at given frequency.

    Parameters
    ----------
    f0 : float
        soliton carrier frequency
    f : array_like
        an array of frequencies

    Returns
    -------
    k : array_like
        dispersive relation for the soliton
    """
    a = fundamental_soliton_amplitude(f0, t0)
    g = gamma(f0)
    return g * a**2 / 2 + beta(f0) + beta1(f0) * (f - f0)
