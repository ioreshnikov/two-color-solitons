"""
This module defines the physical parameters for the fiber used in
Melchert et al. Soliton Molecules with Two Frequencies. Phys. Rev.
Lett. 123, 243905. The original paper uses micrometers and
femtoseconds as fundamental unit of distance and time. Since we are
copying much of the constants from the paper directly, it is safer to
use the same unit system throughout the code.
"""


from scipy.constants import c
from scipy import fft
from scipy.optimize import minimize, root_scalar
from sympy import Symbol
from sympy.utilities.lambdify import lambdify
import numpy

from .helpers import filter_tails, sech, zeros_on_a_grid


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
    profile : array_like
        absorption coefficient profile
    """
    return 1j * profile * u


def gv_matching_frequencies(f1):
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
        a list of frequencies f such that β₁(f1) = β₁(f2)
    """

    d = 1E-2
    if f1 < 2:
        f = numpy.arange(f1 + d, 5.0, d)
    else:
        f = numpy.arange(0.5, f1 - d, d)
    y = beta1(f) - beta1(f1)

    f2 = []
    for f0 in zeros_on_a_grid(f, y):
        res = root_scalar(
            lambda f: beta1(f) - beta1(f1), x0=f0 - d, x1=f0 + d, xtol=1E-4)
        if res.converged:
            f2.append(res.root)

    return f2


def fundamental_soliton_amplitude(f0, t0):
    """
    Calculate fundamental soliton amplitude for the given soliton
    central frequency and the given width.

    Parameters
    ----------
    f0 : float
        soliton carrier frequency
    t0 : float, optional
        soliton width

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
    t0 : float
        soliton width
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


def fundamental_soliton_width(amp, f0):
    """
    Calculate the fundamental soliton width given the central
    frequency and the amplitude. Inverse function to
    `fundamental_soliton_amplitude`.

    Parameters
    ----------
    amp : float
        amplitude of the fundamental soliton in the model fiber
    f0 : float
        soliton carrier frequency

    Returns
    -------
    t0 : float, optional
        soliton width
    """
    return numpy.sqrt(abs(beta2(f0)) / gamma(f0)) / amp


def estimate_soliton_parameters(x, y, a1, a2, t1, t2, f1, f2, w=100):
    """
    Estimate the parameters of the soliton in the field.

    The estimation procedure is a bit involved. We construct a time-domain
    representation as a sum of two solitons each with a non-fixed
    amplitude, width and central frequency. For that trial solution we try
    to minimize the difference between the trial spectrum envelope and a
    true spectrum envelope. This miminization procedure is preformed in
    three passes:

        1. We look for all the parameters -- amplitudes, widths and
           central frequencies of all the solitons. This gives us a really
           good match, but central frequencies of the solitons have
           different group velocities.

        2. Then keeping the amplitudes and the widths of the solitons, we
           try to alter the central frequency of the first soliton, while
           taking the central frequency of the second soliton as a
           group-velocity matched one.

        3. Then keeping the frequencies and the widths of the solitons, we
           finally try to adjust the amplitude of the first soliton.

    This gives more or less nice results.

    Parameters
    ----------
    x : array_like
        coordinate grid
    y : array_like
        a soliton with the radiative tails evaluated on a grid
    a1, a2, t1, t2, f1, f2: float
        initial estimates for amplitude, width and frequency of both solitons
    w : float
        optional window width

    Returns
    -------
    a1, a2, t1, t2, f1, f2: tuple of float
        a 6-tuple of refined soliton parameters -- ampltiudes, temporal widths
        and central frequencies
    """

    y = filter_tails(x, y, w)

    f = fft.fftfreq(len(x), x[1] - x[0])
    f = 2 * numpy.pi * fft.fftshift(f)
    v = fft.fftshift(fft.ifft(y))

    fw = (f > 0.5) & (f < 4.0)
    f = f[fw]
    v = v[fw]

    def envelope_loss(a1, a2, t1, t2, f1, f2):
        y_est = (
            a1 * sech(x/t1) * numpy.exp(-1j * f1 * x) +
            a2 * sech(x/t2) * numpy.exp(-1j * f2 * x))
        v_est = fft.fftshift(fft.ifft(y_est))
        v_est = v_est[fw]

        envelope_loss = abs(abs(v) - abs(v_est)).sum()
        envelope_loss /= abs(v).sum()

        return envelope_loss

    def every_parameter_loss(args):
        return envelope_loss(*args)

    def frequency_loss(args, a1, a2, t1, t2):
        f1, *_ = args
        f2 = gv_matching_frequencies(f1)

        if not f2:
            return 1E6
        else:
            f2 = f2[-1]

        return envelope_loss(a1, a2, t1, t2, f1, f2)

    def ampltidue_loss(args, a2, t1, t2, f1, f2):
        a1, *_ = args
        return envelope_loss(a1, a2, t1, t2, f1, f2, )

    result = minimize(
        every_parameter_loss, (a1, a2, t1, t2, f1, f2))
    a1, a2, t1, t2, f1, f2 = result.x

    result = minimize(
        frequency_loss, (f1, ), args=(a1, a2, f1, f2))
    f1 = result.x
    f2 = gv_matching_frequencies(f1)[-1]

    result = minimize(
        ampltidue_loss, (a1, ), args=(a2, t1, t2, f1, f2))
    a1 = result.x

    return a1, a2, t1, t2, f1, f2
