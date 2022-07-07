"""
Physical parameters of the fiber.

We are using a reference model from Melchert et al. Soliton Molecules with Two
Frequencies. Phys. Rev. Lett. 123, 243905. The original paper uses micrometers
and femtoseconds as fundamental unit of distance and time. Since we are with
the same physical setting, it is safer to adopt the same unit system
throughout the code.
"""


from scipy.constants import c
from sympy import Symbol
from sympy.utilities.lambdify import lambdify


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
    Propagation constant in the model fiber.

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
    Frequency-dependent Kerr nonlinearity coefficient.

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


def kerr_op(z, t, u):  # noqa: D205, D400
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
