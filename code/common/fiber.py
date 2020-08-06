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


cumfs = 1E6 * c / 1E15  # speed of light in μm/fs


# We define the propagation constant symbolically. This allows us to
# automatically calculate the symbolic derivatives which we then wrap
# as numpy functions.

_omega = Symbol(r"\omega")
_numer = 9.654 - 39.739 * _omega + 16.885 * _omega**2 - 2.746 * _omega**3
_denom = 1.000 -  9.496 * _omega +  4.221 * _omega**2 - 0.703 * _omega**3

_beta = _omega / cumfs * _numer / _denom

beta = lambdify(_omega, _beta, 'numpy')
beta.__doc__ = """
Propagation constant in the model fiber

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

beta1 = lambdify(_omega, _beta.diff(_omega), 'numpy')
beta.__doc__ = """
Propagation constant derivative in the model fiber

This function defines the first derivative of the propagation constant
in the model fiber as a function of frequency.

Parameters
----------

f : array_like
    an array of frequencies

Returns
-------
beta'(f) : array_like
    the propagation constant derivative at given frequencies
"""

beta2 = lambdify(_omega, _beta.diff(_omega, 2), 'numpy')
beta2.__doc__ = """
Second derviative of the propagation constant in the model fiber

This function defines the second derivative of the propagation constant
in the model fiber as a function of frequency.

Parameters
----------

f : array_like
    an array of frequencies

Returns
-------
beta''(f) : array_like
    the second derivative of the propagation constant at given frequencies
"""


def gamma(f):
    """
    Gain function spectrum.

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
