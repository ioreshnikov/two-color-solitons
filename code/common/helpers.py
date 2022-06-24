"""Auxiliary functions."""


from scipy import fft
from scipy import signal
from scipy.optimize import minimize, root_scalar

import numpy

from .fiber import beta, beta1, beta2, gamma


def sech(x):
    """
    Calculate hyperbolic secant of the input without overflowing.

    Parameters
    ----------
    x : array_like
        input array

    Returns
    -------
    sech(x) : array_like
        an element-wise hyperbolic secant
    """
    e = numpy.exp(-abs(x))
    return 2*e / (1 + e**2)


# The function below is a copy (with alterations to formatting) of the
# original implementation by O. Melchert in his Ocean Code capsule
# "Code and data for FMAS_MS" that can be found here
#     https://codeocean.com/capsule/0624342


def to_analytic(ut):
    """
    Compute analytic signal.

    Parameters
    ----------
    ut : array_like
        time-domain representation of signal

    Returns
    -------
    zt : array_like
        time-domain representation of analytic signal

    Note
    ----
    - Computes standard discrete-time 'analytic' signal following Eq. (8)
      in sect. IV of Ref. [1] via the two-step procedure:

      1. Compute N-point DFT `uw[:]` using FFT of original `N = ut.size`
         samples: `uw = fft(ut)`
      2. Construct N-point one-sided discrete-time 'analytic' singal trafo:

             zw[0] = uw[0] zw[N/2] = uw[N/2]
             zw[1:N/2] = 2*uw[1:N/2] zw[N/2+1:] = 0j

    - The above procedure ensures the following two properties:
      1. Real part of reverse transform yields the original (real) signal:

                ut = real(ifft(zw))

      2. Real and imaginary component are orthogonal over interval:

                ut = real(ifft(zw))
                utIm = imag(ifft(zw))
                np.sum(ut*utIm) < 1e-8

    References
    ----------
    [1] Computing the Discrete-Time 'Analytic' Signal via FFT
        Marple, S. L. Jr. IEEE Trans. Sig. Proc., 47 (1999) 2600
    """
    zw = fft.ifft(ut)
    zw[1:(zw.size//2)] = 2 * zw[1:(zw.size//2)]
    zw[(zw.size//2)+1:] = 0j
    return fft.fft(zw)


def freqs(x, shift=False):
    """
    Construct corresponding frequency vector.

    Parameters
    ----------
    x : array_like
        coordinate grid
    shift : bool
        if True then perform fftshift on a result array before returning

    Returns
    -------
    k : array_like
        a vector of spatial frequencies
    """
    nx = len(x)
    dx = x[1] - x[0]
    k = 2 * numpy.pi * fft.fftfreq(nx, dx)
    if shift:
        k = fft.fftshift(k)
    return k


def zeros_on_a_grid(x, y, mode="midpoint"):
    """
    Find grid points where a function changes sign.

    Parameters
    ----------
    x : array_like
        coordinate grid
    y : array_like
        function evaluated on a grid
    mode : string, either "midpoint" or "index"
        in midpoint mode one gets the location of a zero approximated as
        a midpoint between the grid points where the sign change occurs;
        in index mode one gets the indices of corresponding points in `x`
    """

    assert mode in ("index", "midpoint"), \
        'only "index" and "midpoint" are allowed for `mode`'

    yp = y[:-1]
    yn = y[1:]

    mask = (yp * yn) < 0

    if mode == "index":
        return numpy.append([False], mask)

    xp = x[:-1][mask]
    xn = x[1:][mask]

    return (xp + xn) / 2


def filter_tails(x, y, w=100):
    """
    Filter the radiative tails using a super-Gaussian window.

    Parameters
    ----------
    x : array_like
        coordinate grid
    y : array_like
        a soliton with the radiative tails evaluated on a grid
    w : float
        optional window width

    Returns
    -------
    y' : array_like
        a filtered solution
    """
    # We construct a super-Gaussian window around the soliton peak.
    # Multiplying the soliton by that window we suppress the radiation
    # tails.
    idx_max = abs(y).argmax()
    x0 = x[idx_max]
    window = numpy.exp(- ((x - x0)/w)**8)
    return window * y


def align_max(x, y):
    """
    Align maximum of signal y with x=0.

    Parameters
    ----------
    x : array_like
        coordinate grid
    y : array_like
        function to be aligned

    Returns
    -------
    y' : array_like
       the same function shifted so that max(y) is at x=0.
    idx_diff : int
       shift distance
    """
    idx_max = abs(y).argmax()
    idx_zero = abs(x).argmin()
    idx_diff = idx_zero - idx_max
    return numpy.roll(y, idx_zero - idx_max), idx_diff


def peaks_close_to(x, y, x0s):
    """
    Find peaks in a vicinity of given points.

    Parameters
    ----------
    x : array_like
        coordinate grid
    y : array_like
        the function
    x0 : array_like
        points of interest

    Note
    ----
    Algorithm parameters (vicinity size, prominence thresholds) work
    for the figures we use in the paper, but otherwise they're not
    tested. It is in no way generic.
    """
    result = []

    for x0 in x0s:
        win = abs(x - x0) < 0.05

        xw = x[win]
        yw = y[win]

        peaks, props = signal.find_peaks(
            yw / yw.max(),
            prominence=(0.5, None))

        if not len(peaks):
            continue

        idx, *_ = peaks
        result.append((xw[idx], yw[idx]))

    return result


def gv_matching_frequencies(f1):
    """
    Given a frequency, find another ones with the same group velocity.

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
            lambda f: beta1(f) - beta1(f1), x0=f0 - d, x1=f0 + d, xtol=1E-3)
        if res.converged:
            f2.append(res.root)

    return f2


def fundamental_soliton_amplitude(f0, t0):
    """
    Fundamental soliton amplitude for a given frequency and width.

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
    Dispersive relation for a fundamental soliton at a given frequency.

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
    Fundametal soliton width at a given frequency.

    Inverse function to `fundamental_soliton_amplitude`.

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


def frame_of_reference(f, f0):
    """
    Dispersive relation for a plane wave at a given frequency.

    Parameters
    ----------
    f : array_like
        frequency grid
    f0 : float
        carrier frequency of the wave

    Note
    ----
    This is useful when we want to calculate dispersive relationship of a
    wave in a frame of reference co-moving with the soliton.
    """
    return beta(f0) + beta1(f0) * (f - f0)


def _envelope_ansatz(x, f, a1, a2, t1, t2, f1, f2):
    nx = len(x)
    dx = x[1] - x[0]
    sf = 1 / (nx * dx)

    return sf * numpy.pi * (
        a1 * t1 * sech(numpy.pi/2 * t1 * (f - f1)) +
        a2 * t2 * sech(numpy.pi/2 * t2 * (f - f2)))


def estimate_soliton_parameters(x, y, a1, a2, t1, t2, f1, f2, w=100, match_f=True):  # noqa: E501
    """
    Estimate the parameters of the soliton in the field.

    The estimation procedure is a bit involved. We construct a time-domain
    representation as a sum of two solitons each with a non-fixed
    amplitude, width and central frequency. For that trial solution we try
    to minimize the difference between the trial spectrum envelope and a
    true spectrum envelope. This miminization procedure is preformed in either
    one or three passes:

        1. We look for all the parameters -- amplitudes, widths and
           central frequencies of all the solitons. This gives us a really
           good match, but central frequencies of the solitons have
           different group velocities.

           If `match_f` is set to False, we return the estimates. Otherwise:

        2. Keeping the amplitudes and the widths of the solitons, we try
           to alter the central frequency of the first soliton, while taking
           the central frequency of the second soliton as a group-velocity
           matched one.

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

    f = freqs(x, shift=True)
    v = fft.fftshift(fft.ifft(y))

    fw = (f > 0.5) & (f < 4.0)
    f = f[fw]
    v = v[fw]

    # Subsample spectrum and frequency if necessary. Around 500 points
    # should be enough for our estimations.
    ssf = len(f) // 500
    f = f[::ssf]
    v = v[::ssf]

    def envelope_loss(a1, a2, t1, t2, f1, f2):
        v_est = _envelope_ansatz(x, f, a1, a2, t1, t2, f1, f2)

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

    if not match_f:
        return a1, a2, t1, t2, f1, f2

    result = minimize(
        frequency_loss, (f1, ), args=(a1, a2, f1, f2))
    f1, *_ = result.x
    f2 = gv_matching_frequencies(f1)[-1]

    result = minimize(
        ampltidue_loss, (a1, ), args=(a2, t1, t2, f1, f2))
    a1, *_ = result.x

    return a1, a2, t1, t2, f1, f2
