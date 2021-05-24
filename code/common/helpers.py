from scipy import fft
from scipy import signal
import numpy


def sech(x):
    """
    Calculates hyperbolic secant of the input without overflowing.

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
# original implementation by O. Melchert in his code Ocean Code
# capsule "Code and data for FMAS_MS" that can be found here
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
    zw[1:int(zw.size/2)] = 2 * zw[1:int(zw.size/2)]
    zw[int(zw.size/2)+1:] = 0j
    return fft.fft(zw)


def zeros_on_a_grid(x, y):
    """
    Given a function evaluated on a grid, find points where the
    function changes sign.

    Parameters
    ----------
    x : array_like
        coordinate grid
    y : array_like
        function evaluated on a grid
    """
    yp = y[:-1]
    yn = y[1:]

    mask = (yp * yn) < 0

    xp = x[:-1][mask]
    xn = x[1:][mask]

    return (xp + xn) / 2


def filter_tails(x, y, w=100):
    """
    Filter the radiative tails from the soliton using a super-Gaussian
    window.

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
