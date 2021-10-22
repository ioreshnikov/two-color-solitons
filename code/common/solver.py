"""
Solvers.

Here we keep the numeric solvers needed for the paper, namely a GNLSE solver
and a solver for the linearization operator spectral problem.
"""


from datetime import timedelta
import logging
import time

from humanize import precisedelta
from scipy import fft
from scipy import integrate
from scipy import linalg
from scipy import sparse
import numpy

from .helpers import freqs


def _report_progress(n, nt, timer_start):
    """Report integration progress to the log."""
    progress = (n / nt)
    timer_current = time.time()
    timer_elapsed = timer_current - timer_start
    timer_left = timer_elapsed * (1 - progress) / progress

    # Format timers to human-readable strings (hence 'hr_').
    hr_elapsed = precisedelta(
        timedelta(seconds=timer_elapsed),
        minimum_unit="seconds",
        format="%02d")
    hr_left = precisedelta(
        timedelta(seconds=timer_left),
        minimum_unit="seconds",
        format="%02d")

    logging.info(
        "  %5.2f%%  %s elapsed  %s left",
        100 * progress, hr_elapsed, hr_left)


class IntegrationResult:
    """
    This is a container class for the integration results.

    Parameters
    ----------
    t : array_like
        an array of time points where the integrated function is
        evaluated
    x : array_like
        the coordinate grid
    k : array_like
        the frequency grid
    u : array_like of shape len(t)×len(x)
        integrated field in time domain evaluated on a grid
    v : array_like of shape len(t)×len(x)
        integrated field in frequency domain evaluated on a grid
    successful : bool
        a flag indicating that the integration was successful
    error_code : int
        return code of the solver if the integration failed

    Note
    ----
    `u` and `v` are output matrices of the integrated field in
    coordinate-domain and spectral-domain representations. Row number
    is position in the time array, column number is the position in
    coordinate/frequency array.
    """

    def __init__(self, t, x, k, u, v, error_code=None):  # noqa: D107
        self.t = t
        self.x = x
        self.k = fft.fftshift(k)
        self.u = u
        self.v = v
        self.error_code = error_code

    @property
    def successful(self):  # noqa: D400
        """Was this run succesful?"""
        return self.error_code is None


def gnlse(t, x, u0, beta, gamma, nonlin, lin=None, dt=None, gpc=None):
    """
    Integrate a GNLSE using the integrating factor method.

    This function integrates a generalized version of nonlinear Schrödinger
    equation

        ∂ₜ ũ = i β(k) ũ(t, k)
             + i γ(k) F{ N(t, x, u(t, x)) }
             + i F{ L(t, x, u(t, x)) },

    where ũ(t, k) is the spectrum of the unknown field, β(k) is a
    dispersive operator, and γ(k) is a gain coefficient. u(t, x) is the
    coordinate-domain representation of the same unknown field, and
    nonolinear operator N(t, x, u(t, x)) is defined in terms of that
    coordinate representation. L(t, x, u(t, x)) is an auxiliary linear
    operator that can be used to introduce an absorbing boundary layer.
    It does not have to have a physical meaning.

    The integration is performed using the integrating factor method as
    described by J.M. Dudley & J.R. Taylor in Chapter 3 of Supercontinuum
    Generation in Optical Fibers, CUP 2010. Instead of integrating the
    original equation we resort to integrating a modified version

        ∂ₜ v = i γ(k) F{ N(t, x, u(t, x)) }
             + i F{ L(t, x, u(t, x)) },

    where v = v(t, k) is the modified spectrum that is defined as

        ũ(t, k) = exp(i β(k) t) v(t, k).

    The modifed equation is supposed to be non-stiff, which allows us to
    apply almost any third-party solver. We chose a scipy-provided wrapper
    of ZVODE solver from ODEPACK, which offers error control and an
    adaptive scheme for step-size selection.

    Parameters
    ----------
    t : array_like
        an array of time points where the integrated function should
        be evaluated
    x : array_like
        the coordinate grid
    u0 : array_like
        the initial condition at t[0]
    beta : callable with the signature of beta(k)
        the dispersive profile as a function of frequency
    gamma : callable with the signature of gamma(k)
        frequency-dependent gain part of the nonlinear operator
    nonlin : callable with the signature of nonlin(t, x, u)
        time-domain part of the nonlinear operator
    lin : callable with the signature of lin(t, x, u)
        time-domain linear operator
    dt : float, optional
        an optional internal time step used for integration
    gpc : callable with the signature of gcp(t, x, k, u, v)
        an optional callback executed after computing up to the next
        grid point

    Returns
    -------
    result : an instance of IntegrationResult
    """
    # Pre-allocate the output matrices.
    nt = len(t)
    nx = len(x)

    u = numpy.zeros((nt, nx), dtype=complex)
    v = numpy.zeros((nt, nx), dtype=complex)

    # Put the initial conditions in time and frequency domains into
    # the output matrices
    v0 = fft.ifft(u0, norm="ortho")
    u[0, :] = u0
    v[0, :] = fft.fftshift(v0)

    # Prepare the frequency scale and evaluate beta on the scale
    k = freqs(x)
    B = beta(k)
    G = gamma(k)

    # Prepare the RHS we feed to the solver
    def rhs(t_, v_):
        # Scale the spectrum by the accumulated phase shift due to
        # the dispersion and transform to coordinate space
        exp = numpy.exp(1j * B * (t_ - t[0]))
        u_ = fft.fft(exp * v_, norm="ortho")

        # Apply nonlinear operator N() and, maybe, linear operator L()
        # as well, transform back to the modified spectrum and return.
        ret = G * fft.ifft(nonlin(t, x, u_), norm="ortho")
        if lin:
            ret += fft.ifft(lin(t, x, u_), norm="ortho")

        return 1j / exp * ret

    # Configure the solver. We pick ZVODE from ODEPACK which provides
    # error control and adaptive stepping out of the box. Since we are
    # working on a modified spectrum, we can hope the problem to be
    # non-stiff, but it's a good idea *not* to impose the choice of
    # the integration method and error tolerances -- this will be
    # estimated by the solver itself.
    ode = integrate.ode(rhs)
    ode.set_integrator("zvode", rtol=1E-6)

    # Those are the internal loop variables. `t_` holds the current
    # integration time, `v_` is the current modified spectrum.
    t_ = t[0]
    v_ = v0

    timer_start = time.time()
    logging.info("Integrating:")

    for n in range(1, nt):
        _report_progress(n, nt, timer_start)

        # Pick a step size. THIS IS NOT THE ADAPTIVE STEP-SIZE CONTROL
        # --- it is done by the solver. But sometimes the time step
        # given by the time grid is too coarse for the given precision
        # and the solver exhaust maximum number of integration steps
        # while integrating the equation over the grid step. Providing
        # a finer time grid will result in larger output matrices, and
        # for long simulations those might not fit in the memory
        # anymore. Therefore we provide an additional work-around:
        # user can specify an intermediate time step that is used
        # internally by the loop but not saved anywhere.
        dt_ = dt or t[n] - t_

        while t_ < t[n]:
            # Set the initial value for the current integration step.
            ode.set_initial_value(v_, t_)

            # Pick the next step -- either dt_ or until the next grid
            # point.
            if t_ + dt_ > t[n]:
                t_ = t[n]
            else:
                t_ += dt_

            # Integrate for one step and check the solver status.
            v_ = ode.integrate(t_)
            if not ode.successful():
                return IntegrationResult(
                    t, x, k, u, v, error_code=ode.get_return_code())

        # Calculate the proper spectrum and then save the spectrum and
        # the coordinate representation into the output matrices.
        exp = numpy.exp(1j * B * (t_ - t[0]))
        u[n, :] = fft.fft(exp * v_, norm="ortho")
        v[n, :] = fft.fftshift(exp * v_)

        if gpc:
            gpc(t, x, k, u, v)

    logging.info("Done!")
    return IntegrationResult(t, x, k, u, v)


def _fftm(n):
    """
    Construct an FFT matrix.

    Parameters
    ----------
    n : int
        matrix size

    Returns
    -------
    FFT : a n×n matrix
        a matrix that acts as fft when multiplied by a vector
    """
    return fft.fft(numpy.eye(n), norm="ortho")


def _ifftm(n):
    """
    Construct an IFFT matrix.

    Note
    ----
    This is exactly FFT⁻¹, but it's faster to construct it directly instead
    of inverting a matrix.

    Parameters
    ----------
    n : int
        matrix size

    Returns
    -------
    IFFT : a n×n matrix
        a matrix that acts as ifft when multiplied by a vector
    """
    return fft.ifft(numpy.eye(n), norm="ortho")


def _flip(n):
    """
    Operator that performs a frequency mirroring.

    Parameters
    ----------
    n : int
        matrix size

    Returns
    -------
    FLIP : a n×n matrix
        a matrix that swaps frequencies, i.e. FLIP f(k) = f(-k)
    """
    assert n % 2 == 1, "only odd-sized vectors are flippable"

    FLIP = numpy.zeros((n, n))
    FLIP[0,  0] = 1
    FLIP[1:, 1:] = numpy.fliplr(numpy.eye(n - 1))

    return FLIP


def _lopm(x, beta, gamma, v1, v2):
    """
    Construct matrix representation of the linearization operator.

    This function constructs a matrix representation for operator L below

        L ψ(k) = β(k) ψ(k) + γ(k) F{ V₁ ψ(x) + V₂ ψ⁺(x) }

    where ψ(x) is an eigenfunction, ψ⁺(x) is the complex conjugate, β(k) is
    the dispersiion operator, γ(k) is an arbitrary function of frequency and
    V₁(x) and V₂(x) are two arbitrary potentials.

    Note
    ----
    ψ(x) is naturally calculated here as IFFT @ ψ(k) and for ψ⁺(x) is
    IFFT @ ψ⁺(-k).

    Parameters
    ----------
    x : array_like
       coordinate grid
    beta : callable with a signature of beta(k)
        the dispersive profile as a function of frequency
    gamma : callable with a signature of gamma(k)
        arbitrary potential scaling as a function of frequency
    v1, v2 : array_like
        potentials as functions of coordinate

    Returns
    -------
    L : nx × nx matrix
        a matrix representation of the linearization operator
    """

    # Construct the frequency grid
    k = freqs(x)

    # Operators β(k) and γ(k) are simple diagonal matrices.
    B = numpy.diag(beta(k))
    G = numpy.diag(gamma(k))

    # The potentials are wrapped into the diagonal matrices
    V1 = numpy.diag(v1)
    V2 = numpy.diag(v2)

    # Finally, we construct the auxiliary fft matrices
    FFT  = _fftm(len(x))  # noqa: E221
    IFFT = _ifftm(len(x))
    FLIP = _flip(len(x))

    # From those we define the so-called L₁ and L₂, which are simply the
    # operators that act on the eigenvector itself and the conjugated version,
    # respectively.
    L1 = B + G @ IFFT @ V1 @ FFT
    L2 =     G @ IFFT @ V2 @ FFT @ FLIP  # noqa: E222

    return numpy.block([
        [+ L1.real + L2.real, - L1.imag + L2.imag],
        [+ L1.imag + L2.imag, + L1.real - L2.real]
    ])


def _unroll(x):
    """
    Unroll a complex-valued vector of size n into a real-valued of size 2n.

    Parameters
    ----------
    x : array_like
        a complex-valued vector of length n

    Returns
    -------
    y : array_like
        a real-valued vector of length 2n
    """
    return numpy.append(x.real, x.imag)


def _roll(x):
    """
    Roll the previously unrolled version.

    Parameters
    ----------
    x : array_like
        a real-valued vector of length 2n

    Returns
    -------
    y : array_like
        a complex-valued vector of length n
    """
    n = len(x) // 2
    return x[:n] + 1j * x[n:]


def _order(p):
    """
    A function used for eigenstate ordering.
    """
    v = p[0]
    return abs(v), v


def lopsp(x, beta, gamma, v1, v2):
    """
    Solve a spectral problem for a linearization operator.

        β(k) ψ(k) + γ(k) F{ V₁ ψ(x) + V₂ ψ⁺(x) } = λ ψ(x)

    where ψ(x) is an eigenfunction, ψ⁺(x) is the complex conjugate, β̂(i∂ₓ) is
    the dispersion operator, γ(k) is an arbitrary function of frequency and
    V₁(x) and V₂(x) are two arbitrary potentials.

    Spectral problems like that usually arise when one analyses the spectrum
    of linear modes that exist in a refractive-index potential induced by a
    nonlinear wave, e.g. when analysing soliton stability or dealing with
    generation and scattering of weak dispersive waves by a soliton.

    Parameters
    ----------
    x : array_like
       coordinate grid
    beta : callable with a signature of beta(k)
        the dispersive profile as a function of frequency
    gamma : callable with a signature of gamma(k)
        arbitrary potential scaling as a function of frequency
    v1, v2 : array_like
        potentials as functions of coordinate

    Returns
    -------
    [(eigenvalue, eigenvector)]: array_like
        an list of eigenvalue-eigenvector pair
    """

    L = _lopm(x, beta, gamma, v1, v2)
    vals_, vecs_ = linalg.eig(L)

    vecs_ = _roll(vecs_)
    vals, vecs = [], []
    for n in range(len(vals_)):
        vals.append(vals_[n])
        vecs.append(vecs_[:, n])

    return sorted(zip(vals, vecs), key=_order)


def lopmodes(x, beta, gamma, v1, v2, lambda0, nmodes=10):
    """
    Find linearization operator modes roughly corresponding to the given
    eigenvalue.

    Parameters
    ----------
    x : array_like
       coordinate grid
    beta : callable with a signature of beta(k)
        the dispersive profile as a function of frequency
    gamma : callable with a signature of gamma(k)
        arbitrary potential scaling as a function of frequency
    v1, v2 : array_like
        potentials as functions of coordinate
    lambda0 : real
        the given eigenvalue
    nmodes : int
        number of modes to find
    """

    # Construct the linearization operator
    L = _lopm(x, beta, gamma, v1, v2)

    # We are looking for the eigenvalues using shift-invert method
    I = numpy.eye(2 * len(x))  # noqa: E741
    A = linalg.inv(L - lambda0 * I)

    # Find the eigenstates and convert them into the usual representation
    vals_, vecs_ = sparse.linalg.eigs(A, k=nmodes)
    vecs_ = _roll(vecs_)
    vals_ = 1/vals_ + lambda0

    vals = []
    vecs = []
    for n in range(nmodes):
        vals.append(vals_[n])
        vecs.append(vecs_[:, n])

    return sorted(zip(vals, vecs), key=_order)


def linsfeq(
    t, x,
    u10, u20, k1, k2, pi0,
    beta, beta1, beta2, gamma,
    dt=None, gpc=None, terms=()
):
    """
    Integrate linearized equation for the scattered field.

    This function integrates linearized equation for the scattered field in
    the form of

        ∂ₜ ψ̃ₛ = i β(k) ψ̃ₛ(t, k)
             + i γ(k) F{ 2 |U₁ + U₂|² ψₛ + (U₁ + U₂)² ψₛ⁺ }
       (7.1) + i [β(k) - β₁(k)] Ũ₁(t, k)
       (7.2) + i [β(k) - β₂(k)] Ũ₂(t, k)
       (8.1) + i γ(k) F{ U₁² U₂⁺ }
       (8.2) + i γ(k) F{ U₂² U₁⁺ }
      (12.1) + i γ(k) F{ 2 |U₁|² ψᵢ }
      (12.2) + i γ(k) F{ 2 |U₂|² ψᵢ }
      (13.1) + i γ(k) F{ 2 U₁ U₂⁺ ψᵢ }
      (13.2) + i γ(k) F{ 2 U₂ U₁⁺ ψᵢ }
      (14.1) + i γ(k) F{ U₁² ψᵢ⁺ }
      (14.2) + i γ(k) F{ U₂² ψᵢ⁺ }
        (15) + i γ(k) F{ 2 U₁ U₂ ψᵢ⁺ }

    The notation follows `gnlse`. Numbers next to the terms refer to the
    resonance conditions from the manuscript.

    Parameters
    ----------
    t : array_like
        an array of time points where the integrated function should
        be evaluated
    x : array_like
        the coordinate grid
    u10 : array_like
        the initial condition for U₁ at t[0]
    u20 : array_like
        the initial condition for U₁ at t[0]
    k1 : float
        wavenumber of U₁ component
    k2 : float
        wavenumber of U₂ component
    pi0 : array_like
        the initial condition for ψᵢ at t[0]
    beta : callable with the signature of beta(k)
        the dispersive profile as a function of frequency
    gamma : callable with the signature of gamma(k)
        frequency-dependent gain part of the nonlinear operator
    dt : float, optional
        an optional internal time step used for integration
    gpc : callable with the signature of gcp(t, x, k, u, v)
        an optional callback executed after computing up to the next
        grid point
    terms : tuple of strings
        resonance terms to include in the simulation

    Returns
    -------
    result : an instance of IntegrationResult
    """

    # Pre-allocate the output matrices.
    nt = len(t)
    nx = len(x)

    u = numpy.zeros((nt, nx, 2), dtype=complex)
    v = numpy.zeros((nt, nx, 2), dtype=complex)

    # Put the initial conditions in time and frequency domains into
    # the output matrices
    v0 = numpy.zeros(nx, dtype=complex)
    u[0, :, 0] = v0
    v[0, :, 0] = v0

    # Prepare the frequency scale and evaluate beta on the scale
    k = freqs(x)
    B = beta(k)
    DB1 = beta(k) - beta1(k)
    DB2 = beta(k) - beta2(k)
    G = gamma(k)

    # Prepare the incident radiation spectrum
    qi0 = fft.ifft(pi0, norm="ortho")
    u[0, :, 1] = pi0
    v[0, :, 1] = qi0

    # Prepare the RHS we feed to the solver
    def rhs(t_, v_):
        # Scale the spectrum by the accumulated phase shift due to
        # the dispersion and transform to coordinate space
        exp = numpy.exp(1j * B * (t_ - t[0]))

        u_ = fft.fft(exp * v_, norm="ortho")
        pi = fft.fft(exp * qi0, norm="ortho")

        u1 = u10 * numpy.exp(1j * k1 * (t_ - t[0]))
        u2 = u20 * numpy.exp(1j * k2 * (t_ - t[0]))

        du_ = numpy.zeros_like(u_)

        if "pot.1" in terms or "pot" in terms:
            du_ += 2 * abs(u1 + u2)**2 * u_

        if "pot.2" in terms or "pot" in terms:
            du_ += (u1 + u2)**2 * u_.conjugate()

        if "8.1" in terms or "8" in terms:
            du_ += u1**2 * u2.conjugate()

        if "8.2" in terms or "8" in terms:
            du_ += u2**2 * u1.conjugate()

        if "12.1" in terms or "12" in terms:
            du_ += 2 * abs(u1)**2 * pi

        if "12.2" in terms or "12" in terms:
            du_ += 2 * abs(u2)**2 * pi

        if "13.1" in terms or "13" in terms:
            du_ += 2 * u1 * u2.conjugate() * pi

        if "13.2" in terms or "13" in terms:
            du_ += 2 * u2 * u1.conjugate() * pi

        if "14.1" in terms or "14" in terms:
            du_ += u1**2 * pi.conjugate()

        if "14.2" in terms or "14" in terms:
            du_ += u2**2 * pi.conjugate()

        if "15" in terms:
            du_ += 2 * u1 * u2 * pi.conjugate()

        dv_ = G * fft.ifft(du_, norm="ortho")

        if "7.1" in terms or "7" in terms:
            v1 = fft.ifft(u1, norm="ortho")
            dv_ += DB1 * v1

        if "7.2" in terms or "7" in terms:
            v2 = fft.ifft(u2, norm="ortho")
            dv_ += DB2 * v2

        dv_ = 1j / exp * dv_
        return dv_

    # Configure the solver. We pick ZVODE from ODEPACK which provides
    # error control and adaptive stepping out of the box. Since we are
    # working on a modified spectrum, we can hope the problem to be
    # non-stiff, but it's a good idea *not* to impose the choice of
    # the integration method and error tolerances -- this will be
    # estimated by the solver itself.
    ode = integrate.ode(rhs)
    ode.set_integrator("zvode", rtol=1E-6)

    # Those are the internal loop variables. `t_` holds the current
    # integration time, `v_` is the current modified spectrum.
    t_ = t[0]
    v_ = v0

    timer_start = time.time()
    logging.info("Integrating:")

    for n in range(1, nt):
        _report_progress(n, nt, timer_start)

        # Pick a step size. THIS IS NOT THE ADAPTIVE STEP-SIZE CONTROL
        # --- it is done by the solver. But sometimes the time step
        # given by the time grid is too coarse for the given precision
        # and the solver exhaust maximum number of integration steps
        # while integrating the equation over the grid step. Providing
        # a finer time grid will result in larger output matrices, and
        # for long simulations those might not fit in the memory
        # anymore. Therefore we provide an additional work-around:
        # user can specify an intermediate time step that is used
        # internally by the loop but not saved anywhere.
        dt_ = dt or t[n] - t_

        while t_ < t[n]:
            # Set the initial value for the current integration step.
            ode.set_initial_value(v_, t_)

            # Pick the next step -- either dt_ or until the next grid
            # point.
            if t_ + dt_ > t[n]:
                t_ = t[n]
            else:
                t_ += dt_

            # Integrate for one step and check the solver status.
            v_ = ode.integrate(t_)
            if not ode.successful():
                return IntegrationResult(
                    t, x, k, u, v, error_code=ode.get_return_code())

        # Calculate the proper spectrum and then save the spectrum and
        # the coordinate representation into the output matrices.
        exp = numpy.exp(1j * B * (t_ - t[0]))

        u[n, :, 0] = fft.fft(exp * v_, norm="ortho")
        v[n, :, 0] = fft.fftshift(exp * v_)

        u[n, :, 1] = fft.fft(exp * qi0, norm="ortho")
        v[n, :, 1] = fft.fftshift(exp * qi0)

        if gpc:
            v10 = fft.fftshift(fft.ifft(u10, norm="ortho"))
            v20 = fft.fftshift(fft.ifft(u20, norm="ortho"))
            gpc(
                t, x, k,
                u10 + u20 + u[:, :, 0] + u[:, :, 1],
                v10 + v20 + v[:, :, 0] + v[:, :, 1])

    logging.info("Done!")
    return IntegrationResult(t, x, k, u, v)
