import logging
import time

from scipy import fft
from scipy import integrate
import numpy


def _report_progress(n, nt, timer_start):
    """
    Report integration progress to the log.
    """

    progress = (n / nt)
    timer_current = time.time()
    timer_elapsed = timer_current - timer_start
    timer_eta = timer_elapsed * (1 - progress) / progress

    logging.debug(
        "  %5.2f%%  %6.1fs elapsed  %6.1fs left",
        100 * progress, timer_elapsed, timer_eta)


def gnlse(t, x, u0, beta, nonlin, dt=None):
    """
    Integrate a GNLSE using the integrating factor method.

    Parameters
    ----------
    t : array_like
        an array of time points where the integrated function should
        be evaluated
    x : array_like
        the coordinate grid
    u0 : array_like
        the initial condition at t[0]
    beta : callable with the signature of beta(f)
        the dispersive profile as a function of frequency
    nonlin : callable with the signature of nonlin(t, x, u)
        the nonlinear operator in the equation

    Returns
    -------
    (u, v) : a tuple of complex-valued matrices of shape len(t)×len(x)
        output matrices of the integrated field in coordinate-domain
        and spectral-domain representations. Row number is position in
        the time array, column number is the position in coordinate/frequency
        array.

    Raises
    ------
    RuntimeError
        when an integration step fails. Please examine the log for
        the error code
    """

    # Pre-allocate the output matrices.
    nt = len(t)
    nx = len(x)

    u = numpy.zeros((nt, nx), dtype=complex)
    v = numpy.zeros((nt, nx), dtype=complex)

    # Put the initial conditions in time and frequency domains into
    # the output matrices
    v0 = fft.ifft(u0)
    u[0, :] = u0
    v[0, :] = fft.fftshift(v0)

    # Prepare the frequency scale and evaluate beta on the scale
    f = 2 * numpy.pi * fft.fftfreq(nx, x[1] - x[0])
    D = beta(f)

    # Prepare the RHS we feed to the solver
    def rhs(t_, v_):
        # Scale the spectrum by the accumulated phase shift due to
        # the dispersion and transform to coordinate space
        exp = numpy.exp(1j * D * t_)
        u_ = fft.fft(exp * v_)

        # Apply the nonlinear operator, transform back to the modified
        # spectrum and return
        return 1j / exp * fft.ifft(nonlin(t, x, u_))

    # Configure the solver. We pick ZVODE from ODEPACK which provides
    # error control and adaptive stepping out of the box. Since we are
    # working on a modified spectrum, we can assume the problem to be
    # non-stiff and we pick 4-th order implicit Adams integrator as the
    # stepper function.
    ode = integrate.ode(rhs, jac=None)  # XXX: Should we pass a Jacobian?
    ode.set_integrator(
        "zvode",
        method="adams",
        atol=1E-10,
        rtol=1E-5,
        order=4)

    # Those are the internal loop variables. `t_` holds the current
    # integration time, `v_` is the current modified spectrum.
    t_ = 0
    v_ = v0

    timer_start = time.time()
    logging.debug("Integrating:")

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
                raise RuntimeError(
                    "integration failed with return code {}"
                    .format(ode.get_return_code()))

        # Calculate the proper spectrum and then save the spectrum and
        # the coordinate representation into the output matrices.
        exp = numpy.exp(1j * D * t_)
        u[n, :] = fft.fft(exp * v_)
        v[n, :] = fft.fftshift(exp * v_)

    logging.debug("%.2fs elapsed", time.time() - timer_start)
    return u, v
