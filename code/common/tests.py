import random

from unittest import TestCase
import numpy

from .solver import gnlse_propagate
from .helpers import sech, to_analytic


class SolitonTestCases(TestCase):
    """
    In this test case we check our solver on various known problems of
    soliton propagation. Those are mostly the input/output equality
    tests.
    """

    def test_fundamental_nlse_soliton_propagation(self):
        """
        The most straightforward test case. Take a classic NLSE and
        launch a fundamental soliton. Verify that the output envelope
        is the same as the input envelope.
        """

        def sod(f):
            return - 1/2 * f**2

        def gamma(f):
            return 1

        def kerr(t, x, u):
            return abs(u)**2 * u

        t = numpy.linspace(0, 10, 1000)
        x = numpy.linspace(-20, +20, 1000)
        u0 = 1 / numpy.cosh(x)

        result = gnlse_propagate(t, x, u0, sod, gamma, kerr)
        self.assertTrue(result.successful)

        u = result.u
        self.assertTrue(all(u[0, :] == u0))
        self.assertLess((abs(u[0, :]) - abs(u[-1, :])).max(), 1E-3)


class HelpersTestCase(TestCase):
    """
    In here we test the miscellaneous helper functions.
    """

    def test_sech_is_inverse_of_cosh(self):
        """
        Compare our sech() implementation to the inverse of hyperbolic
        cosine.
        """
        x = numpy.linspace(-100, +100, int(1E4))

        error = 1 / numpy.cosh(x) - sech(x)
        self.assertAlmostEqual(max(abs(error)), 0, places=10)

    def test_analytic_function_properties(self):
        """
        Analytic signals should satisfy the following properties:
        1. Real part of reverse transform yields the original (real) signal.
        2. Real and imaginary component are orthogonal over interval.
        """

        # Construct a random Gaussian function
        amp = 10 * random.random()
        sigma = 10 * random.random()
        omega = 20 * (random.random() - 0.5)
        x0 = 10 * (random.random() - 0.5)
        # ^^^ generates negative frequencies half the time

        x = numpy.linspace(-50, +50, int(1E4))
        u = (
            amp * numpy.exp(- (x - x0)**2 / sigma**2)
                * numpy.exp(-1j * omega * x))
        u = u.real
        z = to_analytic(u).real

        self.assertAlmostEqual(max(abs(u - z.real)), 0, places=10)
        self.assertAlmostEqual(sum(z.real * z.imag), 0, places=10)
