from unittest import TestCase
import numpy
from .solver import gnlse


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
            return -1/2 * f**2

        def kerr(t, x, u):
            return abs(u)**2 * u

        t = numpy.linspace(0, 10, 1000)
        x = numpy.linspace(-20, +20, 1000)
        u0 = 1 / numpy.cosh(x)

        u, v = gnlse(t, x, u0, sod, kerr)

        self.assertTrue(all(u[0, :] == u0))
        self.assertLess((abs(u[0, :]) - abs(u[-1, :])).max(), 1E-3)
