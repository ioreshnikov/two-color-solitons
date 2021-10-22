import random
import time
import logging

from unittest import TestCase, skip
import numpy
from scipy import fft

from .solver import gnlse, _fftm, _ifftm, _flip, lopsp, linsfeq
from .helpers import align_max, freqs, sech, to_analytic


seed = int(time.time())
random.seed(seed)

print("seed is {}".format(seed))


logging.basicConfig(level=logging.DEBUG)


class SolitonTestCase(TestCase):
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

        result = gnlse(t, x, u0, sod, gamma, kerr)
        self.assertTrue(result.successful)

        u = result.u
        self.assertTrue(all(u[0, :] == u0))
        self.assertLess((abs(u[0, :]) - abs(u[-1, :])).max(), 1E-3)

    def test_second_order_soliton_propagation(self):
        """
        This is a slightly more complicated test. We take a
        second-order soliton and let it propagate for 10 periods. We
        check that the soliton shape is restored after each period.
        """

        def sod(f):
            return - 1/2 * f**2

        def gamma(f):
            return 1

        def kerr(t, x, u):
            return abs(u)**2 * u

        period = numpy.pi / 2

        t = numpy.linspace(0, 10 * period, 1000)  # one period = 100 points
        x = numpy.linspace(-20, +20, 1000)
        u0 = 2 / numpy.cosh(x)

        result = gnlse(t, x, u0, sod, gamma, kerr)
        self.assertTrue(result.successful)

        # Pick the initial profile and check it after each period.
        u = result.u
        for n in range(10):
            un = u[n * 100, :]
            self.assertLess((abs(u0) - abs(un)).max(), 1E-3)

    def test_cherenkov_radiation_in_tod_setting(self):
        """
        In here we test that a simple bright soliton launched into a
        third-order dispersion medium develops a radiation tail. The
        maximum of radiation we expect to see near a certain
        frequency predicted by the resonance condition.
        """

        beta3 = 0.15

        def tod(f):
            return - 1/2 * f**2 + 1/6 * beta3 * f**3

        def gamma(f):
            return 1

        def kerr(t, x, u):
            return abs(u)**2 * u

        t = numpy.linspace(0, 10, 1000)
        x = numpy.linspace(-20, +20, 1000)
        u0 = 3 / numpy.cosh(3*x)

        result = gnlse(t, x, u0, tod, gamma, kerr)
        self.assertTrue(result.successful)

        f = result.k
        v = result.v

        # We take a look at the output spectrum
        vout = v[-1, :]

        # We extract the radiation part of the spectrum (it should be
        # somewhere right to f = 10)
        frad = f[f > 10]
        vrad = vout[f > 10]

        # We look for the maximum radiation frequency.
        fidx = numpy.argmax(abs(vrad))
        fmax = frad[fidx]

        # We expect the maximum spectrum to be within 10% from
        # 3/beta3.
        self.assertLess(abs(fmax - 3/beta3) / fmax, 0.1)

    def test_radiation_stays_suppressed_by_a_linear_absorber(self):
        """
        In here we test that Cherenkov radiation from a bright soliton
        stays completely trapped by an absorbing boundary layer and
        that the presence of the boundary layer does not affect
        soliton behaviour.

        The test itself is simple: we first launch the soliton without
        the absorber and the with the absorber.
        """

        beta3 = 0.15

        def sod(f):
            return - 1/2 * f**2 + 1/6 * beta3 * f**3

        def gamma(f):
            return 1

        def kerr(t, x, u):
            return abs(u)**2 * u

        t = numpy.linspace(0, 10, 1000)
        x = numpy.linspace(-20, +20, 1000)
        u0 = 3 / numpy.cosh(3*x)

        # Additional non-reflective absorbing layer near the
        # boundaries of the computational domain
        absorbing_profile = 50 * sech(x + 18) + 50 * sech(x - 18)

        def absorption(t, x, u):
            return 1j * absorbing_profile * u

        result_without_absorber = gnlse(t, x, u0, sod, gamma, kerr)
        result_with_absorber = gnlse(t, x, u0, sod, gamma, kerr, absorption)
        self.assertTrue(result_without_absorber.successful)
        self.assertTrue(result_with_absorber.successful)

        uwa = result_without_absorber.u[-1, :]
        uw = result_with_absorber.u[-1, :]

        inner = abs(x) < 10  # inside and relatively far away from the boundary layer
        outer = abs(x) > 18  # outside the boundary layer

        # The difference in the inner region should be relatively
        # small: it's mostly the wrap-around Cherenkov radiation.
        # Let's say that that difference should be less than 5% of the
        # soliton peak.
        self.assertLess(
            abs(uw[inner] - uwa[inner]).max(),
            0.05 * abs(uw[inner]).max())

        # The outer region in the case with absorbing layer should
        # contain approximately nothing. We define nothing as less
        # than 0.1% of the original Cherenkov radiation.
        self.assertLess(
            abs(uw[outer]).max(),
            1E-3 * abs(uwa[outer]).max())


def random_gaussian(x):
    """
    Generate a random Gaussian pulse.

    Note
    ----
    Assyming x is symmetric.
    """
    a = random.uniform(1, 10)
    c = random.uniform(x.min() / 4, x.max() / 4)
    d = random.uniform(1, 10)
    f = random.uniform(-1, +1)
    u = a * numpy.exp(-(x - c)**2 / d**2) * numpy.exp(-1j * f * x)
    return u


class HelpersTestCase(TestCase):
    """
    In here we test the miscellaneous helper functions. Unfortunatelly this
    doesn't cover everything from the helpers module.
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

        x = numpy.linspace(-50, +50, 2**10)
        u = random_gaussian(x)
        u = u.real
        z = to_analytic(u).real

        self.assertAlmostEqual(max(abs(u - z.real)), 0, places=10)
        self.assertAlmostEqual(sum(z.real * z.imag), 0, places=10)

    def test_align_max_aligns_max(self):
        """
        Test that align_max actually puts the maximum of the function close to
        the origin.
        """
        x = numpy.linspace(-50, +50, 2**10 + 1)
        # ^^^ odd number of points to include 0
        yrandom = random_gaussian(x)
        ycenter, _ = align_max(x, yrandom)

        idx = abs(ycenter).argmax()
        self.assertEqual(x[idx], 0)


class FFTPrimitivesTestCase(TestCase):
    """
    In here we test whether our understanding of convolutions and Fourier
    transforms actually make sense.
    """

    def test_fft_scaling_with_parceval_theorem(self):
        """
        We test that Parseval theorem holds in case of FFT.

        Note
        ----
        One might wonder -- why would you need to test the `numpy.ifft`? You
        wouldn't. This is mostly us trying to figure out whether _our_
        understanding of FFTs scaling is correct. And while we are at it why
        not do it as a test.
        """

        x = numpy.linspace(-100, +100, 1024)
        y = random_gaussian(x)

        z = fft.ifft(y, norm="ortho")

        sx = sum(abs(y)**2)
        sk = sum(abs(z)**2)

        self.assertAlmostEqual(sx, sk, places=6)

    def test_fftm_is_fft(self):
        """
        We test that FFT matrix as returned by `_fftm` is identical to taking
        the standard FFT transform.
        """
        x = numpy.linspace(-100, +100, 1024)
        y = random_gaussian(x)

        FFT = _fftm(len(x))
        z1 = FFT @ y
        z2 = fft.fft(y, norm="ortho")

        self.assertAlmostEqual(abs(z1 - z2).max(), 0, places=6)

    def test_ifftm_is_ifftm(self):
        """
        We test that inverse fft is produced by IFFT.
        """
        x = numpy.linspace(-100, +100, 1024)
        y = random_gaussian(x)

        IFFT = _ifftm(len(x))
        z1 = IFFT @ y
        z2 = fft.ifft(y, norm="ortho")

        self.assertAlmostEqual(abs(z1 - z2).max(), 0, places=6)

    def test_flipm_flips_the_frequencies(self):
        """
        We test that FLIP matrix indeed flips frequencies.
        """
        n = 2 * random.randint(128, 512) - 1
        x = numpy.linspace(-10, +10, n)
        k1 = freqs(x)

        FLIP = _flip(n)
        k2 = FLIP @ k1

        self.assertTrue((k1 == -k2).all())


class SpectralProblemTestCase(TestCase):
    """
    In here we test the spectral problem solver.
    """

    def test_beta_block_performs_differentiation(self):
        """
        In here we test that the block we are using to represent dispersion
        operator actually performs differentiation.
        """

        def beta(k):
            b = -1/2 * k**2
            return b

        f0 = random.uniform(-2, +2)
        X0 = abs(2 * numpy.pi / f0)

        x = numpy.linspace(-2 * X0, 2 * X0, 1024)
        dx = x[1] - x[0]

        px = numpy.exp(-x**2) * numpy.sin(f0 * x)

        k = freqs(x)
        pk = fft.ifft(px, norm="ortho")

        du_analytic = (
            + (-1 + 2 * x**2) * numpy.exp(-x**2) * numpy.sin(f0 * x)
            - 2 * f0 * x * numpy.exp(-x**2) * numpy.cos(f0 * x)
            - 1/2 * f0**2 * numpy.exp(-x**2) * numpy.sin(f0 * x))

        du_gradient = 1/2 * numpy.gradient(numpy.gradient(px, dx), dx)
        du_spectral1 = fft.fft(beta(k) * pk, norm="ortho")

        B = numpy.diag(beta(k))
        FFT = _fftm(len(x))
        IFFT = _ifftm(len(x))
        du_spectral2 = FFT @ B @ IFFT @ px

        self.assertTrue(
            numpy.allclose(du_analytic, du_gradient, atol=dx))
        # ^^^ In theory the difference is around O(Δx²), but that is not
        #     an upper bound. So we use Δx in the test.
        self.assertTrue(
            numpy.allclose(du_analytic, du_spectral1, atol=dx**2))
        self.assertTrue(
            numpy.allclose(du_analytic, du_spectral2, atol=dx**2))
        # ^^^ Those two are actually exponential

    def test_v1_v2_blocks_perform_multiplication(self):
        """
        In here we test that the blocks correponding to multiplication by a
        potential actually perform multiplcation by potential.
        """

        x = numpy.linspace(-50, +50, 2**10 - 1)
        u = random_gaussian(x) + random_gaussian(x)
        v1 = 2 * abs(u)**2
        v2 = u**2

        px = random_gaussian(x)
        pk = fft.ifft(px, norm="ortho")

        V1 = numpy.diag(v1)
        V2 = numpy.diag(v2)

        FFT = _fftm(len(x))
        IFFT = _ifftm(len(x))
        FLIP = _flip(len(x))

        v1pk1 = IFFT @ V1 @ FFT @ pk
        v2pk1 = IFFT @ V2 @ FFT @ FLIP @ pk.conjugate()

        v1pk2 = fft.ifft(v1 * px, norm="ortho")
        v2pk2 = fft.ifft(v2 * px.conjugate(), norm="ortho")

        self.assertTrue(numpy.allclose(v1pk1, v1pk2, atol=1E-6))
        self.assertTrue(numpy.allclose(v2pk1, v2pk2, atol=1E-6))

    def test_eigenvector_orthogonality(self):
        """
        Find the spectrum of a sech(x) potential and then check that the
        eigenfunctions are orthonormal.
        """

        def beta(f):
            return -1/2 * f**2

        def gamma(f):
            return numpy.ones_like(f)

        x = numpy.linspace(-50, +50, 255)
        u = sech(x)
        v1 = 2 * abs(u)**2
        v2 = u**2

        ls, pks = zip(*lopsp(x, beta, gamma, v1, v2))

        for n in range(len(pks)):
            for m in range(0, n):
                p1 = pks[n]
                p2 = pks[m]

                dp = (p1 * p2.conjugate() + p1.conjugate() * p2).sum()
                self.assertLess(dp, 1E-6)

        for n in range(len(pks)):
            pk = pks[n]
            norm = (pk * pk.conjugate() + pk.conjugate() * pk).sum()
            self.assertLess(1 - norm, 1E-6)

    @skip("Hard to come up with a criteria. Modes are ok, but not everywhere.")
    def test_eigenvectors_by_a_residue(self):
        """
        Find the spectrum of a sech(x) potential in spectral representation,
        substitute it back into the coordinate representation and control the
        residue.
        """

        def beta(f):
            return -1/2 * f**2

        def gamma(f):
            return numpy.ones_like(f)

        x = numpy.linspace(-50, +50, 1023)
        dx = x[1] - x[0]

        u = sech(x)
        v1 = 2 * abs(u)**2
        v2 = u**2

        for l, pk in lopsp(x, beta, gamma, v1, v2):
            # Apply dispersion operator explicitly
            px = fft.fft(pk, norm="ortho")
            bx = 1/2 * numpy.gradient(numpy.gradient(px, dx), dx)

            # Calculate potential term explicitly
            ux = v1 * px + v2 * px.conjugate()

            # Finally, calcualte the residue
            rx = bx + ux - l * px

            # Now we transform to spectral domain
            k  = freqs(x, shift=True)
            pk = fft.fftshift(pk)
            rk = fft.fftshift(fft.ifft(rx, norm="ortho"))

            if not numpy.allclose(abs(rk).max() / abs(pk).max(), 0, atol=5E-2):
                import matplotlib.pyplot as plot

                print(abs(rk).max(), abs(pk).max(), abs(rk).max() / abs(pk).max())

                plot.subplot(2, 1, 1)
                plot.plot(x, u.real, label="u")
                plot.plot(x, l + abs(px), label="px")
                plot.plot(x, l + abs(rx), label="rx")
                plot.legend()

                plot.subplot(2, 1, 2)
                plot.plot(k, abs(pk), label="pk")
                plot.plot(k, abs(rk), label="rk")
                plot.legend()

                plot.show()


class LinearizedEquationTestCase(TestCase):
    """
    In here we test that the linearized equations more or less reproduce
    the radiation in the canonical problems we've checked in the
    SolitonTestCase.
    """

    def test_first_soliton_does_not_radiate(self):
        """
        The soliton of NLSE should not radiate if used as U₁.
        """

        def beta(f):
            return -1/2 * f**2

        def beta0(f):
            return numpy.zeros_like(f)

        def gamma(f):
            return 1

        t = numpy.linspace(0, 10, 1000)
        x = numpy.linspace(-20, +20, 1000)
        soliton = 1 / numpy.cosh(x)
        nothing = numpy.zeros_like(soliton)

        result = linsfeq(
            t, x,
            soliton, nothing, 1/2, 1/2, nothing,
            beta, beta, beta0, gamma,
            terms=("7", "8"))
        self.assertTrue(result.successful)
        self.assertAlmostEqual(abs(result.u).max(), 0.0, places=10)

        result = linsfeq(
            t, x,
            nothing, soliton, 1/2, 1/2, nothing,
            beta, beta0, beta, gamma,
            terms=("7", "8"))
        self.assertTrue(result.successful)
        self.assertAlmostEqual(abs(result.u).max(), 0.0, places=10)

    def test_linearized_equations_reproduce_cherenkov_radiation(self):
        """
        In here we integrate a linearized system corresponding to the soliton
        propagating in setting with a third-order dispersion. We then look at
        the spectrum away from the soliton, find the spectral line
        corresponding to Cherenkov radiation and then we check that the
        radiation freuency corresponds to the theoretical prediction.
        """

        beta3 = 0.15

        def beta(f):
            return - 1/2 * f**2 + 1/6 * beta3 * f**3

        def beta2(f):
            return - 1/2 * f**2

        def gamma(f):
            return 1

        def kerr(t, x, u):
            return abs(u)**2 * u

        t = numpy.linspace(0, 5, 500)
        x = numpy.linspace(-20, +20, 1000)
        u0 = 3 / numpy.cosh(3*x)

        nothing = numpy.zeros_like(x)
        k1 = 3**2 / 2

        canonical_result = gnlse(t, x, u0, beta, gamma, kerr)
        self.assertTrue(canonical_result.successful)

        for u1, u2 in [(u0, nothing), (nothing, u0)]:
            linearized_result = linsfeq(
                t, x,
                u1, u2, k1, k1, nothing,
                beta, beta2, beta2, gamma,
                terms=("7", "8"))
            self.assertTrue(linearized_result.successful)

            cv0 = canonical_result.v[0, :]
            cv1 = canonical_result.v[-1, :]
            lv1 = linearized_result.v[-1, :]

            # The most meaningful test we can do here is take the spectrum
            # right of f=15 and check that there is a maximum corresponding to
            # the Cherenkov radiation.

            f = freqs(x, shift=True)

            fw = f > 15
            v1 = linearized_result.v[-1, :]
            idx = abs(v1[fw]).argmax()
            fmax = f[fw][idx]

            self.assertLess(abs(fmax - 3/beta3) / fmax, 0.10)

    def test_linearized_equations_reproduce_scattering(self):
        """
        Integrate a linearized equation approximating a DW scattering process
        and check whether the scattering frequencies satisfy the resonance
        condition.
        """

        beta3 = 0.15

        def beta(f):
            return - 1/2 * f**2 + 1/6 * beta3 * f**3

        def beta2(f):
            return - 1/2 * f**2

        def gamma(f):
            return 1

        t = numpy.linspace(0, 16, 500)
        x = numpy.linspace(-40, +40, 1000)
        f = freqs(x, shift=True)
        u0 = 3 / numpy.cosh(3*x)

        f0 = 2 / beta3
        fi = f0 + 2
        pi = 0.01 * numpy.exp(- ((x + 15)/5)**2 - 1j * fi * x)

        nothing = numpy.zeros_like(x)
        k1 = 3**2 / 2

        for u1, u2 in [(u0, nothing), (nothing, u0)]:
            linearized_result = linsfeq(
                t, x,
                u1, u2, k1, k1, pi,
                beta, beta2, beta2, gamma,
                terms=("12"), dt=0.01)
            self.assertTrue(linearized_result.successful)

            v1 = linearized_result.v[-1, :]
            idx1 = abs(v1[f < f0]).argmax()
            idx2 = abs(v1[f > f0]).argmax()
            f1 = f[f < f0][idx1]
            f2 = f[f > f0][idx2]
            b1 = beta(f1)
            b2 = beta(f2)
            self.assertLess(abs(b1 - b2) / abs(b1), 1E-2)
