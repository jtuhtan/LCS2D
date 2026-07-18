"""FSLE tests: analytic saddle flow and backend parity."""

import numpy as np
import pytest

from lcs2d import NUMBA_AVAILABLE, VelocityField, fsle_field


def saddle_field(a=0.5, n=101):
    x = np.linspace(-1.0, 1.0, n)
    y = np.linspace(-1.0, 1.0, n)
    mx, my = np.meshgrid(x, y)
    return VelocityField(x, y, a * mx, -a * my, name="saddle")


class TestSaddle:
    def test_matches_analytic_value(self):
        # separation grows as delta * exp(a t), so the threshold
        # alpha*delta is crossed at t* = ln(alpha)/a; the discrete search
        # rounds t* up to a whole step
        a, alpha, dt = 0.5, 2.0, 0.01
        res = fsle_field(saddle_field(a), t_max=20.0, dt=dt, res=0.1,
                         delta=0.01, alpha=alpha, margin=0.55,
                         backend="numpy")
        t_star = np.log(alpha) / a
        expected = np.log(alpha) / (np.ceil(t_star / dt) * dt)
        lam = res.values.ravel()
        np.testing.assert_allclose(lam, expected, rtol=1e-3)

    def test_tau_reported(self):
        res = fsle_field(saddle_field(), t_max=20.0, dt=0.01, res=0.1,
                         delta=0.01, alpha=2.0, margin=0.55,
                         backend="numpy")
        tau = res.extra["tau"]
        assert np.isfinite(tau).all()
        np.testing.assert_allclose(np.log(2.0) / tau, res.values)

    def test_never_separated_marked(self):
        # zero flow: no separation ever
        x = np.linspace(0.0, 1.0, 11)
        field = VelocityField(x, x, np.zeros((11, 11)), np.zeros((11, 11)))
        res = fsle_field(field, t_max=1.0, dt=0.1, res=0.2, delta=0.01,
                         backend="numpy")
        assert (res.values == 0.0).all()
        assert np.isnan(res.extra["tau"]).all()
        assert res.params["n_never_separated"] == res.values.size

    def test_alpha_must_exceed_one(self):
        with pytest.raises(ValueError, match="alpha"):
            fsle_field(saddle_field(), t_max=1.0, alpha=1.0)

    @pytest.mark.skipif(not NUMBA_AVAILABLE, reason="numba not installed")
    def test_backends_agree(self):
        field = saddle_field()
        kw = dict(t_max=20.0, dt=0.01, res=0.1, delta=0.01, alpha=2.0,
                  margin=0.55)
        res_np = fsle_field(field, backend="numpy", **kw)
        res_nb = fsle_field(field, backend="numba", **kw)
        np.testing.assert_array_equal(res_np.values, res_nb.values)
        np.testing.assert_array_equal(np.isnan(res_np.extra["tau"]),
                                      np.isnan(res_nb.extra["tau"]))
