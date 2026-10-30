"""FTLE tests: analytic saddle flow, double-gyre ridge, backend parity."""

import numpy as np
import pytest

from lcs2d import NUMBA_AVAILABLE, VelocityField, double_gyre, ftle_field


def saddle_field(a=0.5, n=101):
    """u = a x, v = -a y: uniform strain, analytic FTLE = a everywhere."""
    x = np.linspace(-1.0, 1.0, n)
    y = np.linspace(-1.0, 1.0, n)
    mx, my = np.meshgrid(x, y)
    return VelocityField(x, y, a * mx, -a * my, name="saddle")


class TestSaddle:
    def test_matches_analytic_value(self):
        a = 0.5
        res = ftle_field(saddle_field(a), t=1.0, dt=0.01, res=0.1,
                         margin=0.55, backend="numpy")
        # bilinear sampling is exact on a linear field and RK4 truncation
        # is ~ (a dt)^5, so agreement should be near machine precision
        np.testing.assert_allclose(res.values, a, rtol=1e-6)

    def test_backward_time_swaps_axes_but_keeps_value(self):
        # reversing time turns the stretching axis into the contracting
        # one; the maximal exponent is a either way
        a = 0.5
        res = ftle_field(saddle_field(a), t=-1.0, dt=0.01, res=0.1,
                         margin=0.55, backend="numpy")
        np.testing.assert_allclose(res.values, a, rtol=1e-6)

    def test_euler_close_but_less_accurate(self):
        a = 0.5
        res = ftle_field(saddle_field(a), t=1.0, dt=0.01, res=0.1,
                         margin=0.55, integrator="euler", backend="numpy")
        err = np.abs(res.values - a).max()
        assert 1e-6 < err < 0.01  # first-order error, still in the ballpark


@pytest.fixture(scope="module")
def gyre_result():
    return ftle_field(double_gyre(nx=201, ny=101), t=8.0, res=0.02,
                      backend="numpy")


class TestDoubleGyre:

    def test_finite_everywhere(self, gyre_result):
        assert np.isfinite(gyre_result.values).all()

    def test_ridge_on_separatrix(self, gyre_result):
        inner = (gyre_result.x > 0.3) & (gyre_result.x < 1.7)
        col_mean = gyre_result.values[:, inner].mean(axis=0)
        ridge_x = gyre_result.x[inner][np.argmax(col_mean)]
        assert abs(ridge_x - 1.0) < 0.1

    def test_ridge_dominates_gyre_cores(self, gyre_result):
        res = gyre_result
        ridge = res.values[:, np.argmin(np.abs(res.x - 1.0))].max()
        core = res.values[np.argmin(np.abs(res.y - 0.5)),
                          np.argmin(np.abs(res.x - 0.5))]
        assert ridge > 5 * abs(core)

    @pytest.mark.skipif(not NUMBA_AVAILABLE, reason="numba not installed")
    def test_backends_agree(self, gyre_result):
        res_nb = ftle_field(double_gyre(nx=201, ny=101), t=8.0, res=0.02,
                            backend="numba")
        np.testing.assert_allclose(res_nb.values, gyre_result.values,
                                   atol=1e-9)


class TestResultObject:
    def test_npz_roundtrip(self, tmp_path):
        res = ftle_field(saddle_field(), t=1.0, dt=0.05, res=0.2,
                         margin=0.55, backend="numpy")
        path = str(tmp_path / "ftle.npz")
        res.save_npz(path)
        from lcs2d import LCSResult
        back = LCSResult.load_npz(path)
        assert back.kind == "ftle"
        np.testing.assert_array_equal(back.values, res.values)
        assert back.params["t"] == 1.0

    def test_csv_export(self, tmp_path):
        res = ftle_field(saddle_field(), t=1.0, dt=0.05, res=0.2,
                         margin=0.55, backend="numpy")
        path = str(tmp_path / "ftle.csv")
        res.save_csv(path)
        table = np.loadtxt(path, delimiter=",", skiprows=1)
        assert table.shape == (res.values.size, 3)
        np.testing.assert_allclose(table[:, 2], res.values.ravel())


class TestSeedGrid:
    def test_seeds_stay_inside_domain(self):
        from lcs2d import seed_grid
        field = double_gyre(nx=41, ny=21)
        for res in (0.3, 0.07, 0.05, 0.02):
            sx, sy = seed_grid(field, res=res)
            assert sx[0] == 0.0 and sy[0] == 0.0
            assert sx[-1] <= 2.0 and sy[-1] <= 1.0
            np.testing.assert_allclose(np.diff(sx), res, rtol=1e-9)

    def test_exact_multiple_keeps_last_seed(self):
        from lcs2d import seed_grid
        sx, sy = seed_grid(double_gyre(nx=41, ny=21), res=0.02)
        assert sx.size == 101 and sy.size == 51
        assert sx[-1] == 2.0 and sy[-1] == 1.0

    @pytest.mark.parametrize("kw", [dict(res=0.0), dict(res=-0.1),
                                    dict(delta=0.0), dict(delta=-1e-3),
                                    dict(dt=0.0), dict(margin=-0.1)])
    def test_bad_parameters_rejected(self, kw):
        with pytest.raises(ValueError):
            ftle_field(saddle_field(), t=1.0, **kw)
