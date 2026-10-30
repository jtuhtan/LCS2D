"""Advection engine tests against flows with exact solutions."""

import numpy as np
import pytest

from lcs2d import NUMBA_AVAILABLE, VelocityField, advect
from lcs2d.advect import _sample_uv_numpy


def make_field(fu, fv, xlim=(-1.0, 1.0), ylim=(-1.0, 1.0), n=101,
               name="test"):
    x = np.linspace(*xlim, n)
    y = np.linspace(*ylim, n)
    mx, my = np.meshgrid(x, y)
    return VelocityField(x, y, fu(mx, my), fv(mx, my), name=name)


BACKENDS = ["numpy"] + (["numba"] if NUMBA_AVAILABLE else [])


class TestSampling:
    def test_exact_on_linear_field(self):
        # bilinear interpolation reproduces any linear function exactly
        field = make_field(lambda x, y: 2.0 * x - y + 0.5,
                           lambda x, y: x + 3.0 * y)
        rng = np.random.default_rng(0)
        px = rng.uniform(-0.99, 0.99, 500)
        py = rng.uniform(-0.99, 0.99, 500)
        us, vs = _sample_uv_numpy(field.u, field.v, field.x[0], field.y[0],
                                  1.0 / field.dx, 1.0 / field.dy, px, py)
        np.testing.assert_allclose(us, 2.0 * px - py + 0.5, atol=1e-12)
        np.testing.assert_allclose(vs, px + 3.0 * py, atol=1e-12)

    def test_zero_outside_domain(self):
        field = make_field(lambda x, y: np.ones_like(x),
                           lambda x, y: np.ones_like(x))
        px = np.array([-5.0, 5.0, 0.0])
        py = np.array([0.0, 0.0, 5.0])
        us, vs = _sample_uv_numpy(field.u, field.v, field.x[0], field.y[0],
                                  1.0 / field.dx, 1.0 / field.dy, px, py)
        assert us.tolist() == [0.0, 0.0, 0.0]
        assert vs.tolist() == [0.0, 0.0, 0.0]

    def test_nan_treated_as_no_flow(self):
        u = np.ones((5, 5))
        u[2, 2] = np.nan
        field = VelocityField(np.arange(5.0), np.arange(5.0), u, u.copy())
        assert field.n_masked == 1
        assert field.u[2, 2] == 0.0


class TestAdvection:
    @pytest.mark.parametrize("backend", BACKENDS)
    def test_uniform_flow_exact(self, backend):
        field = make_field(lambda x, y: np.full_like(x, 0.3),
                           lambda x, y: np.full_like(x, -0.1),
                           xlim=(0, 10), ylim=(0, 10))
        ex, ey = advect(field, np.array([1.0]), np.array([5.0]), t=2.0,
                        dt=0.1, backend=backend)
        np.testing.assert_allclose(ex, [1.6], atol=1e-12)
        np.testing.assert_allclose(ey, [4.8], atol=1e-12)

    @pytest.mark.parametrize("backend", BACKENDS)
    def test_rotation_returns_to_start(self, backend):
        # u = -y, v = x: circular orbits with period 2*pi
        field = make_field(lambda x, y: -y, lambda x, y: x)
        ex, ey = advect(field, np.array([0.5]), np.array([0.0]),
                        t=2.0 * np.pi, dt=0.01, backend=backend)
        np.testing.assert_allclose([ex[0], ey[0]], [0.5, 0.0], atol=1e-6)

    def test_rk4_beats_euler(self):
        field = make_field(lambda x, y: -y, lambda x, y: x)
        p0 = (np.array([0.5]), np.array([0.0]))
        rx, ry = advect(field, *p0, t=2.0 * np.pi, dt=0.01,
                        integrator="rk4", backend="numpy")
        ex, ey = advect(field, *p0, t=2.0 * np.pi, dt=0.01,
                        integrator="euler", backend="numpy")
        err_rk4 = np.hypot(rx[0] - 0.5, ry[0])
        err_euler = np.hypot(ex[0] - 0.5, ey[0])
        assert err_rk4 < 1e-6
        assert err_euler > 100 * err_rk4

    def test_backward_advection_inverts_forward(self):
        field = make_field(lambda x, y: -y, lambda x, y: x)
        fx, fy = advect(field, np.array([0.4]), np.array([0.1]), t=1.0,
                        dt=0.01, backend="numpy")
        bx, by = advect(field, fx, fy, t=-1.0, dt=0.01, backend="numpy")
        np.testing.assert_allclose([bx[0], by[0]], [0.4, 0.1], atol=1e-9)

    @pytest.mark.skipif(not NUMBA_AVAILABLE, reason="numba not installed")
    def test_backends_agree(self):
        field = make_field(lambda x, y: -y + 0.1 * x, lambda x, y: x)
        rng = np.random.default_rng(1)
        px = rng.uniform(-0.5, 0.5, 200)
        py = rng.uniform(-0.5, 0.5, 200)
        nx, ny_ = advect(field, px, py, t=3.0, dt=0.02, backend="numpy")
        bx, by = advect(field, px, py, t=3.0, dt=0.02, backend="numba")
        np.testing.assert_allclose(nx, bx, atol=1e-12)
        np.testing.assert_allclose(ny_, by, atol=1e-12)

    @pytest.mark.parametrize("backend", BACKENDS)
    def test_trajectory_recording(self, backend):
        field = make_field(lambda x, y: np.ones_like(x),
                           lambda x, y: np.zeros_like(x),
                           xlim=(0, 100), ylim=(0, 100))
        traj = advect(field, np.array([0.0]), np.array([50.0]), t=1.0,
                      dt=0.1, record=True, save_every=3, backend=backend)
        # 10 steps: snapshots at steps 0, 3, 6, 9 plus the final step 10
        assert traj.t.shape == (5,)
        assert traj.x.shape == (5, 1)
        np.testing.assert_allclose(traj.t, [0.0, 0.3, 0.6, 0.9, 1.0],
                                   atol=1e-12)
        np.testing.assert_allclose(traj.x[:, 0], traj.t, atol=1e-12)

    def test_particles_freeze_outside(self):
        field = make_field(lambda x, y: np.ones_like(x),
                           lambda x, y: np.zeros_like(x),
                           xlim=(0, 1), ylim=(0, 1))
        ex, _ = advect(field, np.array([0.9]), np.array([0.5]), t=5.0,
                       dt=0.05, integrator="euler", backend="numpy")
        # leaves through x = 1 and freezes just outside
        assert 1.0 <= ex[0] < 1.1

    def test_zero_time_rejected(self):
        field = make_field(lambda x, y: x, lambda x, y: y)
        with pytest.raises(ValueError, match="nonzero"):
            advect(field, np.array([0.0]), np.array([0.0]), t=0.0)


class TestValidation:
    @pytest.mark.parametrize("dt", [0.0, np.nan, np.inf])
    def test_bad_dt_rejected(self, dt):
        field = make_field(lambda x, y: -y, lambda x, y: x)
        with pytest.raises(ValueError, match="dt"):
            advect(field, np.array([0.1]), np.array([0.1]), t=1.0, dt=dt)

    @pytest.mark.parametrize("backend", BACKENDS)
    def test_nonfinite_positions_rejected(self, backend):
        # a NaN position used to index the grid with garbage
        field = make_field(lambda x, y: -y, lambda x, y: x)
        with pytest.raises(ValueError, match="finite"):
            advect(field, np.array([np.nan]), np.array([0.1]), t=1.0,
                   backend=backend)
