"""Analytic test fields: the coastal cape-and-bay residual circulation."""

import numpy as np
import pytest

import lcs2d
from lcs2d.advect import _field_args, _sample_uv_numpy
from lcs2d.fields import _coastline


@pytest.fixture(scope="module")
def coast():
    return lcs2d.coastal_bay(res=0.1)


def _at(field, x, y):
    return int(round((y - field.y[0]) / field.dy)), \
        int(round((x - field.x[0]) / field.dx))


def test_land_is_masked(coast):
    mx, my = np.meshgrid(coast.x, coast.y)
    land = my < _coastline(mx)
    assert coast.n_masked == land.sum()
    assert np.all(coast.u[land] == 0.0) and np.all(coast.v[land] == 0.0)


def test_far_field_is_the_alongshore_current(coast):
    j, i = _at(coast, 0.0, 7.5)
    assert coast.u[j, i] == pytest.approx(0.36, rel=0.1)
    assert abs(coast.v[j, i]) < 0.05


def test_divergence_free_offshore(coast):
    dudx = np.gradient(coast.u, coast.dx, axis=1)
    dvdy = np.gradient(coast.v, coast.dy, axis=0)
    mx, my = np.meshgrid(coast.x, coast.y)
    sea = my > _coastline(mx) + 0.5
    scale = np.abs(np.gradient(coast.u, coast.dy, axis=0))[sea].max()
    assert np.abs(dudx + dvdy)[sea].max() < 1e-6 * scale


def test_no_flow_through_the_coast(coast):
    # 200 m offshore (along the normal), the flow runs along the shore
    xs = np.linspace(-5.0, 20.0, 400)
    yc = _coastline(xs)
    slope = np.gradient(yc, xs)
    norm = np.hypot(1.0, slope)
    px, py = xs - 0.2 * slope / norm, yc + 0.2 / norm
    u, v = _sample_uv_numpy(*_field_args(coast), px, py)
    normal = np.abs(v - slope * u) / norm
    along = np.abs(u + slope * v) / norm
    assert normal.max() < 0.25 * along.max()
    assert np.median(normal / along) < 0.1


@pytest.mark.parametrize("x0, y0, sign", [(2.0, 2.5, +1),   # west of cape
                                          (5.6, 2.0, -1),   # lee of cape
                                          (9.4, 2.3, +1)])  # bay gyre
def test_eddy_rotation(coast, x0, y0, sign):
    dvdx = np.gradient(coast.v, coast.dx, axis=1)
    dudy = np.gradient(coast.u, coast.dy, axis=0)
    j, i = _at(coast, x0, y0)
    assert np.sign(dvdx[j, i] - dudy[j, i]) == sign


def test_bad_resolution():
    with pytest.raises(ValueError):
        lcs2d.coastal_bay(res=0.0)
