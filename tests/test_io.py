"""Field loading: npz round-trip, gridded/scattered text, .mat files."""

import numpy as np
import pytest

from lcs2d import VelocityField, double_gyre, load_field


class TestNpz:
    def test_roundtrip(self, tmp_path):
        field = double_gyre(nx=41, ny=21)
        path = str(tmp_path / "gyre.npz")
        field.save(path)
        back = load_field(path)
        np.testing.assert_array_equal(back.u, field.u)
        np.testing.assert_array_equal(back.x, field.x)
        assert back.name == "double gyre"

    def test_missing_key_message(self, tmp_path):
        path = str(tmp_path / "bad.npz")
        np.savez(path, x=np.arange(3.0), y=np.arange(3.0))
        with pytest.raises(KeyError, match="missing array"):
            load_field(path)


class TestText:
    def test_gridded_csv(self, tmp_path):
        field = double_gyre(nx=41, ny=21)
        mx, my = np.meshgrid(field.x, field.y)
        table = np.column_stack([mx.ravel(), my.ravel(),
                                 field.u.ravel(), field.v.ravel()])
        path = str(tmp_path / "field.csv")
        np.savetxt(path, table, delimiter=",", header="x,y,u,v", comments="")
        back = load_field(path)
        # complete grid detected -> no resampling, values exact
        assert back.shape == field.shape
        np.testing.assert_allclose(back.u, field.u, atol=1e-12)

    def test_whitespace_delimited(self, tmp_path):
        field = double_gyre(nx=21, ny=11)
        mx, my = np.meshgrid(field.x, field.y)
        table = np.column_stack([mx.ravel(), my.ravel(),
                                 field.u.ravel(), field.v.ravel()])
        path = str(tmp_path / "field.dat")
        np.savetxt(path, table)
        back = load_field(path)
        np.testing.assert_allclose(back.u, field.u, atol=1e-12)

    def test_scattered_resampled(self, tmp_path):
        # random samples of a linear velocity field: linear resampling
        # reproduces it exactly (inside the convex hull)
        rng = np.random.default_rng(2)
        x = rng.uniform(0.0, 1.0, 4000)
        y = rng.uniform(0.0, 1.0, 4000)
        u = 2.0 * x - y
        v = x + y
        table = np.column_stack([x, y, u, v])
        path = str(tmp_path / "scatter.txt")
        np.savetxt(path, table)
        field = load_field(path, grid_res=0.05)
        mx, my = np.meshgrid(field.x, field.y)
        interior = (mx > 0.1) & (mx < 0.9) & (my > 0.1) & (my < 0.9)
        np.testing.assert_allclose(field.u[interior],
                                   (2.0 * mx - my)[interior], atol=1e-6)


class TestMat:
    def test_column_vectors(self, tmp_path):
        from scipy.io import savemat
        field = double_gyre(nx=41, ny=21)
        mx, my = np.meshgrid(field.x, field.y)
        path = str(tmp_path / "field.mat")
        savemat(path, {"x": mx.ravel()[:, None], "y": my.ravel()[:, None],
                       "u": field.u.ravel()[:, None],
                       "v": field.v.ravel()[:, None]})
        back = load_field(path)
        np.testing.assert_allclose(back.u, field.u, atol=1e-12)

    def test_meshgrid_layout(self, tmp_path):
        from scipy.io import savemat
        field = double_gyre(nx=41, ny=21)
        mx, my = np.meshgrid(field.x, field.y)
        path = str(tmp_path / "grid.mat")
        savemat(path, {"X": mx, "Y": my, "U": field.u, "V": field.v})
        back = load_field(path, mat_vars=("X", "Y", "U", "V"))
        np.testing.assert_allclose(back.u, field.u, atol=1e-12)

    def test_helpful_error_on_wrong_names(self, tmp_path):
        from scipy.io import savemat
        path = str(tmp_path / "odd.mat")
        savemat(path, {"foo": np.arange(3.0)})
        with pytest.raises(KeyError, match="--mat-vars"):
            load_field(path)


class TestValidation:
    def test_nonuniform_grid_rejected(self):
        x = np.array([0.0, 1.0, 3.0])
        with pytest.raises(ValueError, match="uniformly spaced"):
            VelocityField(x, x, np.zeros((3, 3)), np.zeros((3, 3)))

    def test_shape_mismatch_rejected(self):
        x = np.arange(4.0)
        with pytest.raises(ValueError, match="shape"):
            VelocityField(x, x, np.zeros((3, 3)), np.zeros((3, 3)))
