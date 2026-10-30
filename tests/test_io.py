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


def _grid_table(field):
    mx, my = np.meshgrid(field.x, field.y)
    return mx.ravel(), my.ravel(), field.u.ravel(), field.v.ravel()


class TestTextEdgeCases:
    def test_trailing_comma_and_extra_columns(self, tmp_path):
        field = double_gyre(nx=21, ny=11)
        path = tmp_path / "trailing.csv"
        with open(path, "w") as fh:
            fh.write("x,y,u,v,\n")
            for row in zip(*_grid_table(field)):
                fh.write(",".join(repr(float(val)) for val in row) + ",\n")
        back = load_field(str(path))
        assert back.shape == field.shape
        np.testing.assert_allclose(back.u, field.u, atol=1e-12)

    def test_nan_velocity_stays_masked_on_grid(self, tmp_path):
        # a masked cell must stay a no-flow node of the original grid,
        # not be dropped and then filled in by scattered resampling
        field = double_gyre(nx=21, ny=11)
        x, y, u, v = _grid_table(field)
        u = u.copy()
        u[5 * 21 + 7] = np.nan
        path = str(tmp_path / "masked.csv")
        np.savetxt(path, np.column_stack([x, y, u, v]), delimiter=",",
                   header="x,y,u,v", comments="")
        back = load_field(path)
        assert back.shape == field.shape
        assert back.n_masked == 1
        assert back.u[5, 7] == 0.0 and back.v[5, 7] == 0.0

    def test_duplicate_point_is_not_a_grid(self, tmp_path):
        # 3 unique x * 2 unique y = 6 rows, but (2, 1) is missing and
        # (1, 1) appears twice: must not be reshaped as a complete grid
        x = np.array([0.0, 1.0, 2.0, 0.0, 1.0, 1.0])
        y = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
        path = str(tmp_path / "dup.txt")
        np.savetxt(path, np.column_stack([x, y, x, y]))
        back = load_field(path, grid_res=0.5)
        # linear data resampled: u == x wherever the hull covers it
        mx, _ = np.meshgrid(back.x, back.y)
        covered = back.u != 0.0
        np.testing.assert_allclose(back.u[covered], mx[covered], atol=1e-12)


class TestMatLayouts:
    def test_vector_axes(self, tmp_path):
        from scipy.io import savemat
        field = double_gyre(nx=41, ny=21)
        path = str(tmp_path / "vec.mat")
        savemat(path, {"x": field.x, "y": field.y, "u": field.u,
                       "v": field.v})
        back = load_field(path)
        np.testing.assert_allclose(back.u, field.u, atol=1e-12)

    def test_ndgrid_layout(self, tmp_path):
        from scipy.io import savemat
        field = double_gyre(nx=41, ny=21)
        mx, my = np.meshgrid(field.x, field.y)
        path = str(tmp_path / "nd.mat")
        savemat(path, {"x": mx.T, "y": my.T, "u": field.u.T,
                       "v": field.v.T})
        back = load_field(path)
        np.testing.assert_allclose(back.u, field.u, atol=1e-12)

    def test_matrix_columns(self, tmp_path):
        # the legacy seedXY / seedUV layout: N x 2 matrices
        from scipy.io import savemat
        field = double_gyre(nx=41, ny=21)
        x, y, u, v = _grid_table(field)
        path = str(tmp_path / "cols.mat")
        savemat(path, {"XY": np.column_stack([x, y]),
                       "UV": np.column_stack([u, v])})
        back = load_field(path, mat_vars=("XY:1", "XY:2", "UV:1", "UV:2"))
        np.testing.assert_allclose(back.u, field.u, atol=1e-12)

    def test_bad_column_rejected(self, tmp_path):
        from scipy.io import savemat
        path = str(tmp_path / "cols.mat")
        savemat(path, {"XY": np.ones((5, 2))})
        with pytest.raises(ValueError, match="column 3"):
            load_field(path, mat_vars=("XY:1", "XY:2", "XY:3", "XY:1"))


class TestScatteredValidation:
    def test_collinear_samples_rejected(self):
        with pytest.raises(ValueError, match="span no area"):
            VelocityField.from_scattered(np.arange(5.0), np.zeros(5),
                                         np.ones(5), np.ones(5))
