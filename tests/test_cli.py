"""End-to-end CLI tests (run through the real argument parser)."""

import os

import numpy as np
import pytest

from lcs2d.cli import main
from lcs2d.result import LCSResult


@pytest.fixture()
def gyre_file(tmp_path):
    path = str(tmp_path / "gyre.npz")
    main(["gyre", "--nx", "101", "--ny", "51", "-o", path, "--quiet"])
    return path


class TestGyre:
    def test_creates_field(self, gyre_file):
        assert os.path.exists(gyre_file)
        with np.load(gyre_file) as data:
            assert data["u"].shape == (51, 101)


class TestFTLE:
    def test_full_run(self, tmp_path, gyre_file):
        out = str(tmp_path / "out" / "ftle")
        main(["ftle", gyre_file, "--t", "5", "--res", "0.05",
              "-o", out, "--quiet"])
        for ext in (".npz", ".csv", ".png"):
            assert os.path.exists(out + ext), f"missing {ext}"
        res = LCSResult.load_npz(out + ".npz")
        assert res.kind == "ftle"
        assert np.isfinite(res.values).all()
        assert res.values.max() > 0

    def test_backward_flag(self, tmp_path, gyre_file):
        out = str(tmp_path / "ftle_bwd")
        main(["ftle", gyre_file, "--t", "5", "--backward", "--res", "0.1",
              "-o", out, "--no-plot", "--no-csv", "--quiet"])
        res = LCSResult.load_npz(out + ".npz")
        assert res.params["t"] == -5.0

    def test_missing_file_polite_error(self, tmp_path, capsys):
        with pytest.raises(SystemExit):
            main(["ftle", str(tmp_path / "nope.npz"), "--t", "5"])


class TestFSLE:
    def test_full_run(self, tmp_path, gyre_file):
        out = str(tmp_path / "fsle")
        main(["fsle", gyre_file, "--t-max", "20", "--res", "0.1",
              "-o", out, "--no-plot", "--quiet"])
        res = LCSResult.load_npz(out + ".npz")
        assert res.kind == "fsle"
        assert (res.values >= 0).all()
        assert "tau" in res.extra


class TestAdvect:
    def test_full_run(self, tmp_path, gyre_file):
        out = str(tmp_path / "traj")
        main(["advect", gyre_file, "--t", "2", "--res", "0.1",
              "--save-every", "10", "-o", out, "--no-plot", "--quiet"])
        with np.load(out + ".npz") as data:
            assert data["x"].shape[0] == data["t"].shape[0]
            assert data["x"].shape[0] > 2


class TestInfo:
    def test_prints_summary(self, gyre_file, capsys):
        main(["info", gyre_file])
        out = capsys.readouterr().out
        assert "double gyre" in out
        assert "suggested dt" in out
