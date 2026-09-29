"""Tests for the synthetic plate generator."""

from __future__ import annotations

import numpy as np
import pandas as pd

from immunoassay.io import read_layout, read_plate
from immunoassay.simulate import AnalyteSpec, NoiseSpec, default_demo_study, simulate_study


def test_study_shape_and_balance():
    s = simulate_study(n_per_group=12, seed=1)
    assert len(s.plates) == 2
    truth = s.truth[s.truth["group"].notna()]
    assert truth["sample_id"].nunique() == 36
    per_plate = [(lay["type"] == "sample").sum() // 2 for _, _, _, lay in s.plates]
    assert max(per_plate) - min(per_plate) <= 1
    for _, _, grid, lay in s.plates:
        assert grid.shape == (8, 12)
        assert set(lay["type"]) == {"standard", "blank", "control", "sample"}


def test_deterministic_given_seed():
    a, b = simulate_study(seed=5), simulate_study(seed=5)
    pd.testing.assert_frame_equal(a.truth, b.truth)
    pd.testing.assert_frame_equal(a.plates[0][2], b.plates[0][2])


def test_noiseless_plate_matches_4pl():
    spec = AnalyteSpec()
    noise = NoiseSpec(cv=0, sd_floor=0, od_scale_cv=0, background=0.0)
    s = simulate_study([spec], noise=noise, seed=0)
    _, _, grid, _ = s.plates[0]
    a1 = float(grid.loc["A", 1])
    expected = spec.d + (spec.a - spec.d) / (1 + (spec.top_standard / spec.c) ** spec.b)
    assert abs(a1 - round(expected, 4)) < 1e-4


def test_edge_effect_raises_outer_wells():
    base = NoiseSpec(cv=0, sd_floor=0, od_scale_cv=0)
    edge = NoiseSpec(cv=0, sd_floor=0, od_scale_cv=0, edge_effect=0.2)
    g0 = simulate_study(noise=base, seed=0).plates[0][2]
    g1 = simulate_study(noise=edge, seed=0).plates[0][2]
    # A3 is an outer-row control well; D5 is interior.
    assert np.isclose(float(g1.loc["A", 3]) / float(g0.loc["A", 3]), 1.2, rtol=1e-3)
    assert np.isclose(float(g1.loc["D", 5]), float(g0.loc["D", 5]))


def test_outlier_injection_recorded():
    s = simulate_study(noise=NoiseSpec(outlier_rate=0.1), seed=3)
    assert len(s.injected_outliers) > 0
    assert set(s.injected_outliers.columns) == {"plate_id", "well"}


def test_write_round_trips_through_parsers(tmp_path):
    s = default_demo_study()
    manifest = pd.read_csv(s.write(tmp_path))
    assert len(manifest) == 4
    assert (tmp_path / "control_ranges.csv").exists()
    row = manifest.iloc[0]
    plate = read_plate(tmp_path / row["plate_file"])
    layout = read_layout(tmp_path / row["layout_file"])
    assert (plate["od_status"] != "missing").sum() >= len(layout)
    assert plate["plate_id"].iloc[0] == row["plate_id"]
    ranges = s.control_ranges()
    assert (ranges["low"] < ranges["nominal"]).all() and (ranges["high"] > ranges["nominal"]).all()
