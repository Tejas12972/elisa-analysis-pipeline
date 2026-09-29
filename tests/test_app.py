"""Smoke-test the Streamlit app headlessly (skipped if streamlit is not installed)."""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest

APP = str(Path(__file__).resolve().parents[1] / "app" / "streamlit_app.py")


def _app() -> AppTest:
    return AppTest.from_file(APP, default_timeout=180).run()


def _select(at: AppTest, prefix: str) -> AppTest:
    radio = at.sidebar.radio[0]
    radio.set_value(next(o for o in radio.options if o.startswith(prefix))).run()
    return at


def test_demo_study_renders_all_tabs():
    at = _app()
    assert not at.exception, at.exception
    assert at.metric[0].value == "4/4"
    assert len(at.tabs) == 7


def test_real_il6_dataset():
    at = _select(_app(), "Real: mouse IL-6")
    assert not at.exception, at.exception
    assert at.metric[0].value == "3/4"


def test_real_softmax_dataset():
    at = _select(_app(), "Real: SoftMax")
    assert not at.exception, at.exception
    assert at.metric[0].value == "5/5"
    assert any("no treatment groups" in i.value for i in at.info)


def test_threshold_change_reruns():
    at = _app()
    at.sidebar.slider[0].set_value(10.0).run()  # stricter CV → more flags
    assert not at.exception, at.exception


def test_upload_mode_offers_sample_files():
    at = _select(_app(), "Upload")
    assert not at.exception, at.exception
    assert any("Upload plate exports" in i.value for i in at.info)
