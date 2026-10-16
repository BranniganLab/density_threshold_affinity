"""Tests for protein-landmark loading and plotting."""

import matplotlib.pyplot as plt
import pytest

from dta.plotting import plot_landmarks_on_panels, plot_protein_landmarks
from dta.protein_landmarks import ProteinLandmark, parse_protein_landmarks
from dta.utils import load_inclusion_coordinates


def test_parse_protein_landmarks_preserves_chain_and_occupancy():
    lines = [
        "# chain/occupancy: r_val theta_val ...\n",
        "A/1: 10.0 20.0 A/2: 11.0 25.0\n",
        "B/1: 12.0 120.0 B/2: 13.0 125.0\n",
    ]

    landmarks = parse_protein_landmarks(lines)

    assert [item.identifier for item in landmarks] == [
        "A/1", "A/2", "B/1", "B/2"
    ]
    assert landmarks[0] == ProteinLandmark("A", "1", 10.0, 20.0)


def test_parse_protein_landmarks_rejects_malformed_rows():
    with pytest.raises(ValueError, match="repeating identifier"):
        parse_protein_landmarks(["A/1: 10.0\n"])


def test_parse_protein_landmarks_rejects_anonymous_legacy_coordinates():
    with pytest.raises(ValueError, match="do not contain chain/occupancy"):
        parse_protein_landmarks(["10.0 20.0 11.0 25.0\n"])


def test_load_inclusion_coordinates_retains_landmark_identity(tmp_path):
    (tmp_path / "Protein_coords_upr.dat").write_text(
        "A/1: 10.0 20.0 A/2: 11.0 25.0\n",
        encoding="utf-8",
    )
    (tmp_path / "Protein_coords_lwr.dat").write_text(
        "B/1: 12.0 120.0\n",
        encoding="utf-8",
    )

    upper, lower = load_inclusion_coordinates(tmp_path)

    assert [item.identifier for item in upper] == ["A/1", "A/2"]
    assert [item.identifier for item in lower] == ["B/1"]


def test_plot_protein_landmarks_colors_and_labels_chains():
    _, ax = plt.subplots(subplot_kw={"projection": "polar"})
    landmarks = [
        ProteinLandmark("A", "1", 10.0, 10.0),
        ProteinLandmark("A", "2", 11.0, 20.0),
        ProteinLandmark("B", "1", 10.0, 100.0),
    ]

    plot_protein_landmarks(
        ax,
        landmarks,
        {"A": "tab:blue", "B": "tab:orange"},
    )

    assert [text.get_text() for text in ax.texts] == ["A", "B"]
    assert len(ax.collections) == 2


def test_plot_landmarks_on_panels_uses_consistent_chain_colors():
    fig, axes = plt.subplots(
        1,
        2,
        subplot_kw={"projection": "polar"},
    )
    panels = [
        [ProteinLandmark("A", "1", 10.0, 10.0)],
        [ProteinLandmark("A", "1", 12.0, 20.0)],
    ]

    plot_landmarks_on_panels(fig, panels)

    assert axes[0].collections[0].get_facecolor().tolist() == \
        axes[1].collections[0].get_facecolor().tolist()
