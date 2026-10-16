"""Protein landmarks used to connect DTA heatmaps to molecular structures."""

from dataclasses import dataclass
from typing import Iterable


__all__ = ["ProteinLandmark", "parse_protein_landmarks"]


@dataclass(frozen=True)
class ProteinLandmark:
    """Location of one named protein feature in polar coordinates.

    ``chain_id`` and ``occupancy`` preserve the identifier written by the Tcl
    analysis, for example ``A/1`` for occupancy 1 in chain A.
    """

    chain_id: str
    occupancy: str
    radius: float
    theta: float

    @property
    def identifier(self) -> str:
        """Return the Tcl-style ``chain/occupancy`` identifier."""
        return f"{self.chain_id}/{self.occupancy}"


def parse_protein_landmarks(lines: Iterable[str]) -> list[ProteinLandmark]:
    """Parse landmarks from the existing Tcl protein-coordinate format.

    Non-comment lines contain repeating groups of ``chain/occupancy:``,
    radius, and angle, such as ``A/1: 20.0 90.0 A/2: 22.0 95.0``.
    """
    landmarks = []
    identifiers = set()

    for line_number, line in enumerate(lines, start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        values = stripped.split()
        if values and "/" not in values[0]:
            raise ValueError(
                "Protein coordinates do not contain chain/occupancy "
                "identifiers. Regenerate them with a current DTA Tcl script "
                "before plotting chain landmarks."
            )
        if len(values) % 3 != 0:
            raise ValueError(
                f"Malformed protein-coordinate line {line_number}: expected "
                "repeating identifier, radius, and theta values."
            )

        groups = zip(values[::3], values[1::3], values[2::3], strict=True)
        for raw_identifier, raw_radius, raw_theta in groups:
            landmark = _parse_landmark(
                raw_identifier,
                raw_radius,
                raw_theta,
                line_number,
            )
            identifier = landmark.identifier
            if identifier in identifiers:
                raise ValueError(
                    f"Duplicate protein landmark {identifier!r} on line "
                    f"{line_number}."
                )

            landmarks.append(landmark)
            identifiers.add(identifier)

    return landmarks


def _parse_landmark(
    raw_identifier: str,
    raw_radius: str,
    raw_theta: str,
    line_number: int,
) -> ProteinLandmark:
    """Parse one identifier-radius-theta group."""
    identifier = raw_identifier.removesuffix(":")
    if identifier.count("/") != 1:
        raise ValueError(
            f"Malformed protein landmark {raw_identifier!r} on line "
            f"{line_number}; expected 'chain/occupancy:'."
        )

    chain_id, occupancy = identifier.split("/")
    if not chain_id or not occupancy:
        raise ValueError(
            f"Malformed protein landmark {raw_identifier!r} on line "
            f"{line_number}; chain and occupancy cannot be empty."
        )

    try:
        radius = float(raw_radius)
        theta = float(raw_theta)
    except ValueError as err:
        raise ValueError(
            f"Invalid coordinates for protein landmark {identifier!r} on "
            f"line {line_number}."
        ) from err

    return ProteinLandmark(chain_id, occupancy, radius, theta)
