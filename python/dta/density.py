#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Load Tcl count outputs and calculate bead densities and density enrichment.

File readers return counts and polar grid geometry. Calculation functions
convert counts to density and combine enrichment across systems; plotting is
handled separately.
"""
from typing import Literal
from pathlib import Path
from collections import namedtuple
import numpy as np
from dta.utils import validate_path, confirm_objs_are_equal, evenly_divided_quotient
from dta.bin_logic import PolarBinGrid

SysInfo = namedtuple('SysInfo', ['NL', 'NB', 'BoxArea', 'ExpBeadDensity', 'DrDtheta'])


def parse_tcl_dat_file(filepath: str | Path, bulk: bool) -> tuple[np.ndarray, PolarBinGrid, namedtuple]:
    """
    Extract bin counts and bin dimensions from a TCL output .dat file.

    Parameters
    ----------
    filepath : Path or str
        The path to the TCL .dat or .out output file.
    bulk : boolean
        If True, read headerless bulk counts from ``do_get_counts``. If False,
        read per-frame or averaged counts from ``polarDensityBin``.

    Returns
    -------
    counts : ndarray
        If ``bulk`` is True, a flattened one-dimensional integer array.
        Otherwise, a floating-point array of shape
        ``(n_frames, grid.r.n_bins, grid.theta.n_bins)``, indexed by frame,
        radial bin, and angular bin. Averaged files have one frame and may
        contain fractional counts. These are counts, not densities.
    grid : PolarBinGrid or None
        If ``bulk`` is True, return None. Otherwise, return the PolarBinGrid
        object pertaining to this system.
    system_info : SysInfo or None
        Header values ``NL``, ``NB``, ``BoxArea``, ``ExpBeadDensity``, and
        ``DrDtheta``. Returns None for bulk input.

    Raises
    ------
    FileNotFoundError
        If ``filepath`` does not identify an existing file.
    AssertionError
        If the file extension is neither ``.dat`` nor ``.out``.
    ValueError
        If numeric data cannot be read, the header has an unsupported field count,
        or the count data cannot be arranged on the inferred grid.
    IndexError
        If a header field lacks its colon or value, or non-bulk numeric data
        lacks the two-dimensional metadata layout expected by the parser.

    """
    filepath = validate_path(filepath, file=True)
    assert (filepath.suffixes[-1] == '.dat') or (filepath.suffixes[-1] == '.out'), "You must provide the .dat or .out file output from VMD."
    if bulk:
        return np.loadtxt(filepath).astype(int).flatten(), None, None
    unrolled_data = np.loadtxt(filepath, skiprows=1)
    system_info = _parse_system_info(np.loadtxt(filepath, comments=None, max_rows=1, delimiter=',', dtype=str))
    grid = _create_grid_from_dat_file(unrolled_data)
    counts = _package_counts(unrolled_data, grid)
    return counts, grid, system_info


def load_replica_counts(
    root_path: str | Path,
    replicas_list: list[str],
    system_name: str,
    leaflet_id: Literal[1, 2],
    avg: bool = False
) -> list[np.ndarray]:
    """
    Load the counts from all replicas of a single system. Return them as a list.

    Parameters
    ----------
    root_path : str or Path
        The path to directory containing your replica subdirectories.
    replicas_list : list
        At least two replica subdirectory names, in the desired return order.
    system_name : str
        The file stem for your PolarDensityBin outputs.
    leaflet_id : int
        Which leaflet your Site is in. 1=outer and 2=inner.
    avg : bool, optional
        If True, load average counts from ``{system_name}.{leaflet}.avg.dat``.
        If False (default), load per-frame counts from
        ``{system_name}.{leaflet}.dat``. Files are inside each replica directory;
        ``leaflet`` is ``upp`` for 1 and ``low`` for 2.

    Returns
    -------
    list of numpy.ndarray
        Floating-point count arrays in ``replicas_list`` order, each with shape
        ``(n_frames, n_radial_bins, n_angular_bins)``. Averaged files have one
        frame. Counts are not converted to density, and grids are not compared
        across replicas.

    Raises
    ------
    AssertionError
        If ``leaflet_id`` is not 1 or 2, ``root_path`` is not a string or Path
        or does not exist, ``replicas_list`` is not a list with at least two
        entries, or an expected count file is missing.
    ValueError
        If a count file cannot be parsed; see ``parse_tcl_dat_file``.

    """
    assert leaflet_id in [1, 2], "leaflet_id must be 1 (outer) or 2 (inner)."
    assert isinstance(root_path, (str, Path)), "root_path must be a str or Path."
    if isinstance(root_path, str):
        root_path = Path(root_path)
    assert root_path.exists(), f"could not find root_path {root_path}"
    assert isinstance(replicas_list, list), "replicas_list must be a list."
    assert (len(replicas_list) > 1), "Less than 2 replicas found."
    replica_counts_list = []
    leaflet = {1: "upp", 2: "low"}
    for rep in replicas_list:
        if avg:
            fname = root_path.joinpath(rep, f"{system_name}.{leaflet[leaflet_id]}.avg.dat")
        else:
            fname = root_path.joinpath(rep, f"{system_name}.{leaflet[leaflet_id]}.dat")
        assert fname.is_file(), f"could not find file {fname}"
        counts, _, _ = parse_tcl_dat_file(fname, bulk=False)
        replica_counts_list.append(counts)
    return replica_counts_list


def calculate_density(avg_counts: np.ndarray, grid: PolarBinGrid) -> np.ndarray:
    """
    Calculate the average bead density in each bin.

    Parameters
    ----------
    avg_counts : ndarray
        Two-dimensional average bead counts, indexed by radial and then angular
        bin. Shape must be ``(grid.r.n_bins, grid.theta.n_bins)``.
    grid : PolarBinGrid
        The PolarBinGrid object corresponding to this system.

    Returns
    -------
    density : ndarray
        Average counts divided by bin areas, with the same shape as
        ``avg_counts``. Units are beads per squared grid length unit
        (beads per square angstrom for the Tcl output).

    Raises
    ------
    TypeError
        If ``avg_counts`` is not a NumPy array or ``grid`` is not a ``PolarBinGrid``.
    ValueError
        If ``avg_counts`` is not two-dimensional or its shape does not match
        the grid's bin counts.

    """
    if not isinstance(avg_counts, np.ndarray):
        raise TypeError("avg_counts must be an ndarray")
    if len(avg_counts.shape) != 2:
        raise ValueError("avg_counts must be a 2D array.")
    if not isinstance(grid, PolarBinGrid):
        raise TypeError("grid must be a PolarBinGrid.")
    if avg_counts.shape != (grid.r.n_bins, grid.theta.n_bins):
        raise ValueError(f"Shape mismatch. {avg_counts.shape} != {(grid.r.n_bins, grid.theta.n_bins)}.")
    density = avg_counts / grid.bin_areas
    return density


def calculate_density_enrichment(density: np.ndarray, expected_density: float) -> np.ndarray:
    """
    Calculate dimensionless enrichment relative to the expected bead density.

    Parameters
    ----------
    density : ndarray
        Two-dimensional average bead densities, in the same units as
        ``expected_density``.
    expected_density : float or int
        Reference bead density, such as ``SysInfo.ExpBeadDensity`` from
        ``polarDensityBin`` (beads per square angstrom).

    Returns
    -------
    ndarray
        Dimensionless ratios with the same shape as ``density``. A value of 1
        means the bin density equals the expected density.

    Raises
    ------
    AssertionError
        If ``density`` is not a two-dimensional NumPy array or
        ``expected_density`` is not a float or int.

    Notes
    -----
    No positivity or finiteness check is performed on ``expected_density``.
    Zero denominators follow NumPy division behavior, which can produce
    infinities or NaNs and warnings under default NumPy error settings.

    """
    assert isinstance(density, np.ndarray), "density must be a numpy array."
    assert len(density.shape) == 2, "density must be a 2d array."
    assert isinstance(expected_density, (float, int)), "expected_density must be a scalar."
    return density / expected_density


def aggregate_density_enrichment_scores(file_paths: list[str | Path]) -> tuple[np.ndarray, PolarBinGrid]:
    """
    Average density enrichment across systems with matching polar grids.

    Parameters
    ----------
    file_paths : list of str or Path
        Nonempty list of average-count files from ``polarDensityBin``. Each
        file must contain exactly one frame of counts, and all inferred grids
        must be equal.

    Returns
    -------
    all_reps_avg : numpy ndarray
        Dimensionless enrichment of shape
        ``(grid_final.r.n_bins, grid_final.theta.n_bins)``.
        Each system's counts are divided by its bin areas and its own expected
        bead density before taking an equally weighted mean across systems.
        NaNs are ignored per bin; a bin with no non-NaN values remains NaN.
    grid_final : PolarBinGrid
        Shared grid geometry, returned from the first input file for use in
        downstream analysis or plotting.

    Raises
    ------
    TypeError
        If ``file_paths`` is not a list.
    FileNotFoundError
        If an input path does not identify an existing file.
    ValueError
        If ``file_paths`` is empty, an input file cannot be parsed, an input
        contains more than one frame, or the inferred grids differ.

    Notes
    -----
    Systems are weighted equally regardless of their trajectory lengths.
    NumPy emits a warning when a bin contains only NaNs.

    """
    replica_enrichments_list = []
    replica_dims_list = []
    if not isinstance(file_paths, list):
        raise TypeError("file_paths must be a list.")
    if len(file_paths) == 0:
        raise ValueError("file_paths cannot be empty.")
    for rep_path in file_paths:
        rep_path = validate_path(rep_path, file=True)
        counts, grid, system_info = parse_tcl_dat_file(rep_path, bulk=False)
        counts = counts.squeeze(axis=0)
        density_enrichment = calculate_density_enrichment(calculate_density(counts, grid), system_info.ExpBeadDensity)
        replica_enrichments_list.append(density_enrichment)
        replica_dims_list.append(grid)

    confirm_objs_are_equal(replica_dims_list)
    all_reps_avg = np.nanmean(np.stack(tuple(replica_enrichments_list), axis=0), axis=0)
    grid_final = replica_dims_list[0]

    return all_reps_avg, grid_final


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  UTIL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%% #


def _parse_system_info(dat_file_header: list[str]) -> SysInfo:
    """
    Parse system metadata from a current or legacy Tcl output header.

    Parameters
    ----------
    dat_file_header : sequence of str
        Five comma-separated header fields from ``polarDensityBin``, or six
        fields from the legacy format. In the legacy format the fifth field
        is ignored. Each retained field must contain a colon followed by a
        numeric value.

    Returns
    -------
    sysInfo : SysInfo
        Float-valued ``NL`` (lipid count), ``NB`` (bead count), ``BoxArea``
        (square angstroms), ``ExpBeadDensity`` (beads per square angstrom), and
        ``DrDtheta`` (radial width times angular width in radians).

    Raises
    ------
    ValueError
        If the header does not have five or six fields or a value is not numeric.
    IndexError
        If a field lacks the required colon or value.

    """
    # this is done for backwards-compatibility with a discontinued feature
    if len(dat_file_header) == 6:
        NL, NB, BoxArea, ExpBeadDensity, _, DrDtheta = dat_file_header
    elif len(dat_file_header) == 5:
        NL, NB, BoxArea, ExpBeadDensity, DrDtheta = dat_file_header
    else:
        raise ValueError("The header in your .dat file has an unexpected number of entries.")
    NL = _isolate_number_from_header_string(NL)
    NB = _isolate_number_from_header_string(NB)
    BoxArea = _isolate_number_from_header_string(BoxArea)
    ExpBeadDensity = _isolate_number_from_header_string(ExpBeadDensity)
    DrDtheta = _isolate_number_from_header_string(DrDtheta)
    sysInfo = SysInfo(NL, NB, BoxArea, ExpBeadDensity, DrDtheta)
    return sysInfo


def _isolate_number_from_header_string(string: str) -> float:
    """
    Extract the numeric value after the first colon in a header field.

    Parameters
    ----------
    string : str
        A string from the header of the .dat file output by polarDensityBin.

    Returns
    -------
    float
        The number following the colon in the input string.

    Raises
    ------
    IndexError
        If the string has no colon or no value after it.
    ValueError
        If the value cannot be converted to a float.

    Notes
    -----
    Trailing whitespace-separated units are ignored. For density fields such
    as ``Expected Bead Density : 0.02/A^2``, the text after the slash is ignored.

    """
    right_side = string.split(':')[1].strip()
    if "/" in right_side:
        # Expected bead density has a division symbol in it
        return float(right_side.split('/')[0])
    return float(right_side.split()[0])


def _create_grid_from_dat_file(dat_file_contents: np.ndarray) -> PolarBinGrid:
    """
    Infer a regular polar grid from Tcl count-file metadata columns.

    Parameters
    ----------
    dat_file_contents : ndarray
        Two-dimensional numeric data without the header. The first three
        columns give radial lower edge, radial upper edge, and angular bin
        width in degrees. Radial widths are assumed uniform, with angular bins
        spanning a full circle.

    Returns
    -------
    grid : PolarBinGrid
        Grid with radial bounds ``grid.r.lower_bound`` and ``grid.r.upper_bound``
        and bin counts ``grid.r.n_bins`` and ``grid.theta.n_bins``.

    Raises
    ------
    ValueError
        If the radial span or full circle cannot be evenly divided by the
        inferred bin widths, or the inferred geometry fails grid validation.

    """
    r_min = dat_file_contents[0, 0]
    d_r = dat_file_contents[0, 1] - r_min
    r_max = np.max(dat_file_contents[:, 1])
    n_r = evenly_divided_quotient((r_max - r_min), d_r)
    d_theta_degrees = dat_file_contents[0, 2]
    n_theta = evenly_divided_quotient(360, d_theta_degrees)
    grid = PolarBinGrid(r_min, r_max, int(n_r), int(n_theta))
    return grid


def _package_counts(unrolled_data: np.ndarray, grid: PolarBinGrid) -> np.ndarray:
    """
    Package the unrolled data into a 3d array in [time, r, theta] format.

    Parameters
    ----------
    unrolled_data : ndarray
        Two-dimensional numeric data without the header. The first three
        columns contain metadata, followed by one count column per angular bin.
        Rows are grouped by radial bin, with frames in order within each group.
    grid : PolarBinGrid
        Grid supplying radial and angular bin counts. The frame count is
        inferred from the number of rows divided by ``grid.r.n_bins``.

    Returns
    -------
    counts : ndarray
        Floating-point counts of shape
        ``(n_frames, grid.r.n_bins, grid.theta.n_bins)``.

    Raises
    ------
    ValueError
        If the row count is not divisible by the radial bin count, or the
        angular column count does not match the grid.

    Notes
    -----
    A single extra trailing angular column is accepted and removed only when
    every value in it is zero, supporting files from older grid generation.

    """
    n_rows = unrolled_data.shape[0]
    if n_rows % grid.r.n_bins != 0:
        raise ValueError(f"Expected rows to be divisible by {grid.r.n_bins} radial bins, got {n_rows}.")
    nframes = n_rows // grid.r.n_bins

    # chop off columns containing metadata
    unrolled_counts = unrolled_data[:, 3:]
    if unrolled_counts.shape[1] != grid.theta.n_bins:
        if (unrolled_counts.shape[1] == grid.theta.n_bins + 1) and ((unrolled_counts[:, -1] == 0).all()):
            # backwards compatibility for off-by-one grid generation in previous dta version
            unrolled_counts = unrolled_counts[:, :-1]
        else:
            raise ValueError(f"Expected {grid.theta.n_bins} theta columns, got {unrolled_counts.shape[1]}.")

    # 'sideways' because it is in [r, time, theta] format at first
    sideways_counts = np.zeros((grid.r.n_bins, nframes, grid.theta.n_bins))
    for i in range(grid.r.n_bins):
        sideways_counts[i, :, :] = unrolled_counts[(nframes * i):(nframes * (i + 1)), :]

    # swap axes to put it in [time, r, theta] format
    counts = np.swapaxes(sideways_counts, 0, 1)

    return counts
