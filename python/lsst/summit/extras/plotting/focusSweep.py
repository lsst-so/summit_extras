# This file is part of summit_extras.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from __future__ import annotations

__all__ = [
    "collectSweepData",
    "fitSweepParabola",
    "plotSweepParabola",
]

from typing import TYPE_CHECKING, Any

import numpy as np
from matplotlib.figure import Figure

from lsst.summit.utils.dateTime import efdTimestampToAstropy
from lsst.summit.utils.efdUtils import getMostRecentRowWithDataBefore

if TYPE_CHECKING:
    from astropy.table import Table

    from lsst.daf.butler import DimensionRecord

PLATESCALE = 0.2  # arcsec / pixel


def collectSweepData(records: list[DimensionRecord], consDbClient: Any, efdClient: Any) -> Table:
    """Build a focus-sweep table of hexapod positions and PSF metrics.

    Queries the consolidated database for quick-look PSF measurements
    for the given visits, then joins the most recent compensated
    camera and M2 hexapod positions from the EFD for each visit.

    Parameters
    ----------
    records : `list` [`lsst.daf.butler.DimensionRecord`]
        Visit records spanning a focus sweep.
    consDbClient : `ConsDbClient`
        Consolidated database client.
    efdClient : `lsst_efd_client.EfdClient`
        Engineering facilities database client.

    Returns
    -------
    table : `astropy.table.Table`
        Table containing per-visit PSF metrics (``sigma``, ``T``,
        ``e1``, ``e2``, ``fwhm``), sequence number, and hexapod
        positions in columns prefixed ``cam_`` and ``m2_`` (``x``,
        ``y``, ``z``, ``u``, ``v``, and ``age`` in seconds).
    """
    visitString = ",".join(str(r.id) for r in records)
    instrument = records[0].instrument
    data = consDbClient.query(
        "SELECT "
        "visit_id as visit_id, "
        "n_inputs as n_inputs, "
        "psf_sigma_median as sigma, "
        "psf_ixx_median as ixx, "
        "psf_iyy_median as iyy, "
        "psf_ixy_median as ixy "
        f"from cdb_{instrument.lower()}.visit1_quicklook WHERE visit_id in ({visitString}) "
        "ORDER BY visit_id"
    )
    data["T"] = data["ixx"] + data["iyy"]
    data["e1"] = (data["ixx"] - data["iyy"]) / data["T"]
    data["e2"] = 2 * data["ixy"] / data["T"]
    data["fwhm"] = np.sqrt(np.log(256)) * PLATESCALE * data["sigma"]

    # Add placeholder seqNum column
    data["seqNum"] = np.nan

    # Placeholder columns for EFD data
    for prefix in ["cam_", "m2_"]:
        for k in ["x", "y", "z", "u", "v", "age"]:
            data[prefix + k] = np.nan

    for row in data:
        rIdx = [r.id for r in records].index(row["visit_id"])
        record = records[rIdx]
        row["seqNum"] = record.seq_num
        for pid, prefix in zip(["MTHexapod:1", "MTHexapod:2"], ["cam_", "m2_"]):
            try:
                efdData = getMostRecentRowWithDataBefore(
                    efdClient,
                    "lsst.sal.MTHexapod.logevent_compensatedPosition",
                    timeToLookBefore=record.timespan.begin,
                    maxSearchNMinutes=3,
                    where=lambda df: df["private_identity"] == pid,
                )
            except ValueError:
                row[prefix + "age"] = np.nan
                for k in ["x", "y", "z", "u", "v"]:
                    row[prefix + k] = np.nan
            else:
                age = record.timespan.begin - efdTimestampToAstropy(efdData["private_efdStamp"])
                row[prefix + "age"] = age.sec
                for k in ["x", "y", "z", "u", "v"]:
                    row[prefix + k] = efdData[k]
    return data


def inferSweepVariable(data: Table) -> str:
    """Heuristically identify which hexapod axis is being swept.

    Compares each hexapod axis's overall RMS against the RMS residual
    from a linear fit vs. sequence number. The axis whose total RMS
    is largest relative to its linear-fit residual is taken to be the
    variable being actively stepped.

    Parameters
    ----------
    data : `astropy.table.Table`
        Table of sweep hexapod motions and PSF measurements, as
        produced by `collectSweepData`.

    Returns
    -------
    varName : `str` or `None`
        Name of the inferred active hexapod variable (e.g.
        ``"cam_z"``), or `None` if inference failed.
    """
    # Examine the ratio of RMS hexapod values to RMS hexapod residuals from a
    # linear fit against seqNum.  If removing the linear term significantly
    # reduces the RMS, that's a good sign this is the active variable.
    stats = {}
    for prefix in ["cam_", "m2_"]:
        for k in ["x", "y", "z", "u", "v"]:
            hexapodValue = data[prefix + k].value.astype(float)
            seqNum = data["seqNum"]
            coefs = np.polyfit(seqNum, hexapodValue, 1)
            resids = np.polyval(coefs, seqNum) - hexapodValue
            stdResids = np.nanstd(resids)
            if stdResids == 0:
                stats[prefix + k] = np.nan
            else:
                stats[prefix + k] = np.nanstd(hexapodValue) / stdResids
    statMax = -np.inf
    varName = None
    for vName, stat in stats.items():
        if stat > statMax:
            varName = vName
            statMax = stat
    if varName is None:
        raise ValueError("Failed to infer swept variable from hexapod data")
    return varName


def fitSweepParabola(data: Table, varName: str) -> dict[str, Any]:
    """Fit a parabola to FWHM vs. the swept hexapod variable.

    Fits ``fwhm`` as a quadratic function of ``data[varName]`` and
    returns the vertex, extremum, fit residual RMS, their
    uncertainties, plus the RMS of the e1 and e2 ellipticity
    components across the sweep.

    Parameters
    ----------
    data : `astropy.table.Table`
        Table of sweep data from `collectSweepData`.
    varName : `str`
        Column name to use as the independent variable for the
        parabolic fit (e.g. ``"cam_z"``).

    Returns
    -------
    fitDict : `dict` [`str`, `object`]
        Dict with keys ``vertex``, ``extremum``, ``rms``,
        ``vertexUncertainty``, ``extremumUncertainty``, ``e1Rms``,
        ``e2Rms``, and the quadratic ``coefs`` array.
    """
    fwhms = data["fwhm"]
    e1s = data["e1"]
    e2s = data["e2"]
    xs = data[varName]
    coefs, cov = np.polyfit(xs, fwhms, 2, cov=True)
    a, b, c = coefs
    vertex = -b / (2 * a)
    resids = np.polyval(coefs, xs) - fwhms
    rms = np.sqrt(np.mean(np.square(resids)))
    extremum = np.polyval(coefs, vertex)

    da = np.sqrt(cov[0, 0])
    db = np.sqrt(cov[1, 1])
    covAB = cov[0, 1]

    # Uncertainty propagation to the vertex x-coordinate (xv = -b/(2a)):
    # d(xv)/da = b/(2a^2), d(xv)/db = -1/(2a).
    vertexUncertainty = np.sqrt(
        (db / (2 * a)) ** 2 + (b * da / (2 * a**2)) ** 2 - (b / (2 * a**2)) * (covAB / a)
    )
    # Uncertainty propagation to the extremum f(xv) = c - b^2/(4a), using
    # the full 3x3 covariance matrix. Partial derivatives evaluated at the
    # vertex: df/da = vertex^2, df/db = vertex, df/dc = 1.
    extremumVariance = (
        vertex**4 * cov[0, 0]
        + vertex**2 * cov[1, 1]
        + cov[2, 2]
        + 2 * vertex**3 * cov[0, 1]
        + 2 * vertex**2 * cov[0, 2]
        + 2 * vertex * cov[1, 2]
    )
    extremumUncertainty = np.sqrt(extremumVariance)

    e1Rms = np.sqrt(np.mean(np.square(e1s)))
    e2Rms = np.sqrt(np.mean(np.square(e2s)))

    return dict(
        vertex=vertex,
        extremum=extremum,
        rms=rms,
        vertexUncertainty=vertexUncertainty,
        extremumUncertainty=extremumUncertainty,
        e1Rms=e1Rms,
        e2Rms=e2Rms,
        coefs=coefs,
    )


def plotSweepParabola(
    data: Table,
    varName: str,
    fitDict: dict[str, Any],
    saveAs: str | None = None,
    figAxes: tuple[Figure, np.ndarray[Any, np.dtype[np.object_]]] | None = None,
) -> None:
    """Plot a focus sweep: hexapod motions, PSF metrics, and the fit.

    Produces a 3x4 grid with the camera and M2 hexapod axis positions
    vs. sequence number on the left, the FWHM and ellipticity
    components vs. sequence number and vs. the swept variable on the
    right, and the fitted parabola overlaid on the FWHM panel.

    Parameters
    ----------
    data : `astropy.table.Table`
        Table of sweep data from `collectSweepData`.
    varName : `str`
        Hexapod variable that was swept (e.g. ``"cam_z"``).
    fitDict : `dict` [`str`, `object`]
        Fit results as returned by `fitSweepParabola`.
    saveAs : `str`, optional
        If provided, save the figure to this file.
    figAxes : `tuple` [`matplotlib.figure.Figure`, `numpy.ndarray`], optional
        Pre-existing figure and 3x4 axes grid to plot into. If not
        provided, a new figure and grid are created.
    """
    xs = data[varName]

    if figAxes is None:
        fig = Figure(figsize=(12, 9))
        axes = fig.subplots(nrows=3, ncols=4)
    else:
        fig, axes = figAxes

    camZAx, m2ZAx, *_ = axes[0]
    camXyAx, m2XyAx, fwhmSeqAx, fwhmVarAx = axes[1]
    camRAx, m2RAx, ellipSeqAx, ellipVarAx = axes[2]

    fig.delaxes(axes[0, 2])
    fig.delaxes(axes[0, 3])

    seqNum = data["seqNum"]
    camZAx.scatter(seqNum, data["cam_z"], c="r")
    camXyAx.scatter(seqNum, data["cam_x"], c="c", label="x")
    camXyAx.scatter(seqNum, data["cam_y"], c="m", label="y")
    camXyAx.legend()
    camRAx.scatter(seqNum, data["cam_u"], c="c", label="Rx")
    camRAx.scatter(seqNum, data["cam_v"], c="m", label="Ry")
    camRAx.legend()

    m2ZAx.scatter(seqNum, data["m2_z"], c="r")
    m2XyAx.scatter(seqNum, data["m2_x"], c="c", label="x")
    m2XyAx.scatter(seqNum, data["m2_y"], c="m", label="y")
    m2XyAx.legend()
    m2RAx.scatter(seqNum, data["m2_u"], c="c", label="Rx")
    m2RAx.scatter(seqNum, data["m2_v"], c="m", label="Ry")
    m2RAx.legend()

    fwhmSeqAx.scatter(seqNum, data["fwhm"], c="r")
    ellipSeqAx.scatter(seqNum, data["e1"], c="c", label="e1")
    ellipSeqAx.scatter(seqNum, data["e2"], c="m", label="e2")
    ellipSeqAx.legend()

    var = data[varName]
    fwhmVarAx.scatter(var, data["fwhm"], c="r")
    xlim = fwhmVarAx.get_xlim()
    xs = np.linspace(xlim[0], xlim[1], 100)
    ys = np.polyval(fitDict["coefs"], xs)
    fwhmVarAx.plot(xs, ys, c="k")

    ellipVarAx.scatter(var, data["e1"], c="c", label="e1")
    ellipVarAx.scatter(var, data["e2"], c="m", label="e2")
    ellipVarAx.legend()

    label = varName.replace("_", " ")
    label = label.replace("u", "Rx")
    label = label.replace("v", "Ry")
    unit = "deg" if "R" in label else "µm"

    for ax in [fwhmVarAx, fwhmSeqAx]:
        ax.set_ylabel("fwhm [arcsec]")

    for ax in [camZAx, m2ZAx]:
        ax.set_ylabel("z [µm]")
    camZAx.set_title("Camera")
    m2ZAx.set_title("M2")

    for ax in [camXyAx, m2XyAx]:
        ax.set_ylabel("x or y [µm]")

    for ax in [camRAx, m2RAx]:
        ax.set_ylabel("Rx or Ry [deg]")

    for ax in [camZAx, camXyAx, camRAx, m2ZAx, m2XyAx, m2RAx, fwhmSeqAx, ellipSeqAx]:
        ax.set_xlabel("seqnum")
        ax.set_xlim(min(seqNum) - 0.5, max(seqNum) + 0.5)

    for ax in [fwhmVarAx, ellipVarAx]:
        ax.set_xlabel(label + "[" + unit + "]")

    for ax in [ellipSeqAx, ellipVarAx]:
        ax.set_ylim(-0.21, 0.21)
        ax.axhline(0, c="k")
        ax.set_ylabel("e1 or e2")

    # Print useful info in the top right
    kwargs: dict[str, Any] = dict(fontsize=10, ha="left", fontfamily="monospace")
    xtext = 0.6
    fig.text(xtext, 0.94, "FWHM fit", **kwargs)
    fig.text(xtext, 0.92, "--------", **kwargs)
    fig.text(xtext, 0.90, f"vertex:    {fitDict['vertex']:.3f} {unit}", **kwargs)
    fig.text(
        xtext,
        0.88,
        f"extremum:  {fitDict['extremum']:.3f} arcsec",
        **kwargs,
    )
    fig.text(xtext, 0.86, f"RMS resid: {fitDict['rms']:.3f} arcsec", **kwargs)

    fig.text(xtext, 0.80, "Ellipticity spread", **kwargs)
    fig.text(xtext, 0.78, "------------------", **kwargs)
    fig.text(xtext, 0.76, f"e1 RMS: {fitDict['e1Rms']:.3f}", **kwargs)
    fig.text(xtext, 0.74, f"e2 RMS: {fitDict['e2Rms']:.3f}", **kwargs)

    fig.tight_layout()
    if saveAs is not None:
        fig.savefig(saveAs)
