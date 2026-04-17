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
    "makeFocalPlaneFWHMPlot",
    "getFwhmValues",
]


from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

from lsst.afw.cameraGeom import FIELD_ANGLE, Camera
from lsst.obs.lsst import LsstCam

if TYPE_CHECKING:
    from lsst.afw.table import ExposureCatalog


def getFwhmValues(visitSummary: ExposureCatalog) -> dict[int, float]:
    """Extract per-detector FWHM values from a visit summary table.

    Iterates over all LSSTCam detectors and, for those that appear in
    the visit summary, converts the PSF Gaussian sigma (in pixels) to
    a FWHM in arcseconds using the 0.2 arcsec/pixel plate scale.

    Parameters
    ----------
    visitSummary : `lsst.afw.table.ExposureCatalog`
        The visit summary table containing PSF sigma values.

    Returns
    -------
    fwhmValues : `dict` [`int`, `float`]
        A dictionary mapping detector IDs to FWHM values in
        arcseconds. Detectors missing from the summary are omitted.
    """
    camera = LsstCam().getCamera()
    detectors: list[int] = [det.getId() for det in camera]

    fwhmValues: dict[int, float] = {}
    for detectorId in detectors:
        row = visitSummary[visitSummary["id"] == detectorId]

        if len(row) > 0:
            psfSigma = row["psfSigma"][0]
            fwhm = psfSigma * 2.355 * 0.2  # Convert to microns (0.2"/pixel)
            fwhmValues[detectorId] = float(fwhm)

    return fwhmValues


def makeFocalPlaneFWHMPlot(
    fig: plt.Figure,
    ax: plt.Axes,
    fwhmValues: dict[int, float],
    camera: Camera,
    vmin: float | None = None,
    vmax: float | None = None,
    saveAs: str = "",
    title: str = "",
) -> None:
    """Plot a per-detector FWHM map across the focal plane.

    Each detector is drawn as its projected polygon in field-angle
    coordinates and colored by its FWHM value; the mean, median, and
    standard deviation of the FWHM across detectors are annotated in
    the top-right corner of the axes. If ``vmin``/``vmax`` are not
    supplied, the full range of the supplied FWHM values is used.

    Parameters
    ----------
    fig : `matplotlib.figure.Figure`
        The figure object to plot on.
    ax : `matplotlib.axes.Axes`
        The axes to plot on.
    fwhmValues : `dict` [`int`, `float`]
        Dictionary mapping detector IDs to FWHM values in arcseconds.
    camera : `lsst.afw.cameraGeom.Camera`
        The camera geometry object used to look up detector polygons.
    vmin : `float`, optional
        Minimum value for the color map. Defaults to
        ``nanmin(fwhmValues)``.
    vmax : `float`, optional
        Maximum value for the color map. Defaults to
        ``nanmax(fwhmValues)``.
    saveAs : `str`, optional
        If provided, save the figure to this file path.
    title : `str`, optional
        Suptitle for the plot. If empty, no title is set.
    """
    # If vmin and vmax are None, use the min and max of the FWHM values
    if vmin is None:
        vmin = np.nanmin(list(fwhmValues.values()))
    if vmax is None:
        vmax = np.nanmax(list(fwhmValues.values()))

    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap("viridis")

    for detectorId, fwhm in fwhmValues.items():
        detector = camera.get(detectorId)
        corners = detector.getCorners(FIELD_ANGLE)
        cornersDeg = np.rad2deg(corners)

        x = cornersDeg[:, 0]
        y = cornersDeg[:, 1]

        x = np.append(x, x[0])
        y = np.append(y, y[0])

        color = cmap(norm(fwhm))
        ax.fill(x, y, color=color, edgecolor="gray", linewidth=0.5)

        # Compute center of detector for label
        xCenter = np.mean(cornersDeg[:, 0])
        yCenter = np.mean(cornersDeg[:, 1])

        ax.text(xCenter, yCenter, f"{fwhm:.2f}", color="white", fontsize=10, ha="center", va="center")

    # Calculate statistics
    fwhmValuesList = list(fwhmValues.values())
    meanFwhm = np.nanmean(fwhmValuesList)
    medianFwhm = np.nanmedian(fwhmValuesList)
    stdFwhm = np.nanstd(fwhmValuesList)

    statsText = f"Mean: {meanFwhm:.2f}''\nMedian: {medianFwhm:.2f}''\nStd: {stdFwhm:.2f}''"
    ax.text(
        0.98,
        0.98,
        statsText,
        transform=ax.transAxes,
        fontsize=10,
        va="top",
        ha="right",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
    )

    # Add colorbar
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label("FWHM (arcsec)")

    plt.xlabel("Field Angle Y [deg]")
    plt.ylabel("Field Angle X [deg]")
    plt.axis("equal")
    plt.grid(True, alpha=0.3)
    if title:
        fig.suptitle(title, fontsize=18)

    fig.tight_layout()
    if saveAs:
        fig.savefig(saveAs)
