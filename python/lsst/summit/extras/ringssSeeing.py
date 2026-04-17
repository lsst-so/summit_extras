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

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING
from zoneinfo import ZoneInfo

import pandas as pd
from astropy.time import Time, TimeDelta
from matplotlib.dates import DateFormatter, num2date

from lsst.summit.utils.dateTime import getDayObsEndTime, getDayObsStartTime
from lsst.summit.utils.efdUtils import getEfdData, getMostRecentRowWithDataBefore
from lsst.utils.plotting.figures import make_figure

HAS_EFD_CLIENT = True
try:
    from lsst_efd_client import EfdClient
except ImportError:
    HAS_EFD_CLIENT = False

if TYPE_CHECKING:
    from datetime import datetime

    from matplotlib.figure import Figure
    from pandas import DataFrame, Series

    from lsst.daf.butler import Butler, DataCoordinate, DimensionRecord

RINGSS_TOPIC = "lsst.sal.ESS.logevent_ringssMeasurement"


@dataclass
class SeeingConditions:
    """A single snapshot of RINGSS seeing conditions.

    Values are either taken directly from a RINGSS EFD row or linearly
    interpolated between two rows bracketing a requested time. All
    numeric attributes default to NaN when the underlying column was
    missing.

    Attributes
    ----------
    timestamp : `astropy.time.Time`
        The time to which these conditions apply. For an interpolated
        instance this is the midpoint between the two source rows.
    eRMS : `float`
        Residual defocus signature RMS reported by RINGSS.
    flux : `float`
        Measured flux from the RINGSS target star, in instrumental
        units.
    fwhmFree : `float`
        Free-atmosphere seeing FWHM, in arcsec.
    fwhmScintillation : `float`
        Scintillation-derived total seeing FWHM, in arcsec.
    fwhmSector : `float`
        Profile-weighted sector seeing FWHM, in arcsec.
    hrNum : `float`
        RINGSS star HR number.
    tau0 : `float`
        Atmospheric coherence time, in ms.
    theta0 : `float`
        Isoplanatic angle, in arcsec.
    totalVariance : `float`
        Total variance of the wavefront signal.
    profile0m : `float`
        Turbulence integral at the ground layer (0 m).
    profile250m : `float`
        Turbulence integral at 250 m.
    profile500m : `float`
        Turbulence integral at 500 m.
    profile1000m : `float`
        Turbulence integral at 1 km.
    profile2000m : `float`
        Turbulence integral at 2 km.
    profile4000m : `float`
        Turbulence integral at 4 km.
    profile8000m : `float`
        Turbulence integral at 8 km.
    profile16000m : `float`
        Turbulence integral at 16 km.
    windSpeed : `float`
        Wind speed at the turbulence layer, in m/s.
    zenithDistance : `float`
        Zenith distance of the RINGSS target at the measurement time,
        in degrees.
    """

    timestamp: Time

    eRMS: float = float("nan")
    flux: float = float("nan")
    fwhmFree: float = float("nan")
    fwhmScintillation: float = float("nan")
    fwhmSector: float = float("nan")
    hrNum: float = float("nan")
    tau0: float = float("nan")
    theta0: float = float("nan")
    totalVariance: float = float("nan")
    profile0m: float = float("nan")
    profile250m: float = float("nan")
    profile500m: float = float("nan")
    profile1000m: float = float("nan")
    profile2000m: float = float("nan")
    profile4000m: float = float("nan")
    profile8000m: float = float("nan")
    profile16000m: float = float("nan")
    windSpeed: float = float("nan")
    zenithDistance: float = float("nan")

    @property
    def isoparalacticAngle(self) -> float:
        """Isoparalactic angle, in arcsec.

        Returns
        -------
        angle : `float`
            The isoparalactic angle.

        Raises
        ------
        NotImplementedError
            Always raised; there is no corresponding RINGSS field
            available in the EFD topic yet.
        """
        raise NotImplementedError("isoparalacticAngle is not available from the RINGSS EFD topic")

    @property
    def starName(self) -> str:
        """HD-prefixed name of the RINGSS target star.

        Returns
        -------
        name : `str`
            The star name, e.g. ``"HD12345"``.
        """
        return f"HD{self.starName}"

    @property
    def seeing(self) -> float:
        """Scintillation-derived total seeing FWHM, in arcsec.

        Returns
        -------
        seeing : `float`
            Alias for `fwhmScintillation`.
        """
        return self.fwhmScintillation

    @property
    def seeing2(self) -> float:
        """Profile-weighted seeing FWHM, in arcsec.

        Returns
        -------
        seeing : `float`
            Alias for `fwhmSector`.
        """
        return self.fwhmSector

    @property
    def freeAtmosphericSeeing(self) -> float:
        """Free-atmosphere seeing FWHM, in arcsec.

        Returns
        -------
        seeing : `float`
            Alias for `fwhmFree`.
        """
        return self.fwhmFree

    @property
    def groundLayerSeeing(self) -> float:
        """Ground-layer seeing FWHM in arcsec.

        Returns
        -------
        seeing : `float`
            The seeing contribution from the ground layer.

        Raises
        ------
        NotImplementedError
            Always raised; the conversion from the turbulence
            integral at 0 m to arcsec is not yet implemented.
        """
        raise NotImplementedError("I need the transormation to get this value in arcsec")

    def __init__(self, rows: list[Series]) -> None:
        # Ensure we have valid data
        assert len(rows) >= 1, "Must provide some data!"
        assert len(rows) <= 2, "Provided data must be either or two rows."

        # Extract timestamp from index
        timestamps = [row.name for row in rows]
        self.timestamp = Time(timestamps[0])  # Use the first timestamp as default

        if len(rows) == 1:
            # If only one row, use its values directly
            row = rows[0]

            self.eRMS = row.get("eRMS", float("nan"))
            self.flux = row.get("flux", float("nan"))
            self.fwhmFree = row.get("fwhmFree", float("nan"))
            self.fwhmScintillation = row.get("fwhmScintillation", float("nan"))
            self.fwhmSector = row.get("fwhmSector", float("nan"))
            self.hrNum = row.get("hrNum", float("nan"))
            self.tau0 = row.get("tau0", float("nan"))
            self.theta0 = row.get("theta0", float("nan"))
            self.totalVariance = row.get("totalVariance", float("nan"))
            self.profile0m = row.get("turbulenceProfiles0", float("nan"))
            self.profile250m = row.get("turbulenceProfiles1", float("nan"))
            self.profile500m = row.get("turbulenceProfiles2", float("nan"))
            self.profile1000m = row.get("turbulenceProfiles3", float("nan"))
            self.profile2000m = row.get("turbulenceProfiles4", float("nan"))
            self.profile4000m = row.get("turbulenceProfiles5", float("nan"))
            self.profile8000m = row.get("turbulenceProfiles6", float("nan"))
            self.profile16000m = row.get("turbulenceProfiles7", float("nan"))
            self.windSpeed = row.get("wind", float("nan"))
            self.zenithDistance = row.get("zenithDistance", float("nan"))
        else:
            # Interpolate between two rows
            t1, t2 = timestamps
            # Calculate the midpoint timestamp for interpolation
            t = t1 + (t2 - t1) / 2
            self.timestamp = Time(t)

            # Weight for interpolation (0 to 1)
            w = (t - t1) / (t2 - t1)

            row1, row2 = rows[0], rows[1]

            # Helper function for interpolation
            def interpolate(col: str) -> float:
                if col in row1 and col in row2:
                    v1, v2 = row1.get(col, float("nan")), row2.get(col, float("nan"))
                    return v1 + (v2 - v1) * w
                return float("nan")

            self.eRMS = interpolate("eRMS")
            self.flux = interpolate("flux")
            self.fwhmFree = interpolate("fwhmFree")
            self.fwhmScintillation = interpolate("fwhmScintillation")
            self.fwhmSector = interpolate("fwhmSector")
            self.hrNum = interpolate("hrNum")
            self.tau0 = interpolate("tau0")
            self.theta0 = interpolate("theta0")
            self.totalVariance = interpolate("totalVariance")
            self.profile0m = interpolate("turbulenceProfiles0")
            self.profile250m = interpolate("turbulenceProfiles1")
            self.profile500m = interpolate("turbulenceProfiles2")
            self.profile1000m = interpolate("turbulenceProfiles3")
            self.profile2000m = interpolate("turbulenceProfiles4")
            self.profile4000m = interpolate("turbulenceProfiles5")
            self.profile8000m = interpolate("turbulenceProfiles6")
            self.profile16000m = interpolate("turbulenceProfiles7")
            self.windSpeed = interpolate("wind")
            self.zenithDistance = interpolate("zenithDistance")

    def __repr__(self) -> str:
        return (
            f"SeeingConditions @ {self.timestamp.isot}\n"
            f"  Seeing          = XXX define me!\n"
            f'  Free Atm Seeing = {self.fwhmFree:.2f}"\n'
            # XXX replace this with self.groundLayerSeeing
            f"  Ground layer    = {self.profile0m}\n"
            f"  Wind speed    = {self.windSpeed:.2f}m/s\n"
        )


class RingssSeeingMonitor:
    """Query and plot RINGSS seeing data from the EFD.

    Parameters
    ----------
    efdClient : `lsst_efd_client.EfdClient`
        Client used to query the EFD for RINGSS data.
    warningThreshold : `float`, optional
        Interpolation interval, in seconds, above which a warning is
        logged when returning interpolated seeing conditions.
    errorThreshold : `float`, optional
        Interpolation interval, in seconds, above which an interpolated
        value is refused outright. Also used as the search half-window
        when looking for the nearest RINGSS measurement.
    """

    def __init__(
        self, efdClient: EfdClient, warningThreshold: float = 300, errorThreshold: float = 600
    ) -> None:
        self.warningThreshold = warningThreshold
        self.errorThreshold = errorThreshold
        self.log = logging.getLogger(__name__)
        self.client = efdClient

    def getSeeingAtTime(self, time: Time) -> SeeingConditions:
        """Get the seeing conditions at a specific time.

        Parameters
        ----------
        time : `Time`
            The time at which to get the seeing conditions.

        Returns
        -------
        seeing : `SeeingConditions`
            The seeing conditions at the specified time.
        """
        begin = time - TimeDelta(self.errorThreshold, format="sec")
        end = time + TimeDelta(self.errorThreshold, format="sec")
        data = getEfdData(self.client, RINGSS_TOPIC, begin=begin, end=end, warn=False)
        if data.empty:
            raise ValueError(
                f"Failed to find a RINGSS seeing measurement within"
                f" {self.errorThreshold / 60:.1f} minutes of {time.isot}"
            )

        # Ensure the timestamp is sorted
        data.sort_index(inplace=True)

        # Check if the *exact* time exists - seems unlikely, but need to check
        if time in data.index:
            row = data.loc[time]
            return SeeingConditions(
                rows=[row],
            )

        # Convert the astropy Time object to a timezone-aware pandas Timestamp
        time_datetime = pd.Timestamp(time.datetime).tz_localize("UTC")

        # Use the timezone-aware timestamp for comparison
        earlier = (
            data[data.index < time_datetime].iloc[-1] if not data[data.index < time_datetime].empty else None
        )
        later = (
            data[data.index > time_datetime].iloc[0] if not data[data.index > time_datetime].empty else None
        )

        if later is None and earlier is not None and (time - earlier.name).sec < self.errorThreshold:
            self.log.info("Returning the last available value.")
            return SeeingConditions(
                rows=[earlier],
            )

        if earlier is None or later is None:
            raise ValueError("Cannot interpolate: insufficient data before or after the requested time.")

        # Check time difference: to log warnings/raise as necessary
        earlierTime = earlier.name
        laterTime = later.name
        interval = (laterTime - earlierTime).seconds

        if interval > self.errorThreshold:
            raise ValueError(
                f"Requested time {time.isot} would require interpolating between values more "
                f"than {self.errorThreshold} apart: {interval:.2f} seconds."
            )
        if interval > self.warningThreshold:
            self.log.warning(
                f"Interpolating between values more than {self.warningThreshold / 60:.1f} mins apart."
            )

        return SeeingConditions(
            rows=[earlier, later],
        )

    def getMostRecentTimestamp(self) -> Time:
        """Get the most recent timestamp for which seeing data is available.

        Returns
        -------
        time : `Time`
            The most recent timestamp with seeing data.
        """
        now = Time.now() + TimeDelta(10, format="sec")
        row = getMostRecentRowWithDataBefore(
            self.client,
            RINGSS_TOPIC,
            now,
            maxSearchNMinutes=self.errorThreshold / 60,
        )
        return Time(row.name)

    def getMostRecentSeeing(self) -> SeeingConditions:
        """Get the most recent seeing conditions.

        Returns
        -------
        seeing : `SeeingConditions`
            The most recent seeing conditions.

        Raises
        ------
        ValueError
            If no data is available in the EFD within the error threshold.
        """
        now = Time.now() + TimeDelta(10, format="sec")
        try:
            row = getMostRecentRowWithDataBefore(
                self.client,
                RINGSS_TOPIC,
                now,
                maxSearchNMinutes=self.errorThreshold / 60,
            )
        except ValueError as e:
            raise ValueError("Could not get SeeingConditions - no data available in the EFD.") from e

        # XXX add validity checks here

        return SeeingConditions([row])

    def getSeeingForDataId(self, butler: Butler, dataId: DataCoordinate) -> SeeingConditions:
        """Get the seeing conditions for a specific data ID.

        Parameters
        ----------
        butler : `Butler`
            The Butler instance to query.
        dataId : `DataCoordinate`
            The data ID for which to get the seeing conditions.

        Returns
        -------
        seeing : `SeeingConditions`
            The seeing conditions for the specified data ID.
        """
        (expRecord,) = butler.registry.queryDimensionRecords("exposure", dataId=dataId)
        return self.getSeeingForExpRecord(expRecord)

    def getSeeingForExpRecord(self, expRecord: DimensionRecord) -> SeeingConditions:
        """Get the seeing conditions for a specific exposure record.

        Parameters
        ----------
        expRecord : `DimensionRecord`
            The exposure record for which to get the seeing conditions.

        Returns
        -------
        seeing : `SeeingConditions`
            The seeing conditions for the specified exposure record.
        """
        midPoint = expRecord.timespan.begin + TimeDelta(expRecord.exposure_time / 2, format="sec")
        return self.getSeeingAtTime(midPoint)

    def plotSeeingForDayObs(
        self, dayObs: int, addMostRecentBox: bool = True, fig: Figure | None = None
    ) -> Figure:
        """Plot the seeing conditions for a specific day observation.

        Parameters
        ----------
        dayObs : `int`
            The dayObs to plot the seeing for, in YYYYMMDD format.
        addMostRecentBox : `bool`, optional
            Whether to add a box with the most recent seeing conditions.
        fig : `Figure`, optional
            The figure to plot on. If None, a new figure will be created.

        Returns
        -------
        fig : `Figure`
            The figure with the plotted seeing conditions.
        """
        startTime = getDayObsStartTime(dayObs)
        endTime = getDayObsEndTime(dayObs)
        data = getEfdData(self.client, RINGSS_TOPIC, begin=startTime, end=endTime, warn=False)
        fig = self.plotSeeing(data, addMostRecentBox=addMostRecentBox, fig=fig)
        return fig

    def plotSeeing(
        self, dataframe: DataFrame, addMostRecentBox: bool = True, fig: Figure | None = None
    ) -> Figure:
        """Plot the seeing conditions from a DataFrame.

        Parameters
        ----------
        dataframe : `DataFrame`
            The DataFrame containing the seeing conditions data.
        addMostRecentBox : `bool`, optional
            Whether to add a box with the most recent seeing conditions.
        fig : `Figure`, optional
            The figure to plot on. If None, a new figure will be created.

        Returns
        -------
        fig : `Figure`
            The figure with the plotted seeing conditions.
        """
        ls = "-"
        ms = "o"
        df = dataframe

        if df.empty:
            raise ValueError("No data to plot for the given time range.")

        seeings = [SeeingConditions([row]) for _, row in df.iterrows()]

        if fig is None:
            fig = make_figure(figsize=(18, 10))
            ax1 = fig.add_subplot(111)
        else:
            fig.clear()
            ax1 = fig.add_subplot(111)

        utc = ZoneInfo("UTC")
        chile_tz = ZoneInfo("America/Santiago")

        # Function to convert UTC to Chilean time
        def offset_time_aware(utc_time: datetime) -> datetime:
            # Ensure the time is timezone-aware in UTC
            if utc_time.tzinfo is None:
                utc_time = utc_time.replace(tzinfo=utc)
            return utc_time.astimezone(chile_tz)

        df.index = pd.DatetimeIndex([t for t in df.index])

        ax1.plot([seeing.fwhmFree for seeing in seeings], "b", label='Free atmos. seeing"', ls=ls, marker=ms)
        ax1.plot([seeing.seeing for seeing in seeings], "r", label='Seeing "', ls=ls, marker=ms)
        ax1.plot(
            [seeing.seeing2 for seeing in seeings], "g", label='Profile adjusted seeing "', ls=ls, marker=ms
        )

        ax2 = ax1.twiny()
        ax2.set_xlim(ax1.get_xlim())

        # Format both axes to show only time
        ax1.xaxis.set_major_formatter(DateFormatter("%H:%M:%S"))
        ax2.xaxis.set_major_formatter(DateFormatter("%H:%M:%S"))

        # Apply the timezone-aware offset to the top axis ticks
        ax2.set_xticks(ax1.get_xticks())
        offset_ticks = [offset_time_aware(num2date(tick)) for tick in ax1.get_xticks()]
        ax2.set_xticklabels([tick.strftime("%H:%M:%S") for tick in offset_ticks])

        ax1.set_ylim(0, 1.1 * max([s.seeing2 for s in seeings]))
        ax1.set_xlabel("Time (UTC)")
        ax2.set_xlabel("Time (Chilean Time)")
        ax1.set_ylabel("Seeing (arcsec)")
        ax1.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5, alpha=0.7)
        ax1.xaxis.grid(True, which="major", linestyle="--", linewidth=0.5, alpha=0.7)

        # Update legend with larger font size
        ax1.legend(loc="lower left", fontsize=14)

        # Calculate current seeing and age of data
        if addMostRecentBox:
            currentSeeing = seeings[-1]
            justTime = currentSeeing.timestamp.isot.split("T")[1].split(".")[0]

            text = f'Current Seeing: {currentSeeing.seeing:.2f}"\n' f"Last updated @ {justTime} UTC"
            ax1.text(
                0.05,
                0.95,
                text,
                transform=ax1.transAxes,
                fontsize=14,
                verticalalignment="top",
                horizontalalignment="left",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.7),
            )

        return fig
