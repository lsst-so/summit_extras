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

from lsst.summit.extras.imageSorter import TAGS, ImageSorter


def _idTrans(dataIdDictOrTuple: dict | tuple[int, int]) -> tuple[int, int]:
    """Convert a dataId to the internal ``(dayObs, seqNum)`` tuple form.

    Parameters
    ----------
    dataIdDictOrTuple : `dict` or `tuple` [`int`, `int`]
        A dataId expressed either as a dict with ``dayObs`` and ``seqNum``
        keys or as a ``(dayObs, seqNum)`` tuple.

    Returns
    -------
    dataId : `tuple` [`int`, `int`]
        The dataId as a ``(dayObs, seqNum)`` tuple.

    Raises
    ------
    RuntimeError
        Raised if ``dataIdDictOrTuple`` is neither a `dict` nor a `tuple`.
    """
    if isinstance(dataIdDictOrTuple, tuple):
        return dataIdDictOrTuple
    elif isinstance(dataIdDictOrTuple, dict):
        return (dataIdDictOrTuple["dayObs"], dataIdDictOrTuple["seqNum"])
    else:
        raise RuntimeError(f"Failed to parse dataId {dataIdDictOrTuple}")


class Annotations:
    """Interface for reading annotations written by `ImageSorter`.

    Loads the tag and note dictionaries from an annotations file and
    provides lookup helpers keyed by dataId.

    Parameters
    ----------
    filename : `str`
        Path to the annotations file produced by `ImageSorter`.
    """

    def __init__(self, filename: str):
        self.filename = filename
        self.tags, self.notes = self._load(filename)

    def _load(self, filename: str) -> tuple[dict, dict]:
        """Load tags and notes from the specified annotations file.

        Parameters
        ----------
        filename : `str`
            Path to the annotations file.

        Returns
        -------
        tags : `dict`
            Mapping of dataId tuples to tag strings.
        notes : `dict`
            Mapping of dataId tuples to note strings.
        """
        tags, notes = ImageSorter.loadAnnotations(filename)
        return tags, notes

    def getTags(self, dataId: dict | tuple[int, int]) -> str | None:
        """Get the tags for the specified dataId.

        Parameters
        ----------
        dataId : `dict` or `tuple` [`int`, `int`]
            The dataId to look up.

        Returns
        -------
        tags : `str` or `None`
            The tag string for the dataId. An empty string means the image
            was examined but no tags were set; `None` means the image was
            not examined.
        """
        return self.tags.get(_idTrans(dataId), None)

    def getNotes(self, dataId: dict | tuple[int, int]) -> str | None:
        """Get the notes for the specified dataId.

        Parameters
        ----------
        dataId : `dict` or `tuple` [`int`, `int`]
            The dataId to look up.

        Returns
        -------
        notes : `str` or `None`
            The note string for the dataId, or `None` if no notes exist.
        """
        return self.notes.get(_idTrans(dataId), None)

    def hasTags(self, dataId: dict | tuple[int, int], flags: str) -> bool | None:
        """Check whether a dataId has all of the specified tags.

        Parameters
        ----------
        dataId : `dict` or `tuple` [`int`, `int`]
            The dataId to look up.
        flags : `str`
            String of single-character tag flags to test for. The
            comparison is case-insensitive.

        Returns
        -------
        hasAll : `bool` or `None`
            `True` if all requested tags are present, `False` if any are
            missing, or `None` if the dataId has not been examined.
        """
        tag = self.getTags(dataId)
        if tag is None:  # not just 'if tag' because '' is not the same as None but both are falsy
            return None
        return all(i in tag for i in flags.upper())

    def getListOfCheckedData(self) -> list:
        """Return the sorted list of all examined dataIds.

        Returns
        -------
        dataIds : `list` [`tuple` [`int`, `int`]]
            Sorted list of ``(dayObs, seqNum)`` tuples that have been
            examined.
        """
        return sorted(list(self.tags.keys()))

    def getListOfDataWithNotes(self) -> list:
        """Return the sorted list of all dataIds that have notes.

        Returns
        -------
        dataIds : `list` [`tuple` [`int`, `int`]]
            Sorted list of ``(dayObs, seqNum)`` tuples with notes
            attached.
        """
        return sorted(list(self.notes.keys()))

    def isExamined(self, dataId: dict) -> bool:
        """Check whether the dataId has been examined.

        Parameters
        ----------
        dataId : `dict` or `tuple` [`int`, `int`]
            The dataId to look up.

        Returns
        -------
        examined : `bool`
            `True` if the dataId has an entry in the tag dictionary.
        """
        return _idTrans(dataId) in self.tags

    def printTags(self) -> None:
        """Print the list of tag definitions used by `ImageSorter`."""
        print(TAGS)

    def getIdsWithGivenTags(self, tags: str, exactMatches: bool = False) -> list:
        """Return dataIds that match the specified tag string.

        Parameters
        ----------
        tags : `str`
            String of single-character tag flags to match. The comparison
            is case-insensitive.
        exactMatches : `bool`, optional
            If `True`, only return dataIds whose tag string contains
            exactly the requested tags and no others. If `False`
            (default), return all dataIds whose tags are a superset of
            the requested tags.

        Returns
        -------
        dataIds : `list` [`tuple` [`int`, `int`]]
            List of matching ``(dayObs, seqNum)`` tuples.
        """
        if exactMatches:
            wanted = set(tags.upper())
            return [dId for (dId, tag) in self.tags.items() if set(tag) == wanted]
        else:
            return [dId for (dId, tag) in self.tags.items() if all(t in tag for t in tags.upper())]
