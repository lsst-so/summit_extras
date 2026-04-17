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

import filecmp
import hashlib
import logging
import os
import pickle
import sys
from typing import Any

import astropy
import numpy as np
from astropy.io import fits

# redirect logger to stdout so that logger messages appear in notebooks too
logging.basicConfig(level=logging.INFO, handlers=[logging.StreamHandler(sys.stdout)])
logger = logging.getLogger("headerFunctions")


def loadHeaderDictsFromLibrary(libraryFilename: str) -> tuple[dict, dict]:
    """Load cached header and hash dicts from a pickle library file.

    Parameters
    ----------
    libraryFilename : `str`
        Path of the library pickle file to load from.

    Returns
    -------
    headersDict : `dict`
        A dict, keyed by filename, with values being the full primary
        header, in the same format built by `buildHashAndHeaderDicts`.
        Empty if the file does not exist or fails to load.
    dataDict : `dict`
        A dict, keyed by filename, with values being hashes of the data
        sections, in the same format built by `buildHashAndHeaderDicts`.
        Empty if the file does not exist or fails to load.
    """
    try:
        with open(libraryFilename, "rb") as pickleFile:
            headersDict, dataDict = pickle.load(pickleFile)

        if len(headersDict) != len(dataDict):
            print("Loaded differing numbers of entries for the header and data dicts.")
            print(f"{len(headersDict)} vs {len(dataDict)}")
            print("Something has gone badly wrong - your library seems corrupted!")
        else:
            print(f"Loaded {len(headersDict)} values from pickle files")
    except Exception as e:
        if not os.path.exists(libraryFilename):
            print(
                f"{libraryFilename} not found. If building the header dicts for the first time this"
                " is to be expected.\nOtherwise you've misspecified the path to you library!"
            )
        else:
            print(f"Something more sinister went wrong loading headers from {libraryFilename}:\n{e}")
        return {}, {}

    return headersDict, dataDict


def _saveToLibrary(libraryFilename: str, headersDict: dict, dataDict: dict) -> None:
    """Persist header and data-hash dicts to a pickle library file.

    Parameters
    ----------
    libraryFilename : `str`
        Path of the library pickle file to write.
    headersDict : `dict`
        Dict of primary headers keyed by filename.
    dataDict : `dict`
        Dict of data-section hashes keyed by filename.
    """
    try:
        with open(libraryFilename, "wb") as dumpFile:
            pickle.dump((headersDict, dataDict), dumpFile, pickle.HIGHEST_PROTOCOL)
    except Exception:
        print("Failed to write pickle file! Here's a debugger so you don't lose all your work:")
        import ipdb as pdb

        pdb.set_trace()


def _findKeyForValue(
    dictionary: dict, value: Any, warnOnCollision: bool = True, returnCollisions: bool = False
) -> Any:
    """Return the key or keys that map to a given value.

    Parameters
    ----------
    dictionary : `dict`
        The dictionary to search.
    value : `object`
        The value to look up.
    warnOnCollision : `bool`, optional
        If `True`, log a warning when the value is held by more than one
        key.
    returnCollisions : `bool`, optional
        If `True`, return the full list of matching keys. Otherwise
        return only the first match.

    Returns
    -------
    keys : `object` or `list`
        A single matching key, or the full list of matching keys if
        ``returnCollisions`` is `True`.
    """
    listOfKeys = [k for (k, v) in dictionary.items() if v == value]
    if warnOnCollision and len(listOfKeys) > 1:
        logger.warning("Found multiple keys for value! Returning only first.")
    if returnCollisions:
        return listOfKeys
    if not listOfKeys:
        return None
    return listOfKeys[0]


def _hashFile(fileToHash: Any, dataHdu: int | str, sliceToUse: slice) -> str:
    """Return a hash of a 2D region of a FITS HDU's data array.

    Kept as a small helper so that hashing multiple HDUs (e.g. to cope
    with HDUs filled with zeros) can be added straightforwardly.

    Parameters
    ----------
    fileToHash : `astropy.io.fits.HDUList`
        The opened FITS file.
    dataHdu : `int` or `str`
        Index or extension name of the HDU holding the pixel data.
    sliceToUse : `slice`
        Slice applied to both axes of the data array before hashing.

    Returns
    -------
    hashStr : `str`
        Hex-encoded SHA-256 digest of the sliced data.
    """
    # Cast to a fixed dtype so an all-zero array hashes to ZERO_HASH
    # regardless of the HDU's native pixel dtype.
    data = fileToHash[dataHdu].data[sliceToUse, sliceToUse].astype(np.int32).tobytes()
    h = _hashData(data)
    return h


def _hashData(data: np.ndarray) -> str:
    """Return a hex SHA-256 hash of the given array's bytes.

    Parameters
    ----------
    data : `numpy.ndarray` or `bytes`
        The data to hash.

    Returns
    -------
    hashStr : `str`
        Hex-encoded SHA-256 digest, suitable for storing in a dict.
    """
    h = hashlib.sha256(data).hexdigest()  # hex because we want it readable in the dict
    return h


ZERO_HASH = _hashData(np.zeros((100, 100), dtype=np.int32))


def buildHashAndHeaderDicts(
    fileList: list[str], dataHdu: int | str = "Segment00", libraryLocation: str | None = None
) -> tuple[dict, dict]:
    """Build dicts of primary headers and data-section hashes.

    Data is hashed using the hard-coded 100x100 corner of the pixel
    array, i.e. ``file[dataHdu].data[0:100, 0:100]``. If a library
    location is supplied, existing entries are loaded from there and
    newly computed entries are written back.

    Parameters
    ----------
    fileList : `list` [`str`]
        The fully-specified paths of the files to scrape.
    dataHdu : `int` or `str`, optional
        The HDU to use for the pixel data to hash.
    libraryLocation : `str`, optional
        Path to a pickle file used to cache results across runs. If
        provided, existing entries are loaded before processing and
        newly computed entries are written back.

    Returns
    -------
    headersDict : `dict`
        A dict, keyed by filename, with values being the full primary
        header.
    dataDict : `dict`
        A dict, keyed by filename, with values being hashes of the
        file's data section.
    """
    headersDict: dict[str, Any] = {}
    dataDict: dict[str, str] = {}

    if libraryLocation:
        headersDict, dataDict = loadHeaderDictsFromLibrary(libraryLocation)

    # don't load files we already know about from the library
    filesToLoad = [f for f in fileList if f not in headersDict.keys()]

    s = slice(0, 100)
    for filenum, filename in enumerate(filesToLoad):
        if len(filesToLoad) > 1000 and filenum % 1000 == 0:
            if libraryLocation:
                logger.info(f"Processed {filenum} of {len(filesToLoad)} files not loaded from library...")
            else:
                logger.info(f"Processed {filenum} of {len(fileList)} files...")
        with fits.open(filename) as f:
            try:
                headersDict[filename] = f[0].header
                h = _hashFile(f, dataHdu, s)
                if h in dataDict.values():
                    collision = _findKeyForValue(dataDict, h, warnOnCollision=False)
                    logger.warning(
                        f"Duplicate file (or hash collision!) for files {filename} and " f"{collision}!"
                    )
                    if filecmp.cmp(filename, collision):
                        logger.warning("Filecmp shows files are identical")
                    else:
                        logger.warning(
                            "Filecmp shows files differ - "
                            "likely just zeros for data (or a genuine hash collision!)"
                        )

                dataDict[filename] = h
            except Exception:
                logger.warning(f"Failed to load {filename} - file is likely corrupted.")

    # we have always added to this, so save it back over the original
    if libraryLocation and len(filesToLoad) > 0:
        _saveToLibrary(libraryLocation, headersDict, dataDict)

    # have to pare these down, as library loaded could be a superset
    headersDict = {k: headersDict[k] for k in fileList if k in headersDict.keys()}
    dataDict = {k: dataDict[k] for k in fileList if k in dataDict.keys()}

    return headersDict, dataDict


def sorted(inlist: set | list, replacementValue: str = "<BLANK VALUE>") -> list:
    """Sort a list, coercing all values to strings and handling blanks.

    Replaces any ``astropy.io.fits.card.Undefined`` entries with
    ``replacementValue`` so they sort consistently alongside mixed
    string and integer values. This deliberately shadows the built-in
    `sorted` in this module's namespace.

    Parameters
    ----------
    inlist : `set` or `list`
        Values to sort.
    replacementValue : `str`, optional
        String used in place of undefined FITS card values.

    Returns
    -------
    output : `list` [`str`]
        Sorted list of string values.
    """
    from builtins import sorted as _sorted

    output = [
        str(x) if not isinstance(x, astropy.io.fits.card.Undefined) else replacementValue for x in inlist
    ]
    output = _sorted(output)
    return output


def keyValuesSetFromFiles(
    fileList: list[str],
    keys: list[str],
    joinKeys: list[str],
    noWarn: bool = False,
    printResults: bool = True,
    libraryLocation: str | None = None,
    printPerFile: bool = False,
) -> Any:
    """Get the set of values seen for a set of header keys across files.

    Parameters
    ----------
    fileList : `list` [`str`]
        The fully-specified paths of the files to scrape.
    keys : `list` [`str`]
        The header keys to scrape.
    joinKeys : `list` [`str`]
        List of keys to concatenate with ``+`` when scraping. For
        example, a header with ``FILTER1 = SDSS_u`` and
        ``FILTER2 = NB_640nm`` yields ``SDSS_u+NB_640nm``. Useful when
        looking for the actually observed combinations rather than the
        Cartesian product of the individual sets. Pass an empty list to
        skip join processing.
    noWarn : `bool`, optional
        If `True`, suppress warnings about keys missing from a header.
    printResults : `bool`, optional
        If `True`, print the collected values (and files with all-zero
        data) to stdout.
    libraryLocation : `str`, optional
        Path to a cached header library; passed through to
        `buildHashAndHeaderDicts`.
    printPerFile : `bool`, optional
        If `True`, print each ``filename<tab>key<tab>value`` triple as
        it is scraped. Prompts for confirmation when the output would
        exceed 200 lines.

    Returns
    -------
    kValues : `dict` [`str`, `set`] or `None`
        Mapping of each requested key to the set of values seen. `None`
        if ``keys`` is empty.
    joinedValues : `set` [`str`]
        Set of joined value strings. Only returned if ``joinKeys`` is
        non-empty.
    """
    print(f"Scraping headers from {len(fileList)} files...")
    if printPerFile and (len(fileList) * len(keys) > 200):
        print(f"You asked to print headers per-file, for {len(fileList)} files x {len(keys)} keys.")
        cont = input("Are you sure? Press y to continue, anything else to quit:")
        if not cont or cont.lower()[0] != "y":
            exit()

    headerDict, hashDict = buildHashAndHeaderDicts(fileList, libraryLocation=libraryLocation)

    kValues: dict[str, set[Any]] | None
    if keys:  # necessary so that -j works on its own
        kValues = {k: set() for k in keys}
    else:
        keys = []
        kValues = None

    joinedValues: set[str] = set()

    for filename in headerDict.keys():
        header = headerDict[filename]
        for key in keys:
            if key in header:
                assert kValues is not None
                kValues[key].add(header[key])
                if printPerFile:
                    print(f"{filename}\t{key}\t{header[key]}")
                    if len(keys) > 1 and key == keys[-1]:
                        # newline between files if multikey
                        print()
            else:
                if not noWarn:
                    logger.warning(f"{key} not found in header of {filename}")

        if joinKeys:
            jVals = None
            # Note that CCS doesn't leave values blank, it misses the whole
            # card out for things like FILTER2 when not being used
            jVals = [header[k] if k in header else "<missing card>" for k in joinKeys]

            # However, we do ALSO get blank cards to, so:
            # substitute <BLANK_VALUE> when there is an undefined card
            # because str(v) will give the address for each blank value
            # too, meaning each blank card looks like a different value
            joinedValues.add(
                "+".join(
                    [
                        str(v) if not isinstance(v, astropy.io.fits.card.Undefined) else "<BLANK_VALUE>"
                        for v in jVals
                    ]
                )
            )

    if printResults:
        # Do this first because it's messy
        zeroFiles = _findKeyForValue(hashDict, ZERO_HASH, warnOnCollision=False, returnCollisions=True)
        if zeroFiles:
            print("\nFiles with zeros for data:")
        for filename in zeroFiles:
            print(f"{filename}")

        if kValues is not None:
            for key in kValues.keys():
                print(f"\nValues found for header key {key}:")
                print(f"{sorted(kValues[key])}")

        if joinKeys:
            print(f"\nValues found when joining {joinKeys}:")
            print(f"{sorted(joinedValues)}")

    if joinKeys:
        return kValues, joinedValues

    return kValues


def compareHeaders(filename1: str, filename2: str) -> None:
    """Compare the headers of two FITS files in detail.

    First, the two files are confirmed to have the same pixel data (by
    hashing the first 100x100 pixels in HDU 1) to ensure the files
    really should be compared.

    The function then prints:

    - keys that appear in A but not in B
    - keys that appear in B but not in A
    - keys common to both files, split into those with identical
      values and those whose values differ (showing the differing
      values side-by-side)
    - whether the post-PDU HDU ``EXTNAME`` ordering differs between
      the two files

    Parameters
    ----------
    filename1 : `str`
        Full path to the first file to compare.
    filename2 : `str`
        Full path to the second file to compare.
    """
    assert isinstance(filename1, str)
    assert isinstance(filename2, str)

    headerDict1, hashDict1 = buildHashAndHeaderDicts([filename1])
    headerDict2, hashDict2 = buildHashAndHeaderDicts([filename2])

    if hashDict1[filename1] != hashDict2[filename2]:
        print("Pixel data was not the same - did you really mean to compare these files?")
        print(f"{filename1}\n{filename2}")
        cont = input("Press y to continue, anything else to quit:")
        if not cont or cont.lower()[0] != "y":
            exit()

    # you might think you don't want to always call sorted() on the key sets
    # BUT otherwise they seem to be returned in random order each time you run
    # and that can be crazy-making

    h1 = headerDict1[filename1]
    h2 = headerDict2[filename2]
    h1Keys = list(h1.keys())
    h2Keys = list(h2.keys())

    commonKeys = set(h1Keys)
    commonKeys = commonKeys.intersection(h2Keys)

    keysInh1NotInh2 = sorted([_ for _ in h1Keys if _ not in h2Keys])
    keysInh2NotInh1 = sorted([_ for _ in h2Keys if _ not in h1Keys])

    print(f"Keys in {filename1} not in {filename2}:\n{keysInh1NotInh2}\n")
    print(f"Keys in {filename2} not in {filename1}:\n{keysInh2NotInh1}\n")
    print(f"Keys in common:\n{sorted(commonKeys)}\n")

    # put in lists so we can output neatly rather than interleaving
    identical = []
    differing = []
    for key in commonKeys:
        if h1[key] == h2[key]:
            identical.append(key)
        else:
            differing.append(key)

    assert len(identical) + len(differing) == len(commonKeys)

    if len(identical) == len(commonKeys):
        print("All keys in common have identical values :)")
    else:
        print("Of the common keys, the following had identical values:")
        print(f"{sorted(identical)}\n")
        print("Common keys with differing values were:")
        for key in sorted(differing):
            d = "<blank card>".ljust(25)
            v1 = str(h1[key]).ljust(25) if not isinstance(h1[key], astropy.io.fits.card.Undefined) else d
            v2 = str(h2[key]).ljust(25) if not isinstance(h2[key], astropy.io.fits.card.Undefined) else d
            print(f"{key.ljust(8)}: {v1} vs {v2}")

    # Finally, check the extension naming has the same ordering.
    # We have to touch the files again, which is pretty lame
    # but not doing so would require the header builder to know about
    # file pairings or return extra info, and that's not ideal either,
    # and also not worth the hassle to optimise as this is only
    # ever for a single file, not bulk file processing
    numbering1, numbering2 = [], []
    with fits.open(filename1) as f1, fits.open(filename2) as f2:
        for hduF1, hduF2 in zip(f1[1:], f2[1:]):  # skip the PDU
            if "EXTNAME" in hduF1.header and "EXTNAME" in hduF2.header:
                numbering1.append(hduF1.header["EXTNAME"])
                numbering2.append(hduF2.header["EXTNAME"])

    if numbering1 != numbering2:
        print("\nSection numbering differs between files!")
        for s1, s2 in zip(numbering1, numbering2):
            print(f"{s1.ljust(12)} vs {s2.ljust(12)}")
    if len(numbering1) != len(numbering2):
        print("The length of those lists was also DIFFERENT! Presumably a non-image HDU was interleaved.")
