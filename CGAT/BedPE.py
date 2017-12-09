'''bedPE.py - Tools for working with bedPEPE files (https://bedPEtools.readthedocs.io/en/latest/content/general-usage.html)
=========================================

This module contains methods for working with :term:`bedPE`
formatted files.

.. note::
   Another way to access the information in :term:`bedPE` formatted
   files is through pysam_.

The principal class is :class:`bedPE` to represent :term:`bedPE` formatted
entries.  The method :func:`iterate` iterates over a bedPE file and is
aware of UCSC track information that might be embeded in the
file.

Reference
---------

'''
import re
import numpy
import bisect
import itertools

from CGAT import NCL as ncl
from CGAT import IOTools as IOTools

Headers = [
    "contig1", "start1", "end1",
    "contig2", "start2", "end2",
    "name", "score", "strand1",
    "strand2",
    "additional_info"]


class bedPE(object):
    """an interval in bedPE format.

    Coordinates are represented as 0-based, half-open intervals.

    Fields in the record can be accessed as attributes or through
    a dictionary type access::

       print b.contig()
       print b["contig"]

    bedPE-formatted records can have a variable number of columuns
    with a minimum of 3. Accessing an optional attribute that is not present
    will raise an IndexError.

    Attributes
    ----------
    contig1 : string
       Chromosome/contig of the first interval.
    start1 : int
       Start position of the first interval.
    end1 : int
       End position of the first interval.
    contig2 : string
       Chromosome/contig of the second interval.
    start2 : int
       Start position of the second interval.
    end2 : int
       End position of the second interval.
    name : string
       Name of the whole interval (optional).
    score : string
       Score associated with interval (optional).
    strand1 : char
       Strand of the first interval (optional).
    strand2 : char
       Strand of the second interval (optional).
    additional_info : list
       Tab-separated list of additional, user-defined fields (optional).
    """

    map_key2field = {'name': 0,
                     'score': 1,
                     'strand1': 2,
                     'strand2': 3,
                     'additional_info': 4}

    default_value = "."

    def __init__(self):
        self.contig1 = None
        self.start1 = 0
        self.end1 = 0
        self.contig2 = None
        self.start2 = 0
        self.end2 = 0
        self.fields = [] # Optional fields
        self.track = None

    def __str__(self):
        return "\t".join((self.contig1, str(self.start1),
                          str(self.end1), (self.contig2),
                          str(self.start2), str(self.end2)) 
                         + tuple(map(str, self.fields)))

    def copy(self):
        '''Returns a new bedPE object that is a copy of this one'''

        new_entry = bedPE()
        new_entry.__dict__ = self.__dict__.copy()
        return new_entry



    def __contains__(self, key):
        return self.map_key2field[key] < len(self.fields)

    def __getitem__(self, key):
        return self.fields[self.map_key2field[key]]

    def __setitem__(self, key, value):
        try:
            position = self.map_key2field[key]
        except IndexError:
            raise IndexError("Unknown key: %s" % key)

        try:
            self.fields[position] = value
        except IndexError:

            self.fields.extend([self.default_value]
                               * (position - len(self.fields) + 1))

            self.fields[position] = value

    def __getattr__(self, key):
        try:
            return self.fields[self.map_key2field[key]]
        except IndexError:
            return None

    @property
    def columns(self):
        '''return number of columns in bedPE-entry.'''
        return 6 + len(self.fields)


class Track(object):
    """bedPE track information."""

    def __init__(self, line):
        r = re.compile('([^ =]+) *= *("[^"]*"|[^ ]*)')

        self._d = {}
        for k, v in r.findall(line[:-1]):
            if v[:1] == '"':
                self._d[k] = v[1:-1]
            else:
                self._d[k] = v

        self._line = line[:-1]

    def __str__(self):
        return self._line

    def __getitem__(self, key):
        return self._d[key]

    def __setitem__(self, key, val):
        self._d[key] = val


def iterator(infile):
    """iterate over a :term:`bedPE` formatted file.

    Comments and empty lines are ignored. The iterator is
    :term:`track` aware and will set the ``track`` attribute for the
    bedPE objects it yields.

    Arguments
    ---------
    infile : File

    Yields
    ------
    bedPE
       :class:`bedPE` object

    """

    track = None
    for line in infile:
        if line.startswith("track"):
            track = Track(line)
            continue
        # ignore comments
        if line.startswith("#"):
            continue
        # ignore empty lines (in order to parse pseudo bedPE files)
        if line.strip() == "":
            continue

        b = bedPE()
        # split at tab (bedPE standard, do not split at space as this will split
        # the name field)
        data = line[:-1].split("\t")
        try:
            b.contig1, b.start1, b.end1, b.contig2, b.start2, b.end2 = data[0], int(data[1]), int(data[2]), data[3], int(data[4]), int(data[5])
        except IndexError:
            raise ValueError("parsing error in line '%s'" % line[:-1])
        b.fields = data[6:]
        b.track = track
        yield b


def bedPE_iterator(infile):
    """Deprecated, use :func:`iterator`."""
    return iterator(infile)


def setName(iterator):
    """yield bedPE entries in which name is set to the record number if
    unset.

    Yields
    ------
    bedPE
       :class:`bedPE` object
    """
    for i, bedPE in enumerate(iterator):
        if "name" not in bedPE:
            bedPE.name = str(i)
        yield bedPE


def grouped_iterator(iterator):
    """yield bedPE results grouped by track.

    Note that the iterator supplied needs to be sorted by the track
    attribute. This is usually the case in :term:`bedPE` formatted
    files.

    Yields
    ------
    bedPE
       :class:`bedPE` object

    """
    return itertools.groupby(iterator, lambda x: x.track)


def getNumColumns(filename):
    '''return number of fields in bedPE-file by looking at the first
    entry.

    Returns
    -------
    ncolumns : int
       The number of columns. If the file is empty, 0 is returned.
    '''
    with IOTools.openFile(filename) as inf:
        for line in inf:
            if line.startswith("#"):
                continue
            if line.startswith("track"):
                continue
            return len(line[:-1].split("\t"))
    return 0
