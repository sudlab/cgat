##########################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
##########################################################################
"""FastaIterator.py - Iteration over fasta files
================================================

This module provides a simple iterator of Fasta formatted files.  The
difference to the biopython iterator is that the iterators in this
module skip over comment lines starting with "#".

.. note::
   Another way to access the information in :term:`fasta` formatted
   files is through pysam_.

Reference
---------

"""
import subprocess
import os


class FastaRecord:
    """a :term:`fasta` record.

    Attributes
    ----------
    title: string
       the title of the sequence

    sequence: string
       the sequence

    fold : int
       the number of bases per line when writing out
    """

    def __init__(self, title, sequence, fold=False):

        self.title = title
        self.sequence = sequence
        self.fold = fold

    def __str__(self):
        ''' str method for writing out'''

        if self.fold:
            seq = [self.sequence[i:i+self.fold]
                   for i in range(0, len(self.sequence), self.fold)]
        else:
            seq = (self.sequence,)

        return ">%s\n%s" % (self.title, "\n".join(seq))


class FastaIterator:
    '''a iterator of :term:`fasta` formatted files.

    Yields
    ------
    FastaRecord

    '''

    def __init__(self, f, *args, **kwargs):
        self.iterator = iterate(f)

    def __iter__(self):
        return self

    def next(self):
        return self.iterator.next()


def iterate(infile, comment="#", fold=False):
    '''iterate over fasta data in infile

    Lines before the first fasta record are
    ignored (starting with ``>``) as well as
    lines starting with the comment character.

    Parameters
    ----------
    infile : File
        the input file
    comment : char
        comment character
    fold : int
        the number of bases before line split when writing out

    Yields
    ------
    FastaRecord
    '''

    h = infile.readline()[:-1]

    if not h:
        raise StopIteration

    # skip everything until first fasta entry starts
    while h[0] != ">":
        h = infile.readline()[:-1]
        if not h:
            raise StopIteration
        continue

    h = h[1:]
    seq = []

    for line in infile:

        if line.startswith(comment):
            continue

        if line.startswith('>'):
            yield FastaRecord(h, ''.join(seq), fold)

            h = line[1:-1]
            seq = []
            continue

        seq.append(line[:-1])

    yield FastaRecord(h, ''.join(seq), fold)


def iterate_together(*args):
    """iterate synchronously over one or more fasta files.

    The iteration finishes once any of the files is exhausted.

    Arguments
    ---------

    :term:`fasta`-formatted files to be iterated upon

    Yields
    ------
    tuple
       a tuple of :class:`FastaRecord` corresponding to
       the current record in each file.
    """

    iterators = [FastaIterator(x) for x in args]

    while 1:
        yield [x.next() for x in iterators]


def count(filename):
    '''count number of sequences in fasta file.

    This method uses the ``grep`` utility to count
    lines starting with ``>``.

    Arguments
    ---------
    filename : string
        The filename

    Raises
    ------
    OSError
        If the file does not exist

    Returns
    -------
    int
        The number of sequences in the file.
    '''

    if filename.endswith(".gz"):
        statement = "zcat %s | grep -c '>'" % filename
    else:
        statement = "cat %s | grep -c '>'" % filename

    if not os.path.exists(filename):
        raise OSError("file '%s' does not exist" % filename)

    # grep returns error if no match is found
    try:
        return subprocess.check_output(statement, shell=True)
    except subprocess.CalledProcessError:
        return 0

class FastaWindow:

    def __init__(self, name, seq, offset):

        self.name = name
        self.sequence = seq
        self.offset = offset


class BufferedIterator:
    ''' For dealing with Fasta files whose records are too big to fit in memory.
    returns a record in chunks. The size of chunks is determined by size. 
    The sizes are only approximate: sequence is added a line at a time .If
    overlap is set, then the chunks overlap. Each iteration returns a FastaWindow
    object with the sequence, the name it came from and the offset from zero
    to count the position from '''

    def iterate(self,infile):

        h = infile.readline()[:-1]

        if not h:
            raise StopIteration

        # skip everything until first fasta entry starts
        while h[0] != ">":
            h = infile.readline()[:-1]
            if not h:
                raise StopIteration
            continue

        h = h[1:]
        seq = ""
        offset = 0

        for line in infile:
            
            if line.startswith(comment):
                continue

            if line.startswith('>'):
                yield FastaWindow(h, seq, offset)
                h = line[1:-1]
                offset = 0
                seq = ""

            seq += line[:-1]
            
            if len(seq) >= self.size:
                yield FastaWindow(h, seq, offset)
                offset = offset + len(seq) - self.overlap
                seq = seq[-self.overlap:]
                             
        yield FastaWindow(h, seq, offset)

    def __init__(self, infile, size, overlap):
        self.size = size
        self.overlap = overlap
        self.mIterator = self.iterate(infile)

    def __iter__(self):

        return self

    def __next__(self):

        return self.mIterator.next()

