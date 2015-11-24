'''
psl2wiggle_stats.py - intersect psl and wiggle files to compute stats
=====================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python psl2wiggle_stats.py --help

Type::

   python psl2wiggle_stats.py --help

for command line help.

Command line options
--------------------

'''

import sys
import glob

import CGAT.Experiment as E
import CGAT.Stats as Stats
import CGAT.Blat as Blat

import CGAT.Wiggle as Wiggle
import alignlib_lite


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: psl2wiggle_stats.py 2781 2009-09-10 11:33:14Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome.")

    parser.add_option("--wiggle-file", dest="wiggle_files", type="string",
                      help="glob expression for wiggle files [%default].")

    parser.add_option("--column-prefix", dest="prefix", type="string",
                      help="prefix to add to contig names before lookup [%default].")

    parser.add_option("-z", "--from-zipped", dest="from_zipped", action="store_true",
                      help="input is zipped.")

    parser.add_option("--test", dest="test", type="int",
                      help="test - stop after # rows of parsing [%default].")

    parser.add_option("--with-values", dest="with_values", action="store_true",
                      help="output values in last column [%default].")

    parser.set_defaults(wiggle_files="*.data.bz2",
                        from_zipped=False,
                        prefix="",
                        with_values=False,
                        test=None)

    (options, args) = E.Start(parser, add_pipe_options=True)

    # open indexed access to wiggles
    wiggle_files = glob.glob(options.wiggle_files)
    if not wiggle_files:
        raise IOError("could not find wiggle files with '%s'" %
                      options.wiggle_files)

    index = Wiggle.WiggleMultiIndexedAccess(wiggle_files,
                                            keep_open=True,
                                            use_cache=False)

    iterator = Blat.BlatIterator(sys.stdin)

    ninput, noutput, nskipped = 0, 0, 0

    options.stdout.write(
        "query\tnali\t%s" % ("\t".join(Stats.DistributionalParameters().getHeaders())))
    if options.with_values:
        options.stdout.write("\tvalues")
    options.stdout.write("\n")

    while 1:

        if options.test and ninput >= options.test:
            break

        match = iterator.next()

        if match is None:
            break

        ninput += 1

        if options.loglevel >= 2:
            options.stdlog.write(str(match) + "\n")

        # psl always matches on the forward strand

        map_genome2query = alignlib_lite.py_makeAlignmentBlocks()
        f = alignlib_lite.py_AlignmentFormatBlat("%i\t%i\t%i\t%i\t%s\t%s\t%s\n" % (
            match.mSbjctFrom,
            match.mSbjctTo,
            match.mQueryFrom,
            match.mQueryTo,
            match.mSbjctBlockStarts,
            match.mQueryBlockStarts,
            match.mBlockSizes))
        f.copy(map_genome2query)

        data = index.get(options.prefix + match.mSbjctId,
                         match.mSbjctFrom,
                         match.mSbjctTo)

        values = []
        for x, vv in data:
            for v in vv:
                if map_genome2query.mapRowToCol(x) >= 0:
                    values.append(v)
                x += 1
        if len(values) == 0:
            nskipped += 1
            continue

        noutput += 1

        if options.loglevel >= 2:
            options.stdlog.write(
                "# %s\n" % ",".join(["%5.3f" % v for v in values]))

        s = Stats.DistributionalParameters(values)
        options.stdout.write("%s\t%i\t%s" % (match.mQueryId,
                                             match.mNMismatches +
                                             match.mNMatches,
                                             str(s)))

        if options.with_values:
            options.stdout.write(
                "\t%s" % ",".join(["%5.3f" % v for v in values]))

        options.stdout.write("\n")

    if options.loglevel >= 1:
        options.stdlog.write(
            "# ninput=%i, noutput=%i, nskipped=%i\n" % (ninput, noutput, nskipped))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
