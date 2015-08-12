##########################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
'''
Expression.py - wrap various differential expression tools
===========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This module provides tools for differential expression analysis
for a variety of methods.

Methods implemented are:

   DESeq
   EdgeR
   ttest

The aim of this module is to run these individual tools and
output a table in a common format.

Usage
-----

Documentation
-------------

Requirements:

* DESeq >= 1.17
* DESeq2 >= 1.5.62
* edgeR >= 3.7.16
* gplots >= 2.14.2
* ggplot2 >= 1.0.0
* reshape >= 0.8.5
* RColorBrewer >= 1.0.5
* grid >= 3.1.1
* limma >= 3.21.18
* samr >= 2.0 (optional)
* siggenes >= 1.39.0 (optional)

Code
----

To do:
--check contrasts against design model

'''

import math
import numpy
import sys
import collections
import itertools
import re
import pandas
import ggplot
import copy
import numpy as np
from scipy.stats import ttest_ind
import matplotlib
import matplotlib.pyplot as plt
import pylab

import rpy2
from rpy2.robjects import r as R
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

try:
    import CGAT.Experiment as E
    import CGAT.IOTools as IOTools
    import CGAT.Stats as Stats
except ImportError:
    import Experiment as E
    import IOTools
    import Stats

# activate pandas/rpy conversion
pandas2ri.activate()

# AH: Only do this on demand, module might not be
#     be able to be imported if there are any issues.
# grdevices = importr('grDevices')


def runDETest(raw_DataFrame,
              design_file,
              outfile,
              de_caller,
              **kwargs):
    ''' provide higher level API to run tools with default setting '''
    if de_caller.lower() == "deseq":
        pass
    else:
        raise ValueError("Unknown caller")


def splitModel(model):
    '''returns the terms in the model'''
    return [x for x in
            re.split("[\.:,~+\\s*]", re.sub("~(\s*0\s*)?", "", model)) if
            len(x) > 0]


def adjustPvalues(p_values):
    '''return a list of BH adjusted pvalues'''
    # import r stats module to adjust pvalues
    stats = importr('stats')
    adj_pvalues = list(stats.p_adjust(FloatVector(p_values), method='BH'))
    return adj_pvalues


def pvaluesToSignficant(p_values, fdr):
    '''return a list of bools for significance'''
    return [int(x < fdr) for x in p_values]


def makeMAPlot(resultsTable, title, outfile):
    '''make an MA plot from a DE results table'''

    plotter = R('''
    function(results, title, outfile){
    suppressMessages(library(ggplot2))

    m_text = element_text(size = 15)

    results=results[order(results$significant),]
    results$significant = factor(results$significant, levels=c(0,1))

    p = ggplot(results[!is.na(results$p_value),],
        aes(x=log(((control_mean + treatment_mean)/2), 10),
            y = log(fold,2), col = significant))
    p = p + geom_point(aes(size=significant))
    p = p + scale_colour_manual(values=c("black", "red"),name="p-adjust < fdr")
    p = p + scale_size_manual(values=c(1,2), name="p-adjust < fdr")
    p = p + xlab("Normalised counts (log10)") + ylab("Fold change (log2)")
    p = p + ggtitle(title)
    p = p + theme(axis.title.x = m_text, axis.title.y = m_text,
                  axis.text.x = m_text, axis.text.y = m_text,
                  legend.title = m_text, legend.text = m_text)

    ggsave(plot = p, file = outfile, width=7, height=7)}''')

    r_resultsTable = pandas2ri.py2ri(resultsTable)

    plotter(r_resultsTable, title, outfile)


class ExperimentalDesign(object):
    """Objects representing experimental designs.

    This class takes an experimental design in tabular
    form and exports several convenience functions and
    attributes.

    `filename_or_table` can be a filename of a tab-separated table
    with the following columns.

    track
       the sample name

    include
       whether or not this sample should be included
       in the design
    
    group
       a label grouping several samples into a group
    
    pair
       for paired tests, a numeric identifier linking
       samples that are paired.
    
    An example of an experimtal design with two groups and paired
    samples is below::

        track   include group     pair
        sample1 1       treatment 1
        sample2 1       treatment 2
        sample3 1       control   1
        sample4 1       control   2

    When providing `filename_or_table`, the `include` column is used
    to directly filter the design to remove any non-included samples.

    Additional columns will be added as factors to the design.

    Alternatively, `filename_or_table` can be a pandas DataFrame with
    sample names as row index and the appropriate columns.

    Attributes
    -----------

    table : pandas DataFrame
       dataframe object describing the design
    group : list
       list of groups in the design
    conditions : list
       group for each sample
    pairs : list
       pair status for each sample
    samples: list
       sample names
    factors: list
       factors for each sample
    has_replicates : bool
       True if at least one group has multiple samples
    has_pairs : bool
       True if design is a paired design

    """

    def __init__(self, filename_or_table):
        # read in table in the constructor for ExpDesign
        # e.g design = ExpDesign(pd.read_csv(...))

        if isinstance(filename_or_table, str):
            self.table = pandas.read_csv(filename_or_table, sep="\t",
                                         index_col=0)
        else:
            self.table = filename_or_table

        assert self.table.shape, "design table is empty"

        # parse the design table. Users probably expect this
        # to happen once data is uploaded.
        self._update()

    def _update(self):
        """parse design file and fill class attributes.

        Call this function whenever self.table changes.
        """

        # remove all entries that should not be included
        self.table = self.table[self.table["include"] != 0]

        # define attributes
        self.conditions = self.table['group'].tolist()
        self.pairs = self.table['pair'].tolist()
        # TS - use OrderedDict to retain order in unique
        self.groups = (list(collections.OrderedDict.fromkeys(
            self.conditions)))
        self.samples = self.table.index.tolist()

        # TS, adapted from JJ code for DESeq2 design tables:
        # if additional columns present, pass to 'factors'
        if len(self.table.columns) > 3:
            self.factors = self.table.iloc[:, 3:]
        else:
            self.factors = None

        # Test if replicates exist, i.e at least one group has multiple samples
        # TS - does this need to be extended to check whether replicates exist
        # for each group?
        max_per_group = max([self.conditions.count(x) for x in self.groups])
        self.has_replicates = max_per_group >= 2

        # Test if pairs exist:
        npairs = len(set(self.pairs))
        has_pairs = npairs == 2

        # ..if so, at least two samples are required per pair
        if has_pairs:
            min_per_pair = min([self.pairs.count(x) for x in set(self.pairs)])
            self.has_pairs = min_per_pair >= 2
        else:
            self.has_pairs = False

    def validate(self, counts, model=None):
        
        missing = set(self.samples).difference(set(counts.table.columns))
        if len(missing) > 0:
            raise ValueError("following samples in design table are missing"
                             " from counts table: %s" % ", ".join(missing))

        if model is not None:

            model_terms = splitModel(model)

            # check all model terms exist
            missing = set(model_terms).difference(
                set(self.factors.columns.tolist()))

            if len(missing) > 0:
                raise ValueError("following terms in the model are missing"
                                 " from the design table: %s" %
                                 ", ".join(missing))

            # check there are at least two values for each level
            for term in model_terms:
                levels = set(self.factors.ix[:, term])
                if len(levels) < 2:
                    raise ValueError("term '%s' in the model has less "
                                     "than two levels (%s) in the "
                                     " design table" %
                                     (term, ", ".join(levels)))

    def restrict(self, counts):
        ''' return design with samples not in counts table removed '''

        self.table = self.table.ix[counts.table.columns, :]

    def revalidate(self, counts, model=None):
        ''' re-validate, i.e post filtering of counts table '''

        if len(set(self.samples).symmetric_difference(
                set(counts.table.columns))) > 0:
            self.restrict(counts)
            self._update()
            self.validate(counts, model)

        else:
            pass

    def firstPairOnly(self):
        '''restrict the design table to the first pair only.

        If unpaired will retain whole design table
        '''

        if not self.pairs:
            self.pairs = self.table['pair'].tolist()
        self.table = self.table.ix[self.table['pair'] == min(self.pairs), ]

    def getSamplesInGroup(self, group):
        """return list of sample names belonging to group."""
        if group not in self.groups:
            raise KeyError("group '%s' not present")
        return self.table[self.table["group"] == group].index.tolist()

    def getGroupForSample(self, sample):
        """return group a sample belongs to"""
        return self.table.loc[sample]["group"]

    def getGroups2Samples(self):
        """return a dictionary mapping a group to samples within the group.

        Returns
        -------
        dict
            with groups as keys and list of samples within a group as values.
        
        """
        groups_to_tracks = {}

        for group in self.groups:
            match_group = (self.table['group'] == group).tolist()
            subset = self.table.iloc[match_group, ]
            groups_to_tracks[group] = subset.index.tolist()
        return groups_to_tracks

    def mapGroupsSuffix(self, shuffle_suffix, keep_suffix):
        '''use suffixes supplied to extract groups from the
        design table and return dictionaries mapping each group to tracks
        for keeping with tracks which should be shuffled
        '''

        groups_to_keep_tracks = {}
        groups_to_spike_tracks = {}

        keep_suffix = keep_suffix.split(",")

        for group in self.groups:
            match_group = (self.table['group'] == group).tolist()
            tmp_design = self.table.iloc[match_group, ]

            groups_to_spike_tracks[group] = [
                x + shuffle_suffix for x in tmp_design.index.tolist()]

            groups_to_keep_tracks[group] = copy.copy(
                groups_to_spike_tracks[group])
            groups_to_keep_tracks[group].extend(
                [x + y for x in tmp_design.index.tolist() for y in keep_suffix])

        return groups_to_keep_tracks, groups_to_spike_tracks


class DEResult(object):
    ''' base class for DE result '''

    def __init__(self, testTable=None):
        self.table = testTable

    def getResults(self):
        ''' post-process results into generic output
        columns are:
        - contrast
        - treatment_name
        - control_name
        - test_id
        - control_mean
        - treatment_mean
        - control_std
        - treatment_std
        - p_value
        - p_value_adj
        - significant
        - l2fold
        - transformed_l2fold
        - fold
        - status
        '''
        pass

    def summariseDEResults(self):
        ''' summarise DE results. Counts instances of possible outcomes'''

        # TS: the summarising is now split by the comparison being made and a
        # dict returned with keys=comparisons, value=E.Counter per comparison

        self.Summary = {}

        control_names = set(self.table['control_name'])
        treatment_names = set(self.table['treatment_name'])

        for control, treatment in itertools.product(control_names,
                                                    treatment_names):

            tmp_table = self.table[self.table['control_name'] == control]
            tmp_table = tmp_table[tmp_table['treatment_name'] == treatment]
            tmp_table.reset_index(inplace=True)

            # check control, treatment combination exists
            n_rows = tmp_table.shape[0]
            if n_rows > 0:

                if control != treatment:
                    label = control + "_" + treatment
                else:
                    label = control
                label = re.sub(":", "_int_", label)

                counts = E.Counter()

                counts.signficant = sum(tmp_table['significant'])
                counts.insignficant = (len(tmp_table['significant']) -
                                       counts.signficant)
                counts.all_over = sum([x > 0 for x in tmp_table['l2fold']])
                counts.all_under = sum([x < 0 for x in tmp_table['l2fold']])
                counts.signficant_over = sum(
                    [tmp_table['significant'][x] == 1 and
                     tmp_table['l2fold'][x] > 1 for x in range(0, n_rows)])
                counts.signficant_under = sum(
                    [tmp_table['significant'][x] == 1 and
                     tmp_table['l2fold'][x] < 1 for x in range(0, n_rows)])

                self.Summary[label] = counts

    def plotMA(self, contrast=None, outfile_prefix=None):
        ''' base function for making a MA plot '''

        makeMAPlot = R('''
        function(df){

        suppressMessages(library(ggplot2))
        suppressMessages(library(grid))

        l_txt = element_text(size=20)

        tmp_df = df[df$contrast=="%(contrast)s",]
        p = ggplot(tmp_df, aes(log((control_mean+treatment_mean)/2,2),
                            -transformed_l2fold,
                            colour=as.factor(significant))) +

        geom_point() + xlab("log2 mean expression") + ylab("log2 fold change")+
        ggtitle("%(contrast)s") +
        scale_colour_manual(name="Significant", values=c("black", "#619CFF")) +
        guides(colour = guide_legend(override.aes = list(size=10)))+
        theme(axis.text.x = l_txt, axis.text.y = l_txt,
              axis.title.x = l_txt, axis.title.y = l_txt,
              legend.title = l_txt, legend.text = l_txt,
              title=l_txt, legend.key.size=unit(1, "cm"))

        suppressMessages(
          ggsave(file="%(outfile_prefix)s_%(contrast)s_MA_plot.png",
          width=10, height=10))}''' % locals())

        makeMAPlot(pandas2ri.py2ri(self.table))

    def plotVolcano(self, contrast=None, outfile_prefix=None):
        ''' base function for Volcano plotting'''

        makeVolcanoPlot = R('''
        function(df){

        suppressMessages(library(ggplot2))
        suppressMessages(library(grid))

        l_txt = element_text(size=20)

        tmp_df = df[df$contrast=="%(contrast)s",]

        p = ggplot(tmp_df, aes(transformed_l2fold, -log(p_value,10),
                   colour=as.factor(significant))) +
        geom_point() + xlab("log2 fold change") + ylab("p-value (-log10)") +
        ggtitle("%(contrast)s") +
        scale_colour_manual(name="Significant", values=c("black", "#619CFF")) +
        guides(colour = guide_legend(override.aes = list(size=10)))+
        theme(axis.text.x = l_txt, axis.text.y = l_txt,
              axis.title.x = l_txt, axis.title.y = l_txt,
              legend.title = l_txt, legend.text = l_txt,
              title=l_txt, legend.key.size=unit(1, "cm"))

        suppressMessages(
          ggsave(file="%(outfile_prefix)s_%(contrast)s_volcano_plot.png",
          width=10, height=10))}''' % locals())

        makeVolcanoPlot(pandas2ri.py2ri(self.table))


class DEExperiment(object):
    ''' base clase for DE experiments '''

    def __init__(self):
        pass

    def __call__(self):
        ''' call DE and generate an initial results table '''
        self.callDifferentialExpression()

    def run(self):
        ''' Custom DE functions '''


class DEExperiment_TTest(DEExperiment):
    '''DECaller object to run TTest on counts data'''

    # TS: to do: deal with genes/regions with zero counts

    def run(self, counts, design, normalise=True,
            normalise_method="deseq-size-factors"):

        # TS: normalisation performed here rather than earlier as
        # the method of normalisation is dependent upon the DE test
        if normalise is True:
            counts.normalise(method=normalise_method)

        df_dict = collections.defaultdict(list)

        for combination in itertools.combinations(design.groups, 2):

            control, treatment = combination
            n_rows = counts.table.shape[0]
            df_dict["control_name"].extend((control,)*n_rows)
            df_dict["treatment_name"].extend((treatment,)*n_rows)
            df_dict["test_id"].extend(counts.table.index.tolist())

            # set all status values to "OK"
            df_dict["status"].extend(("OK",)*n_rows)

            # subset counts table for each combination
            c_keep = [x == control for
                      x in design.conditions]
            control_counts = counts.table.iloc[:, c_keep]
            t_keep = [x == treatment for
                      x in design.conditions]
            treatment_counts = counts.table.iloc[:, t_keep]

            c_mean = control_counts.mean(axis=1)
            df_dict["control_mean"].extend(c_mean)
            df_dict["control_std"].extend(control_counts.std(axis=1))

            t_mean = treatment_counts.mean(axis=1)
            df_dict["treatment_mean"].extend(t_mean)
            df_dict["treatment_std"].extend(treatment_counts.std(axis=1))

            t, prob = ttest_ind(control_counts, treatment_counts, axis=1)
            df_dict["p_value"].extend(prob)

        result = DEResult_TTest(testTable=pandas.DataFrame(df_dict))
        result.table.set_index("test_id", inplace=True)

        return result


class DEResult_TTest(DEResult):

    def getResults(self, fdr):
        ''' post-process test results table into generic results output '''

        # TS - what about zero values?!
        self.table["fold"] = (
            self.table["treatment_mean"] / self.table["control_mean"])

        self.table["p_value_adj"] = adjustPvalues(self.table["p_value"])

        self.table["significant"] = pvaluesToSignficant(
            self.table["p_value_adj"], fdr)

        self.table["l2fold"] = list(numpy.log2(self.table["fold"]))

        # note: the transformed log2 fold change is not transformed for TTest
        self.table["transformed_l2fold"] = self.table["l2fold"]
        self.table["contrast"] = "_vs_".join((self.table['control_name'],
                                              self.table['treatment_name']))


class DEExperiment_edgeR(DEExperiment):
    '''DEExperiment object to run edgeR on counts data

    See page 13 of the EdgeR user guide::

       2. Simply pick a reasonable dispersion value, based on your
       experience with similar data, and use that. Although
       subjective, this is still more defensible than assuming Poisson
       variation. Typical values are dispersion=0.4 for human data,
       dispersion=0.1 for data on genetically identical model
       organisms or dispersion=0.01 for technical replicates.
    '''

    def run(self,
            counts,
            design,
            model=None,
            dispersion=None,
            ref_group=None,
            contrasts=None,
            outfile_prefix=None):

        if not design.has_replicates and dispersion is None:
            raise ValueError("no replicates and no dispersion")

        # create r objects
        r_counts = pandas2ri.py2ri(counts.table)
        r_groups = ro.StrVector(design.conditions)
        r_pairs = ro.StrVector(design.pairs)
        r_has_pairs = ro.default_py2ri(design.has_pairs)
        r_has_replicates = ro.default_py2ri(design.has_replicates)

        if dispersion is not None:
            r_dispersion = ro.default_py2ri(dispersion)
        else:
            r_dispersion = ro.default_py2ri(False)

        if model is not None:
            r_factors_df = pandas2ri.py2ri(design.factors)
        else:
            r_factors_df = ro.default_py2ri(False)

        if ref_group is not None:
            r_ref_group = ro.default_py2ri(ref_group)
        else:
            r_ref_group = ro.default_py2ri(design.groups[0])

        if contrasts is not None:
            raise ValueError("cannot currently handle user defined contrasts")
            r_contrasts = ro.default_py2ri(contrasts)

        E.info('running EdgeR: groups=%s, pairs=%s, replicates=%s, pairs=%s' %
               (design.groups, design.pairs, design.has_replicates,
                design.has_pairs))

        # build DGEList object
        buildDGEList = R('''
        suppressMessages(library('edgeR'))

        function(counts, groups, ref_group, factors_df){

        countsTable = DGEList(counts, group=groups)

        if (factors_df != FALSE){
          for (level in colnames(factors_df)){
            countsTable$samples[level] <- factors_df[level]
        }}

        countsTable$samples$group <- relevel(countsTable$samples$group,
        ref = ref_group)

        countsTable = calcNormFactors(countsTable)

        return(countsTable)}''')

        r_countsTable = buildDGEList(r_counts, r_groups,
                                     r_ref_group, r_factors_df)

        # build design matrix
        if model is None:
            buildDesign = R('''

            function(countsTable, has_pairs){

            if (has_pairs==TRUE) {
              design <- model.matrix( ~pairs + countsTable$samples$group ) }

            else {
              design <- model.matrix( ~countsTable$samples$group ) }

            return(design)
            }''')

            r_design = buildDesign(r_countsTable, r_has_pairs)
        else:
            buildDesign = R('''
            function(factors_df){
            design <- model.matrix(%s, data=factors_df)
            return(design)}''' % model)

            r_design = buildDesign(r_factors_df)

        # TS - for debugging, remove from final version
        E.info("design_table:")
        E.info(r_design)

        # fit model
        fitModel = R('''
        function(countsTable, design, has_replicates, dispersion){

        if (has_replicates == TRUE) {

            # estimate common dispersion
            countsTable = estimateGLMCommonDisp( countsTable, design )

            # estimate trended dispersion
            countsTable <- estimateGLMTrendedDisp( countsTable, design)

            # estimate tagwise dispersion
            countsTable = estimateGLMTagwiseDisp( countsTable, design )

            # fitting model to each tag
            fit = glmFit( countsTable, design ) }

        else {
            # fitting model to each tag
            fit = glmFit(countsTable, design, dispersion=dispersion) }

        return(fit)}''')

        r_fit = fitModel(r_countsTable, r_design,
                         r_has_replicates, r_dispersion)

        E.info("Conducting liklihood ratio tests")

        # TS - if no contrasts are specified, perform LR test on all possible
        # contrasts, otherwise, only perform the contrasts specified
        # TS - Function definition should depend on whether contrasts
        # are specified (keep the decision tree in python)
        # TS - To do:
        # add lrtTest definition for user-supplied contrasts
        if contrasts is None:
            lrtTest = R('''
        function(fit, prefix, countsTable, design){
        suppressMessages(library(reshape2))

        lrt_table_list = NULL

        for(coef in seq(2, length(colnames(design)))){
          lrt = glmLRT(fit, coef = coef)


          lrt_table = lrt$table
          # need to include observations as a seperate column as there will
          # be non-unique
          lrt_table$observation = rownames(lrt_table)
          rownames(lrt_table) <- NULL

          lrt_table_list[[coef]] = lrt_table
          lrt_table_list[[coef]]['contrast'] = colnames(design)[coef]

          dt <- decideTestsDGE(lrt)
          isDE <- as.logical(dt)
          DEnames <- rownames(fit)[isDE]

          contrast = gsub(":", "_interaction_", colnames(design)[coef])
          png(paste0(contrast, "MAplot.png"))
          plotSmear(lrt, de.tags=DEnames, cex=0.35, main=contrast)
          abline(h=c(-1,1), col="blue")
          dev.off()
        }

            lrt_final = do.call(rbind, lrt_table_list)

        return(lrt_final)}''')

            r_lrt_table = lrtTest(r_fit, outfile_prefix,
                                  r_countsTable, r_design)
        else:
            # TS - shouldn't get to here as error thrown earlier if
            # contrasts is not None
            pass

        E.info("Generating output - cpm table")

        # output cpm table
        outputCPMTable = R('''function(countsTable, outfile_prefix){
        suppressMessages(library(reshape2))
        countsTable.cpm <- cpm(countsTable, normalized.lib.sizes=TRUE)
        melted <- melt(countsTable.cpm)
        names(melted) <- c("test_id", "sample", "ncpm")

        # melt columns are factors - convert to string for sorting
        melted$test_id = levels(melted$test_id)[as.numeric(melted$test_id)]
        melted$sample = levels(melted$sample)[as.numeric(melted$sample)]

        # sort cpm table by test_id and sample
        sorted <- melted[with(melted, order(test_id, sample)),]
        gz <- gzfile(paste0(outfile_prefix,"cpm.tsv.gz"), "w" )
        write.table(sorted, file=gz, sep = "\t", row.names=FALSE, quote=FALSE)
        close(gz)}''')

        outputCPMTable(r_countsTable, outfile_prefix)

        result = DEResult_edgeR(testTable=pandas2ri.ri2py(r_lrt_table))

        return result


class DEResult_edgeR(DEResult):

    def getResults(self, fdr, DEtype="GLM"):
        ''' post-process test results table into generic results output '''

        E.info("Generating output - results table")

        df_dict = collections.defaultdict()

        n_rows = self.table.shape[0]

        df_dict["treatment_name"] = self.table['contrast']
        df_dict["control_name"] = self.table['contrast']
        df_dict["test_id"] = self.table['observation']
        df_dict["control_mean"] = self.table['logCPM']
        df_dict["treatment_mean"] = self.table['logCPM']
        df_dict["control_std"] = (0,)*n_rows
        df_dict["treatment_std"] = (0,)*n_rows
        df_dict["p_value"] = self.table['PValue']
        df_dict["p_value_adj"] = adjustPvalues(self.table['PValue'])
        df_dict["significant"] = pvaluesToSignficant(
            df_dict["p_value_adj"], fdr)
        df_dict["l2fold"] = list(numpy.log2(self.table['logFC']))

        # TS: the transformed log2 fold change is not transformed!
        df_dict["transformed_l2fold"] = df_dict["l2fold"]

        # TS: check what happens when no fold change is available
        # TS: may need an if/else in list comprehension. Raise E.warn too?
        df_dict["fold"] = [math.pow(2, float(x)) for
                           x in self.table['logFC']]

        # set all status values to "OK"
        # TS: again, may need an if/else to check...
        df_dict["status"] = ("OK",)*n_rows

        self.table = pandas.DataFrame(df_dict)
        self.table.set_index("test_id", inplace=True)

    def plotMAplot(self, design, outfile_prefix):
        # need to implement edgeR specific MA plot
        raise ValueError("MA plotting is not yet implemented for edgeR")


class DEExperiment_DESeq(DEExperiment):
    '''DEExperiment object to run DEseq on counts data
    '''

    # TS: this is a work in progress!

    def run(self,
            counts,
            design,
            model=None,
            dispersion_method="pooled",
            ref_group=None,
            contrasts=None,
            outfile_prefix=None,
            sharing_mode="maximum",
            fit_type="parametric"):

        # create r objects
        r_counts = pandas2ri.py2ri(counts.table)
        r_groups = ro.StrVector(design.conditions)
        r_pairs = ro.StrVector(design.pairs)
        r_has_pairs = ro.default_py2ri(design.has_pairs)
        r_has_replicates = ro.default_py2ri(design.has_replicates)
        r_dispersion_method = ro.default_py2ri(dispersion_method)
        r_sharing_mode = ro.default_py2ri(sharing_mode)
        r_fit_type = ro.default_py2ri(fit_type)

        if model is not None:
            r_factors_df = pandas2ri.py2ri(design.factors)
        else:
            r_factors_df = ro.default_py2ri(False)

        if ref_group is not None:
            r_ref_group = ro.default_py2ri(ref_group)
        else:
            r_ref_group = ro.default_py2ri(design.groups[0])

        if contrasts is not None:
            raise ValueError("cannot currently handle user defined contrasts")
            r_contrasts = ro.default_py2ri(contrasts)

        if dispersion_method == "pooled" and not design.has_replicates:
            E.warn("cannot use pooled dispersion without replicates, switching"
                   " to 'blind dispersion method")
            r_dispersion_method = ro.default_py2ri("blind")

        E.info('running DESeq: groups=%s, pairs=%s, replicates=%s, pairs=%s' %
               (design.groups, design.pairs, design.has_replicates,
                design.has_pairs))

        E.info("dispersion_method=%s, fit_type=%s, sharing_mode=%s" %
               (dispersion_method, fit_type, sharing_mode))

        # load DESeq
        R('''suppressMessages(library('DESeq'))''')

        buildCountDataSet = R('''
        function(counts, groups, dispersion_method, fit_type, sharing_mode){

        # create new counts data object for deseq
        cds <- newCountDataSet( counts, groups)

        # estimate size factors
        cds <- estimateSizeFactors(cds)
        print(sharing_mode)
        print(dispersion_method)
        print(fit_type)
        # estimate dispersion
        cds <- estimateDispersions(cds, method=dispersion_method,
               fitType=fit_type, sharingMode=sharing_mode)

        }''')

        buildCountDataSet(r_counts, r_groups,  r_dispersion_method,
                          r_fit_type, r_sharing_mode)

        result = None
        # result = DEResult_DESeq(testTable=pandas2ri.ri2py(r_lrt_table))

        return result


class DEResult_DEResult(DEResult):

    # re-write for DESeq
    def getResults(self, fdr, DEtype="pairwise"):
        ''' post-process test results table into generic results output '''

        E.info("Generating output - results table")

        df_dict = collections.defaultdict()

        n_rows = self.table.shape[0]

        if DEtype == "GLM":
            df_dict["treatment_name"] = self.table['contrast']
            df_dict["control_name"] = self.table['contrast']

        else:
            # TS: edgeR is currently only set up to run GLM-based tests
            pass

        df_dict["test_id"] = self.table['observation']
        df_dict["contrast"] = self.table['contrast']
        df_dict["control_mean"] = self.table['logCPM']
        df_dict["treatment_mean"] = self.table['logCPM']
        df_dict["control_std"] = (0,)*n_rows
        df_dict["treatment_std"] = (0,)*n_rows
        df_dict["p_value"] = self.table['PValue']
        df_dict["p_value_adj"] = adjustPvalues(self.table['PValue'])
        df_dict["significant"] = pvaluesToSignficant(
            df_dict["p_value_adj"], fdr)
        df_dict["l2fold"] = list(numpy.log2(self.table['logFC']))

        # TS: the transformed log2 fold change is not transformed!
        df_dict["transformed_l2fold"] = df_dict["l2fold"]

        # TS: check what happens when no fold change is available
        # TS: may need an if/else in list comprehension. Raise E.warn too?
        df_dict["fold"] = [math.pow(2, float(x)) for
                           x in self.table['logFC']]

        # set all status values to "OK"
        # TS: again, may need an if/else to check...
        df_dict["status"] = ("OK",)*n_rows

        self.table = pandas.DataFrame(df_dict)
        self.table.set_index("test_id", inplace=True)

    def plotMAplot(self, design, outfile_prefix):
        # need to implement edgeR specific MA plot
        raise ValueError("MA plotting is not yet implemented for edgeR")


class DEExperiment_DESeq2(DEExperiment):
    '''DEExperiment object to run DESeq2 on counts data'''

    def run(self,
            counts,
            design,
            model=None,
            contrasts=None,
            outfile_prefix=None,
            fdr=0.1):

        # create r objects
        r_counts = pandas2ri.py2ri(counts.table)
        r_groups = ro.StrVector(design.conditions)
        r_pairs = ro.StrVector(design.pairs)
        r_has_pairs = ro.default_py2ri(design.has_pairs)
        r_has_replicates = ro.default_py2ri(design.has_replicates)

        if design.factors is not None:
            r_factors_df = pandas2ri.py2ri(design.factors)
        else:
            r_factors_df = ro.default_py2ri(False)

        if contrasts is not None:
            DEtype = "GLM"

            # if model not included, use the column names from design.factors
            if not model:

                if design.factors is not None:
                    model = "~" + "+".join(design.factors.columns)
                    model_terms = design.factors.columns.values.tolist()

                else:
                    E.warn("need to supply a full model or else "
                           "additional columns in the design table "
                           "which will be taken as the full model")
            else:
                model_terms = [x for x in re.split("[\+~ ]+", model)[1:]
                               if x != "0"]

            r_model = ro.default_py2ri(model)

        else:
            DEtype = "pairwise"

        r_DEtype = ro.default_py2ri(DEtype)

        E.info('running DESeq2: groups=%s, pairs=%s, replicates=%s, pairs=%s, '
               'additional_factors:' %
               (design.groups, design.pairs, design.has_replicates,
                design.has_pairs))
        E.info(design.factors)

        # load DESeq
        R('''suppressMessages(library('DESeq2'))''')

        # build design matrix
        if DEtype == "pairwise":
            buildDesign = R('''
              function(counts, groups){

                design = data.frame(row.names = colnames(counts),
                                    condition = groups)
                return(design)}''')

            r_design = buildDesign(r_counts, r_groups)

            buildCountDataSet = R('''
            function(counts, design){

            for(column in colnames(design)){
              design[[column]] = factor(design[[column]])
            }

            dds <- suppressMessages(DESeqDataSetFromMatrix(
                     countData= counts,
                     colData = design,
                     design = ~condition))

            return(dds)
            }''' % locals())

            r_dds = buildCountDataSet(r_counts, r_design)

            performDifferentialTesting = R('''
            function(dds){
            dds = suppressMessages(DESeq(dds))

            png("%(outfile_prefix)s_dispersion.png")
            plotDispEsts(dds)
            dev.off()

            res = suppressMessages(results(dds, addMLE=TRUE))
            res = as.data.frame(res)

            contrast = "condition"
            res$contrast = contrast
            contrast_levels = levels(dds@colData[[contrast]])

            if(length(contrast_levels)==2){
              res$control = contrast_levels[1]
              res$treatment = contrast_levels[2]}
            else{
              res$control = contrast
              res$treatment = contrast}

            return(res)}''' % locals())

            results = pandas2ri.ri2py(performDifferentialTesting(r_dds))
            results['test_id'] = results.index

        # DEtype == "GLM"
        else:
            r_design = r_factors_df

            buildCountDataSet = R('''
            function(counts, design, model){

            for(column in colnames(design)){
              design[[column]] = factor(design[[column]])
            }

            full_model <- formula("%(model)s")

            dds <- suppressMessages(DESeqDataSetFromMatrix(
                     countData= counts,
                     colData = design,
                     design = full_model))

            return(dds)
            }''' % locals())

            r_dds = buildCountDataSet(r_counts, r_design, r_model)

            results = pandas.DataFrame()

            n = 0
            for contrast in contrasts:
                assert contrast in design.factors.columns, "contrast not found in\
                design factors columns"
                model = [x for x in model_terms if x != contrast]
                model = "~" + "+".join(model)

                performDifferentialTesting = R('''
                function(dds){

                ddsLRT = suppressMessages(
                  DESeq(dds, reduced=formula(%(model)s),betaPrior=TRUE))

                png("%(outfile_prefix)s_dispersion.png")
                plotDispEsts(ddsLRT)
                dev.off()

                contrast_levels = as.vector(levels(dds@colData$%(contrast)s))

                for(levels in combn(contrast_levels, 2, simplify=F)){


                    res = suppressMessages(results(ddsLRT, addMLE=TRUE,
                                  contrast=c("%(contrast)s",
                                  levels[1], levels[2])))

                    png(paste0(c("%(outfile_prefix)s", "%(contrast)s",
                               levels[1], levels[2], "MA.png"), collapse="_"))
                    plotMA(res, alpha=%(fdr)s)
                    dev.off()
                }

                res = as.data.frame(res)
                res$contrast = "%(contrast)s"

                if(length(contrast_levels)==2){
                  tmp_df = data.frame(contrast_levels)
                  res$control = tmp_df[1,1]
                  res$treatment = tmp_df[2,1]
                  }
                else{
                  res$control = "%(contrast)s"
                  res$treatment = "%(contrast)s"}

                return(res)}''' % locals())

                tmp_results = pandas2ri.ri2py(performDifferentialTesting(r_dds))
                tmp_results['test_id'] = tmp_results.index

                # need to set index to sequence of ints to avoid duplications
                n2 = n+tmp_results.shape[0]
                tmp_results.index = range(n, n2)
                n = n2

                results = results.append(tmp_results)

        final_result = DEResult_DESeq2(testTable=results)

        return final_result


class DEResult_DESeq2(DEResult):

    def getResults(self, fdr):
        ''' post-process test results table into generic results output '''

        E.info("Generating output - results table")

        df_dict = collections.defaultdict()

        n_rows = self.table.shape[0]

        df_dict["treatment_name"] = self.table['control']
        df_dict["control_name"] = self.table['treatment']
        df_dict["test_id"] = self.table['test_id']
        df_dict["contrast"] = self.table['contrast']
        df_dict["control_mean"] = self.table['baseMean']
        df_dict["treatment_mean"] = self.table['baseMean']
        df_dict["control_std"] = (0,)*n_rows
        df_dict["treatment_std"] = (0,)*n_rows
        df_dict["p_value"] = self.table['pvalue']
        df_dict["p_value_adj"] = adjustPvalues(self.table['pvalue'])
        df_dict["significant"] = pvaluesToSignficant(
            df_dict["p_value_adj"], fdr)
        df_dict["l2fold"] = self.table['lfcMLE']

        # Transformed l2fold is the shrunken values
        df_dict["transformed_l2fold"] = self.table['log2FoldChange']

        # TS: check what happens when no fold change is available
        # TS: may need an if/else in list comprehension. Raise E.warn too?
        df_dict["fold"] = [math.pow(2, float(x)) for
                           x in self.table['log2FoldChange']]

        # set all status values to "OK"
        # TS: again, may need an if/else to check...
        df_dict["status"] = ("OK",)*n_rows

        self.table = pandas.DataFrame(df_dict)
        # causes errors if multiple instance of same test_id exist, for example
        # if multiple constrasts have been tested
        # self.table.set_index("test_id", inplace=True)


###############################################################################


def buildProbeset2Gene(infile,
                       outfile,
                       database="hgu133plus2.db",
                       mapping="hgu133plus2ENSEMBL"):
    '''build map relating a probeset to an ENSEMBL gene_id'''

    R.library(database)

    # map is a Bimap object
    m = R(mapping)

    result = R.toTable(m)

    outf = open(outfile, "w")
    outf.write("probe_id\tgene_id\n")
    for probeset_id, gene_id in zip(result["probe_id"],
                                    result["ensembl_id"]):
        outf.write("%s\t%s\n" % (probeset_id, gene_id))
    outf.close()

    E.info("written %i mappings to %s: probes=%i, genes=%i" %
           (len(result),
            outfile,
            len(set(result["probe_id"])),
            len(set(result["ensembl_id"]))))

GeneExpressionResult = collections.namedtuple(
    "GeneExpressionResult",
    "test_id treatment_name treatment_mean treatment_std "
    "control_name control_mean control_std "
    "pvalue qvalue l2fold fold transformed_l2fold "
    "significant status")


def writeExpressionResults(outfile, result):
    '''output expression results table.'''
    if outfile == sys.stdout:
        outf = outfile
    else:
        outf = IOTools.openFile(outfile, "w")

    outf.write("%s\n" % "\t".join(GeneExpressionResult._fields))
    for x in sorted(result):
        outf.write("%s\n" % "\t".join(map(str, x)))

    if outf != sys.stdout:
        outf.close()


class WelchsTTest(object):

    '''base class for computing expression differences.
    '''

    def __call__(self,
                 probesets,
                 treatments,
                 controls):

        assert len(probesets) == len(treatments[0])
        assert len(probesets) == len(controls[0])

        nskipped = 0
        results = []

        for probeset, treatment, control in zip(
                probesets, zip(*treatments), zip(*controls)):

            nval1, nval2 = len(treatment), len(control)
            mean1, mean2 = numpy.mean(treatment), numpy.mean(control)
            stddev1, stddev2 = numpy.std(treatment), numpy.std(control)

            try:
                s = Stats.doWelchsTTest(nval1, mean1, stddev1,
                                        nval2, mean2, stddev2,
                                        alpha=0.05)
            except ValueError:
                E.warn(
                    "expressionDifferences: standard deviations are 0 for "
                    "probeset %s - skipped" % probeset)
                nskipped += 1
                continue

            s.mProbeset = probeset
            results.append(s)

        qvalues = Stats.doFDR([x.mPValue for x in results]).mQValues

        for s, qvalue in zip(results, qvalues):
            s.mQValue = qvalue

        return results, nskipped


class SAMR(object):

    '''SAM analysis of microarray data.

    Use the Two-Class Unpaired Case Assuming Unequal Variances.

    This uses the samr library.

    Significant genes are either called at *fdr* or the
    top *ngenes* are returned.

    *treatments* and *control* are arrays of
    arrays of expression values.

    See

    https://stat.ethz.ch/pipermail/bioconductor/2008-July/023251.html

    for an explanation of the differences between siggens SAM
    and Excel SAM. This version is parameterised to reproduce Excel SAM
    by setting::

       var.equal = TRUE
       med = TRUE

    .. note::
        SAM requires log2 scaled expression levels.
    '''

    def __call__(self, probesets,
                 treatments,
                 controls,
                 pattern=None,
                 fdr=0.10,
                 ngenes=None,
                 npermutations=1000,
                 ndelta=10,
                 method="ttest"):

        if ngenes and fdr:
            raise ValueError("either supply ngenes or fdr, but not both.")

        R.library("samr")

        m = numpy.matrix(treatments + controls)
        m = numpy.transpose(m)
        labels = numpy.array([1] * len(treatments) + [2] * len(controls))

        R.assign("x", numpy.array(m))
        R.assign("y", labels)
        R.assign("probesets", probesets)

        data = R(
            '''data=list( x=x, y=y, geneid=1:length(probesets), genenames=probesets, logged2=TRUE)''')
        result = R(
            '''samr.obj<-samr(data,  resp.type="Two class unpaired", nperms=100)''')
        R('''plot(samr.obj, delta=.4)''')


class SAM(object):

    '''SAM analysis of microarray data.

    Use the Two-Class Unpaired Case Assuming Unequal Variances.

    This uses the siggenes library. Note that there is also
    an rsam package at:

    http://rss.acs.unt.edu/Rdoc/library/samr/html/samr.html

    Significant genes are either called at *fdr* or the
    top *ngenes* are returned.

    *treatments* and *control* are arrays of
    arrays of expression values.

    See
    https://stat.ethz.ch/pipermail/bioconductor/2008-July/023251.html

    for an explanation of the differences between siggens SAM
    and Excel SAM. To parameterize the FDR to excel sam, set the
    flag *use_excel_sam*.

    .. note::
        SAM requires log2 scaled expression levels.

    I ran into trouble using this library. I was not able to
    reproduce the same results from the original SAM study getting
    differences in d and in the fdr.

    fold change is treatment / control.

    '''

    def __call__(self, probesets,
                 treatments,
                 controls,
                 pattern=None,
                 fdr=0.10,
                 ngenes=None,
                 npermutations=1000,
                 ndelta=10,
                 method="ttest",
                 use_excel_sam=False,
                 treatment_label="treatment",
                 control_label="control"):

        if ngenes and fdr:
            raise ValueError("either supply ngenes or fdr, but not both.")

        R.library("siggenes")

        m = numpy.matrix(treatments + controls)
        m = numpy.transpose(m)

        E.debug("build expression matrix: %i x %i" % m.shape)

        labels = numpy.array([1] * len(treatments) + [0] * len(controls))
        # 1000 permutations for P-Values of down to 0.0001. Setting this
        # to a high value improved reproducibility of results.

        kwargs = {}
        # kwargs set to replicate excel SAM
        if use_excel_sam:
            kwargs.update(
                {"control":
                 R('''samControl( lambda = 0.5, n.delta = %(ndelta)s) ''' %
                   locals()),
                 "med": True,
                 "var.equal": True})
        else:
            kwargs.update({"control":
                           R('''samControl( n.delta = %(ndelta)s ) ''' %
                             locals())},)

        # the option B needs to be not set if wilc.stat is chosen

        if method == "ttest":
            kwargs["method"] = R('''d.stat''')
            kwargs["B"] = npermutations
        elif method == "wilc":
            kwargs["method"] = R('''wilc.stat''')
        elif method == "cat":
            kwargs["method"] = R('''cat.stat''')
        else:
            raise ValueError("unknown statistic `%s`" % method)

        E.info("running sam with the following options: %s" % str(kwargs))

        a = R.sam(numpy.array(m),
                  labels,
                  gene_names=numpy.array(probesets),
                  **kwargs)

        # E.debug("%s" % str(a))

        R.assign("a", a)

        fdr_data = collections.namedtuple("sam_fdr", (
            "delta", "p0", "false", "significant", "fdr", "cutlow",
            "cutup", "j2", "j1"))
        cutoff_data = collections.namedtuple(
            "sam_cutoff", ("delta", "significant", "fdr"))
        gene_data = collections.namedtuple(
            "sam_fdr", ("row", "dvalue", "stddev", "rawp", "qvalue", "rfold"))

        def _totable(robj):
            '''convert robj to a row-wise table.'''
            s = numpy.matrix(robj)
            t = [numpy.array(x).reshape(-1,) for x in s]
            return t

        # extract the fdr values
        # returns R matrix
        t = _totable(a.do_slot('mat.fdr'))
        assert len(t[0]) == len(fdr_data._fields)
        # for x in t: E.debug( "x=%s" % str(x))
        fdr_values = [fdr_data(*x) for x in t]

        # find d cutoff
        if fdr is not None and fdr > 0:
            s = numpy.matrix(R.findDelta(a, fdr))
            try:
                cutoffs = [cutoff_data(*numpy.array(x).reshape(-1,))
                           for x in s]
                E.debug("sam cutoffs for fdr %f: %s" % (fdr, str(cutoffs)))
                cutoff = cutoffs[-1]
            except TypeError:
                E.debug("could not get cutoff")
                cutoff = None
        elif ngenes:
            s = numpy.matrix(R.findDelta(a, ngenes))
            try:
                cutoffs = [cutoff_data(*numpy.array(x).reshape(-1,))
                           for x in s]
                E.debug("sam cutoffs for fdr %f: %s" % (fdr, str(cutoffs)))
                cutoff = cutoffs[-1]
            except TypeError:
                E.debug("could not get cutoff")
                cutoff = None
        else:
            raise ValueError("either supply ngenes or fdr")

        # collect (unadjusted) p-values and qvalues for all probesets
        pvalues = dict(zip(probesets, R('''a@p.value''')))
        qvalues = dict(zip(probesets, R('''a@q.value''')))

        if pattern:
            outfile = pattern % "sam.pdf"
            R.pdf(outfile)
            if cutoff:
                R.plot(a, cutoff.delta)
            else:
                R.plot(a)
            R['dev.off']()

        siggenes = {}
        significant_genes = set()
        if cutoff is not None:
            E.debug("using cutoff %s" % str(cutoff))

            summary = R('''summary( a, %f )''' % cutoff.delta)

            # summary = R.summary( a, cutoff.delta )
            R.assign("summary", summary)

            significant_genes = set(
                [probesets[int(x) - 1] for x
                 in R('''summary@row.sig.genes''')])
            # E.debug( "significant genes=%s" % str(significant_genes))

            r_result = zip(*_totable(summary.do_slot('mat.sig')))

            if len(r_result) > 0:

                assert len(r_result[0]) == 6, \
                    "expected six columns from siggenes module, got: %s" % \
                    len(r_result[0])

                for x in r_result:
                    if x[4] > fdr:
                        E.warn("%s has qvalue (%f) larger than cutoff, but "
                               "is significant significant." % (str(x), x[4]))

                # except TypeError:
                # only a single value
                #     x = [r_result[y] for y in ("Row", "d.value", "stdev", "rawp", "q.value", "R.fold") ]
                #     if x[4] > fdr:
                #         E.warn( "%s has qvalue (%f) larger than cutoff, but is called significant." % (str(x), x[4]))

                siggenes[probesets[int(x[0]) - 1]] = gene_data(*x)

        else:
            E.debug("no cutoff found - no significant genes.")

        genes = []
        for probeset, treatment, control in zip(
                probesets, zip(*treatments), zip(*controls)):

            mean1, mean2 = numpy.mean(treatment), numpy.mean(control)

            if probeset in siggenes:
                s = siggenes[probeset]
                pvalue = s.rawp
                qvalue = s.qvalue
            else:
                pvalue = pvalues[probeset]
                qvalue = qvalues[probeset]

            significant = (0, 1)[probeset in significant_genes]

            genes.append(GeneExpressionResult._make((probeset,
                                                     treatment_label,
                                                     mean1,
                                                     numpy.std(treatment),
                                                     control_label,
                                                     mean2,
                                                     numpy.std(control),
                                                     pvalue,
                                                     qvalue,
                                                     mean1 - mean2,
                                                     math.pow(
                                                         2, mean1 - mean2),
                                                     math.pow(
                                                         2, mean1 - mean2),
                                                     significant,
                                                     "OK")))

        return genes, cutoff, fdr_values


#########################################################################
#########################################################################
#########################################################################
def loadTagData(tags_filename, design_filename):
    '''load tag data for deseq/edger analysis.

    *Infile* is a tab-separated file with counts.

    *design_file* is a tab-separated file with the
    experimental design with a minimum of four columns::

      track   include group   pair
      CW-CD14-R1      0       CD14    1
      CW-CD14-R2      0       CD14    1
      CW-CD14-R3      1       CD14    1
      CW-CD4-R1       1       CD4     1
      FM-CD14-R1      1       CD14    2
      FM-CD4-R2       0       CD4     2
      FM-CD4-R3       0       CD4     2
      FM-CD4-R4       0       CD4     2

    track
        name of track - should correspond to column header in *infile*
    include
        flag to indicate whether or not to include this data
    group
        group indicator - experimental group
    pair
        pair that sample belongs to (for paired tests)

    Additional columns in design file are taken to contain levels for
    additional factors and may be included for tests that allow multi-factor
    model designs.

    This method creates various R objects:

    countsTable : data frame with counts.
    groups : vector with groups
    pairs  : vector with pairs
    factors : df of additional factors for more complex model designs

    '''

    # Load counts table
    E.info("loading tag data from %s" % tags_filename)

    R('''counts_table = read.table('%(tags_filename)s',
    header=TRUE,
    row.names=1,
    stringsAsFactors=TRUE,
    comment.char='#')''' % locals())

    E.info("read data: %i observations for %i samples" %
           tuple(R('''dim(counts_table)''')))
    E.debug("sample names: %s" % R('''colnames(counts_table)'''))

    # Load comparisons from file
    R('''pheno = read.delim(
    '%(design_filename)s',
    header=TRUE,
    stringsAsFactors=TRUE,
    comment.char='#')''' % locals())

    # Make sample names R-like - substitute - for .
    R('''pheno[,1] = gsub('-', '.', pheno[,1]) ''')
    E.debug("design names: %s" % R('''pheno[,1]'''))

    # Ensure pheno rows match count columns
    pheno = R(
        '''pheno2 = pheno[match(colnames(counts_table),pheno[,1]),,drop=FALSE]''')
    missing = R('''colnames(counts_table)[is.na(pheno2)][1]''')
    if missing:
        E.warn("missing samples from design file are ignored: %s" %
               missing)

    # Subset data & set conditions
    R('''includedSamples <- !(is.na(pheno2$include) | pheno2$include == '0') ''')
    E.debug("included samples: %s" %
            R('''colnames(counts_table)[includedSamples]'''))
    R('''countsTable <- counts_table[ , includedSamples ]''')
    R('''groups <- factor(pheno2$group[ includedSamples ])''')
    R('''conds <- pheno2$group[ includedSamples ]''')
    R('''pairs <- factor(pheno2$pair[ includedSamples ])''')

    # JJ if additional columns present, pass to 'factors'
    R('''if (length(names(pheno2)) > 4) {
           factors <- data.frame(pheno2[includedSamples,5:length(names(pheno2))])
         } else {
           factors <- NA
         }''')

    E.info("filtered data: %i observations for %i samples" %
           tuple(R('''dim(countsTable)''')))


def filterTagData(filter_min_counts_per_row=1,
                  filter_min_counts_per_sample=10,
                  filter_percentile_rowsums=0):
    '''filter tag data.

    * remove rows with fewer than x counts in most highly expressed sample

    * remove samples with fewer than x counts in most highly expressed row

    * remove the lowest percentile of rows in the table, sorted
       by total tags per row
    '''

    # Remove windows with no data
    R('''max_counts = apply(countsTable,1,max)''')
    R('''countsTable = countsTable[max_counts>=%i,]''' %
      filter_min_counts_per_row)
    E.info("removed %i empty rows" %
           tuple(R('''sum(max_counts == 0)''')))
    observations, samples = tuple(R('''dim(countsTable)'''))
    E.info("trimmed data: %i observations for %i samples" %
           (observations, samples))

    # remove samples without data
    R('''max_counts = apply(countsTable,2,max)''')

    empty_samples = tuple(
        R('''max_counts < %i''' % filter_min_counts_per_sample))
    sample_names = R('''colnames(countsTable)''')
    nempty_samples = sum(empty_samples)

    if nempty_samples:
        E.warn("%i empty samples are being removed: %s" %
               (nempty_samples,
                ",".join([sample_names[x]
                          for x, y in enumerate(empty_samples) if y])))
        R('''countsTable <- countsTable[, max_counts >= %i]''' %
          filter_min_counts_per_sample)
        R('''groups <- groups[max_counts >= %i]''' %
          filter_min_counts_per_sample)
        R('''pairs <- pairs[max_counts >= %i]''' %
          filter_min_counts_per_sample)
        R('''if (!is.na(factors)) {factors <- factors[max_counts >= %i,]}''' %
          filter_min_counts_per_sample)
        observations, samples = tuple(R('''dim(countsTable)'''))

    # percentile filtering
    if filter_percentile_rowsums > 0:
        percentile = float(filter_percentile_rowsums) / 100.0
        R('''sum_counts = rowSums( countsTable )''')
        R('''take = (sum_counts > quantile(sum_counts, probs = %(percentile)f))''' %
          locals())
        discard, keep = R('''table( take )''')
        E.info("percentile filtering at level %f: keep=%i, discard=%i" %
               (filter_percentile_rowsums,
                keep, discard))
        R('''countsTable = countsTable[take,]''')

    observations, samples = tuple(R('''dim(countsTable)'''))

    return observations, samples


def groupTagData(ref_group=None):
    '''compute groups and pairs from tag data table.'''

    # Relevel the groups so that the reference comes first
    if ref_group is not None:
        R('''groups <- relevel(groups, ref = "%s")''' % ref_group)

    groups = R('''levels(groups)''')
    pairs = R('''levels(pairs)''')
    factors = R('''factors''')

    # JJ - check whether there are additional factors in design file...
    # warning... isintance(df, rpy.robjects.vectors.Vector) returns True
    if isinstance(factors, rpy2.robjects.vectors.DataFrame):
        E.warn("There are additional factors in design file that are ignored"
               " by groupTagData: %s" % factors.r_repr())
    else:
        # Hack... must be a better way to evaluate r NA instance in python?
        assert len(list(factors)) == 1 and bool(list(factors)[0]) is False, \
            "factors must either be DataFrame or NA in R global namespace"

    # Test if replicates exist - at least one group must have multiple samples
    max_per_group = R('''max(table(groups)) ''')[0]

    has_replicates = max_per_group >= 2

    # Test if pairs exist:
    npairs = R('''length(table(pairs)) ''')[0]
    has_pairs = npairs == 2

    # at least two samples per pair
    if has_pairs:
        min_per_pair = R('''min(table(pairs)) ''')[0]
        has_pairs = min_per_pair >= 2

    return groups, pairs, has_replicates, has_pairs


def plotCorrelationHeatmap(method="correlation"):
    '''plot a heatmap of correlations derived from
    countsTable.
    '''

    if method == "correlation":
        R('''dists <- dist( (1 - cor(countsTable)) / 2 )''')
    else:
        R('''dists <- dist( t(as.matrix(countsTable)), method = '%s' )''' %
          method)

    R('''heatmap( as.matrix( dists ), symm=TRUE )''')


def plotPairs():
    '''requires counts table'''
    # Plot pairs
    R('''panel.pearson <- function(x, y, digits=2, prefix="", cex.cor, ...)
            {
            usr <- par("usr"); on.exit(par(usr))
            par(usr = c(0, 1, 0, 1))
            r <- abs(cor(x, y))
            txt <- format(c(r, 0.123456789), digits=digits)[1]
            txt <- paste(prefix, txt, sep="")
            if(missing(cex.cor)) cex <- 0.6/strwidth(txt)
            x = 0.5;
            y = 0.5;
            if (par("xlog")) { x = 10^x };
            if (par("ylog")) { y = 10^y };
            text(x, y, txt, cex = cex);
            }
       ''')
    try:
        R('''pairs(countsTable,
        lower.panel = panel.pearson,
        pch=".",
        labels=colnames(countsTable),
        log="xy")''')
    except rpy2.rinterface.RRuntimeError, msg:
        E.warn("can not plot pairwise scatter plot: %s" % msg)


def plotPCA(groups=True):
    '''plot a PCA plot from countsTable using ggplot.

    If groups is *True*, the variable ``groups`` is
    used for colouring. If *False*, the groups are
    determined by sample labels.
    '''
    R('''suppressMessages(library(ggplot2))''')
    R('''pca = prcomp(t(countsTable))''')
    # Build factor groups by splitting labels at "."
    R('''colour=groups''')
    R('''shape=0''')
    R('''size=1''')
    if groups is False:
        R('''mm = matrix(
        unlist(sapply(colnames(countsTable),strsplit,'[.]')),
        nrow=length(colnames(countsTable)),
        byrow=T)''')
        nrows, nlevels = R('''dim(mm)''')
        if nlevels > 1:
            R('''colour=mm[,1]''')
        if nlevels > 2:
            R('''shape=mm[,2]''')

    try:
        R('''p1 = ggplot(
        as.data.frame(pca$x),
        aes(x=PC1, y=PC2,
        colour=colour,
        shape=shape,
        label=rownames(pca$x))) \
        + geom_text(size=4, vjust=1) \
        + geom_point()''')
        R('''p2 = qplot(x=PC1, y=PC3,
        data = as.data.frame(pca$x),
        label=rownames(pca$x),
        shape=shape,
        colour=colour)''')
        R('''p3 = qplot(x=PC2, y=PC3,
        data = as.data.frame(pca$x),
        label=rownames(pca$x),
        shape=shape,
        colour=colour)''')
        # TODO: plot all in a multi-plot with proper scale
        # the following squishes the plots
        # R('''source('%s')''' %
        #   os.path.join(os.path.dirname(E.__file__),
        #                "../R",
        #                "multiplot.R"))
        # R('''multiplot(p1, p2, p3, cols=2)''')
        R('''plot(p1)''')
    except rpy2.rinterface.RRuntimeError, msg:
        E.warn("could not plot in plotPCA(): %s" % msg)


def runEdgeR(outfile,
             outfile_prefix="edger.",
             fdr=0.1,
             prefix="",
             dispersion=None,
             ref_group=None
             ):
    '''run EdgeR on countsTable.

    Results are stored in *outfile* and files prefixed by *outfile_prefix*.

    The dispersion is usually measuered from replicates. If there are no
    replicates, you need to set the *dispersion* explicitely.

    See page 13 of the EdgeR user guide::

       2. Simply pick a reasonable dispersion value, based on your
       experience with similar data, and use that. Although
       subjective, this is still more defensible than assuming Poisson
       variation. Typical values are dispersion=0.4 for human data,
       dispersion=0.1 for data on genetically identical model
       organisms or dispersion=0.01 for technical replicates.

    '''

    # load library
    R('''suppressMessages(library('edgeR'))''')

    groups, pairs, has_replicates, has_pairs = groupTagData(ref_group)

    # output heatmap plot
    print "outfile_prefix:"
    print outfile_prefix
    print '%(outfile_prefix)sheatmap.png' % locals()
    R.png('%(outfile_prefix)sheatmap.png' % locals())
    plotCorrelationHeatmap()
    R['dev.off']()

    E.info('running EdgeR: groups=%s, pairs=%s, replicates=%s, pairs=%s' %
           (groups, pairs, has_replicates, has_pairs))

    if has_pairs:
        # output difference between groups
        R.png('''%(outfile_prefix)sbalance_groups.png''' % locals())
        first = True
        for g1, g2 in itertools.combinations(groups, 2):
            R('''a = rowSums( countsTable[groups == '%s'] ) ''' % g1)
            R('''b = rowSums( countsTable[groups == '%s'] ) ''' % g2)
            if first:
                R('''plot(cumsum(sort(a - b)), type = 'l') ''')
                first = False
            else:
                R('''lines(cumsum(sort(a - b))) ''')

        R['dev.off']()

        R('''suppressMessages(library('ggplot2'))''')
        R('''suppressMessages(library('reshape'))''')

        # output difference between pairs within groups
        first = True
        legend = []
        for pair in pairs:
            for g1, g2 in itertools.combinations(groups, 2):
                key = re.sub("-", "_", "pair_%s_%s_vs_%s" % (pair, g1, g2))
                legend.append(key)
                # print R('''colnames( countsTable) ''')
                # print R(''' pairs=='%s' ''' % pair)
                # print R(''' groups=='%s' ''' % g1)
                R('''a = rowSums( countsTable[pairs == '%s' & groups == '%s'])''' % (
                    pair, g1))
                R('''b = rowSums( countsTable[pairs == '%s' & groups == '%s'])''' % (
                    pair, g2))
                R('''c = cumsum( sort(a - b) )''')
                R('''c = c - min(c)''')
                if first:
                    data = R('''d = data.frame( %s = c)''' % key)
                    first = False
                else:
                    R('''d$%s = c''' % key)

        # remove row names (gene idenitifiers)
        R('''row.names(d) = NULL''')
        # add numbers of genes (x-axis)
        R('''d$genes=1:nrow(d)''')

        # merge data for ggplot
        R('''d = melt( d, 'genes', variable_name = 'comparison' )''')

        # plot
        R('''gp = ggplot(d)''')
        R('''pp = gp + \
            geom_line(aes(x=genes,y=value,group=comparison,color=comparison))''')

        try:
            R.ggsave('''%(outfile_prefix)sbalance_pairs.png''' % locals())
            R['dev.off']()
        except rpy2.rinterface.RRuntimeError, msg:
            E.warn("could not plot: %s" % msg)

    # build DGEList object
    # ignore message: "Calculating library sizes from column totals"
    R('''countsTable = suppressMessages(DGEList(countsTable, group=groups))''')

    # Relevel groups to make the results predictable - IMS
    if ref_group is not None:
        R('''countsTable$samples$group <- relevel(countsTable$samples$group,
        ref = "%s")''' % ref_group)
    else:
        # if no ref_group provided use first group in groups
        R('''countsTable$sample$group <- relevel(countsTable$samples$group,
        ref = "%s")''' % groups[0])

    # calculate normalisation factors
    E.info("calculating normalization factors")
    R('''countsTable = calcNormFactors( countsTable )''')
    E.info("output")

    # output MDS plot
    R.png('''%(outfile_prefix)smds.png''' % locals())
    try:
        R('''plotMDS( countsTable )''')
    except rpy2.rinterface.RRuntimeError, msg:
        E.warn("can not plot mds: %s" % msg)
    R['dev.off']()

    # build design matrix
    if has_pairs:
        R('''design = model.matrix(~pairs + countsTable$samples$group)''')
    else:
        R('''design = model.matrix(~countsTable$samples$group)''')

    # R('''rownames(design) = rownames( countsTable$samples )''')
    # R('''colnames(design)[length(colnames(design))] = "CD4" ''' )

    # fitting model to each tag
    if has_replicates:
        # estimate common dispersion
        R('''countsTable = estimateGLMCommonDisp(countsTable, design)''')
        # estimate tagwise dispersion
        R('''countsTable = estimateGLMTagwiseDisp(countsTable, design)''')
        # fitting model to each tag
        R('''fit = glmFit(countsTable, design)''')
    else:
        # fitting model to each tag
        if dispersion is None:
            raise ValueError("no replicates and no dispersion")
        E.warn("no replicates - using a fixed dispersion value")
        R('''fit = glmFit(countsTable, design, dispersion=%f)''' %
          dispersion)

    # perform LR test
    R('''lrt = glmLRT(fit)''')

    E.info("Generating output")

    # output cpm table
    R('''suppressMessages(library(reshape2))''')
    R('''countsTable.cpm <- cpm(countsTable, normalized.lib.sizes=TRUE)''')
    R('''melted <- melt(countsTable.cpm)''')
    R('''names(melted) <- c("test_id", "sample", "ncpm")''')
    # melt columns are factors - convert to string for sorting
    R('''melted$test_id = levels(melted$test_id)[as.numeric(melted$test_id)]''')
    R('''melted$sample = levels(melted$sample)[as.numeric(melted$sample)]''')
    # sort cpm table by test_id and sample
    R('''sorted = melted[with(melted, order(test_id, sample)),]''')
    R('''gz = gzfile("%(outfile_prefix)scpm.tsv.gz", "w" )''' % locals())
    R('''write.table(sorted, file=gz, sep = "\t",
                     row.names=FALSE, quote=FALSE)''')
    R('''close(gz)''')

    # compute adjusted P-Values
    R('''padj = p.adjust(lrt$table$PValue, 'BH')''')

    rtype = collections.namedtuple("rtype", "lfold logCPM LR pvalue")

    # output differences between pairs
    if len(groups) == 2:
        R.png('''%(outfile_prefix)smaplot.png''' % locals())
        R('''plotSmear(countsTable, pair=c('%s'))''' % "','".join(groups))
        R('''abline(h=c(-2, 2), col='dodgerblue') ''')
        R['dev.off']()

    # I am assuming that logFC is the base 2 logarithm foldchange.
    # Parse results and parse to file
    results = []
    counts = E.Counter()

    for interval, data, padj in zip(
            R('''rownames(lrt$table)'''),
            zip(*R('''lrt$table''')),
            R('''padj''')):
        d = rtype._make(data)

        counts.input += 1

        # set significant flag
        if padj <= fdr:
            signif = 1
            counts.significant += 1
            if d.lfold > 0:
                counts.significant_over += 1
            else:
                counts.significant_under += 1
        else:
            signif = 0
            counts.insignificant += 1

        if d.lfold > 0:
            counts.all_over += 1
        else:
            counts.all_under += 1

        # is.na failed in rpy2 2.4.2
        if d.pvalue != R('''NA'''):
            status = "OK"
        else:
            status = "FAIL"

        counts[status] += 1

        try:
            fold = math.pow(2.0, d.lfold)
        except OverflowError:
            E.warn("%s: fold change out of range: lfold=%f" %
                   (interval, d.lfold))
            # if out of range set to 0
            fold = 0

        # fold change is determined by the alphabetical order of the factors.
        # Is the following correct?
        results.append(GeneExpressionResult._make((
            interval,
            groups[1],
            d.logCPM,
            0,
            groups[0],
            d.logCPM,
            0,
            d.pvalue,
            padj,
            d.lfold,
            fold,
            d.lfold,  # no transform of lfold
            str(signif),
            status)))

    writeExpressionResults(outfile, results)

    outf = IOTools.openFile("%(outfile_prefix)ssummary.tsv" % locals(), "w")
    outf.write("category\tcounts\n%s\n" % counts.asTable())
    outf.close()

# needs to put into class
##


def deseqPlotSizeFactors(outfile):
    '''plot size factors - requires cds object.'''
    R.png(outfile)
    R('''par(mar=c(8,4,4,2))''')
    R('''barplot( sizeFactors( cds ), main="size factors", las=2)''')
    R['dev.off']()


def deseqOutputSizeFactors(outfile):
    '''output size factors - requires cds object.'''
    size_factors = R('''sizeFactors( cds )''')
    samples = R('''names(sizeFactors(cds))''')
    with IOTools.openFile(outfile, "w") as outf:
        outf.write("sample\tfactor\n")
        for name, x in zip(samples, size_factors):
            outf.write("%s\t%s\n" % (name, str(x)))


def deseqPlotCorrelationHeatmap(outfile, vsd):
    '''plot a heatmap

    Use variance stabilized data in object vsd.
    Should be 'blind', as then the transform is
    not informed by the experimental design.
    '''

    # rpy2.4.2 - passing of arrays seems to be broken - do it in R
    # dists = R['as.matrix'](R.dist(R.t(R.exprs(vsd))))
    dists = R('''as.matrix(dist(t(exprs(vsd))))''')
    R.png(outfile)
    R['heatmap.2'](
        dists,
        trace='none',
        margin=ro.IntVector((10, 10)))
    R['dev.off']()


def deseqPlotGeneHeatmap(outfile,
                         data,
                         Rowv=False,
                         Colv=False):
    '''plot a heatmap of all genes

    Use variance stabilized data in object vsd.
    Should be 'blind', as then the transform is
    not informed by the experimental design.
    '''
    if len(data) == 0:
        return

    # do not print if not enough values in one
    # direction (single row or column)
    if min(R.dim(data)) < 2:
        return

    R.png(outfile, width=500, height=2000)
    hmcol = R.colorRampPalette(R['brewer.pal'](9, "GnBu"))(100)

    R['heatmap.2'](
        data,
        col=hmcol,
        trace="none",
        dendrogram="none",
        Rowv=Rowv,
        Colv=Colv,
        labRow=False,
        margin=ro.IntVector((5, 5)),
        lhei=ro.IntVector((1, 10)),
        key=False)

    R['dev.off']()


def deseqPlotPCA(outfile, vsd, max_genes=500):
    '''plot a PCA

    Use variance stabilized data in object vsd.
    Should be 'blind', as then the transform is
    not informed by the experimental design.
    '''
    R.png(outfile)
    #  if there are more than 500 genes (after filtering)
    #  use the 500 most variable in the PCA
    #  else use the number of genes
    R('''ntop = ifelse(as.integer(dim(vsd))[1] >= %(max_genes)i,
    %(max_genes)i,
    as.integer(dim(vsd))[1])''' % locals())
    try:
        R('''plotPCA(vsd)''')
    except rpy2.rinterface.RRuntimeError, msg:
        E.warn("can not plot PCA: %s" % msg)
    R['dev.off']()


def deseqPlotPairs(outfile):
    '''requires counts table'''
    # Plot pairs
    R.png(outfile, width=960, height=960)
    plotPairs()
    R['dev.off']()


def deseqPlotPvaluesAgainstRowsums(outfile):
    '''plot pvalues against row sum rank.

    This plot is useful to see if quantile filtering could
    be applied.
    '''

    R('''counts_sum = rowSums( countsTable )''')
    R.png(outfile)
    R('''plot( rank( counts_sum)/length(counts_sum),
               -log10( res$pval),
               pch = 16,
               cex= 0.1)''')

    R('''abline( a=3, b=0, col='red')''')
    R['dev.off']()


def deseqParseResults(control_name, treatment_name, fdr, vsd=False):
    '''parse deseq output.

    retrieve deseq results from object 'res' in R namespace.

    The 'res' object is a dataframe with the following columns (from the
    DESeq manual):

    id: The ID of the observable, taken from the row names of the
          counts slots.

    baseMean: The base mean (i.e., mean of the counts divided by the size
          factors) for the counts for both conditions

    baseMeanA: The base mean (i.e., mean of the counts divided by the size
          factors) for the counts for condition A

    baseMeanB: The base mean for condition B

    foldChange: The ratio meanB/meanA

    log2FoldChange: The log2 of the fold change

    pval: The p value for rejecting the null hypothesis 'meanA==meanB'

    padj: The adjusted p values (adjusted with 'p.adjust( pval,
          method="BH")')

    vsd_log2FoldChange: The log2 fold change after variance stabilization.
          This data field is not part of DESeq proper, but has been added
          in this module in the runDESeq() method.

    Here, 'conditionA' is 'control' and 'conditionB' is 'treatment'
    such that a foldChange of 2 means that treatment is twice
    upregulated compared to control.

    Returns a list of results.

    If vsd is True, the log fold change will be computed from the variance
    stabilized data.

    '''

    results = []
    isna = R["is.na"]

    # Get column names from output and edit
    names = list(R['res'].names)
    m = dict([(x, x) for x in names])
    m.update(dict(
        pval="pvalue",
        baseMeanA="value1",
        baseMeanB="value2",
        id="interval_id",
        log2FoldChange="lfold"))

    rtype = collections.namedtuple("rtype", names)
    counts = E.Counter()

    for data in zip(*R['res']):
        counts.input += 1

        d = rtype._make(data)

        # set significant flag
        if d.padj <= fdr:
            signif = 1
            counts.significant += 1
            if d.log2FoldChange > 0:
                counts.significant_over += 1
            else:
                counts.significant_under += 1
        else:
            signif = 0
            counts.insignificant += 1

        if d.log2FoldChange > 0:
            counts.all_over += 1
        else:
            counts.all_under += 1

        # set lfold change to 0 if both are not expressed
        if d.baseMeanA == 0.0 and d.baseMeanB == 0.0:
            d = d._replace(foldChange=0, log2FoldChange=0)

        if d.pval != R('''NA'''):
            status = "OK"
        else:
            status = "FAIL"

        counts[status] += 1

        counts.output += 1

        # check if our assumptions about the direction of fold change
        # are correct
        assert (d.foldChange > 1) == (d.baseMeanB > d.baseMeanA)

        # note that fold change is computed as second group (B) divided by
        # first (A)
        results.append(GeneExpressionResult._make((
            d.id,
            treatment_name,
            d.baseMeanB,
            0,
            control_name,
            d.baseMeanA,
            0,
            d.pval,
            d.padj,
            d.log2FoldChange,
            d.foldChange,
            d.transformed_log2FoldChange,
            str(signif),
            status)))

    return results, counts


def runDESeq(outfile,
             outfile_prefix="deseq.",
             fdr=0.1,
             prefix="",
             fit_type="parametric",
             dispersion_method="pooled",
             sharing_mode="maximum",
             ref_group=None,
             ):
    '''run DESeq on countsTable.

    Results are stored in *outfile* and files prefixed by *outfile_prefix*.

    The current analysis follows the analysis as outlined in version
    1.14.0

    DESeq ignores any pair information in the design matrix.

    The output is treatment and control. Fold change values are
    computed as treatment divided by control.

    '''

    # load library
    R('''suppressMessages(library('DESeq'))''')
    R('''suppressMessages(library('gplots'))''')
    R('''suppressMessages(library('RColorBrewer'))''')

    groups, pairs, has_replicates, has_pairs = groupTagData(ref_group)

    # Run DESeq
    # Create Count data object
    E.info("running DESeq: replicates=%s" % (has_replicates))
    R('''cds <-newCountDataSet( countsTable, groups) ''')

    # Estimate size factors
    R('''cds <- estimateSizeFactors( cds )''')

    no_size_factors = R('''is.na(sum(sizeFactors(cds)))''')[0]
    if no_size_factors:
        E.warn("no size factors - can not estimate - no output")
        return

    # estimate variance
    if has_replicates:
        E.info("replicates - estimating variance from replicates")
    else:
        E.info("no replicates - estimating variance with method='blind'")
        dispersion_method = "blind"

    E.info("dispersion_method=%s, fit_type=%s, sharing_mode=%s" %
           (dispersion_method, fit_type, sharing_mode))
    R('''cds <- estimateDispersions( cds,
    method='%(dispersion_method)s',
    fitType='%(fit_type)s',
    sharingMode='%(sharing_mode)s')''' % locals())

    # bring into python namespace
    cds = R('''cds''')

    # plot fit - if method == "pooled":
    if dispersion_method == "pooled":
        R.png('%sdispersion_estimates_pooled.png' %
              outfile_prefix)
        R.plotDispEsts(cds)
        R['dev.off']()
    elif not has_replicates:
        # without replicates the following error appears
        # in the rpy2.py2ri conversion:
        #   'dims' cannot be of length 0
        pass
    else:
        dispersions = R('''ls(cds@fitInfo)''')
        for dispersion in dispersions:
            R.png('%sdispersion_estimates_%s.png' %
                  (outfile_prefix, dispersion))
        R.plotDispEsts(cds, name=dispersion)
        R['dev.off']()

    # plot size factors
    deseqPlotSizeFactors('%(outfile_prefix)ssize_factors.png' % locals())

    # output size factors
    deseqOutputSizeFactors("%(outfile_prefix)ssize_factors.tsv" % locals())

    # plot scatter plots of pairs
    deseqPlotPairs('%(outfile_prefix)spairs.png' % locals())

    if dispersion_method not in ("blind",):
        # also do a blind dispersion estimate for
        # a variance stabilizing transform
        R('''cds_blind <- estimateDispersions( cds,
        method='blind',
        fitType='%(fit_type)s',
        sharingMode='%(sharing_mode)s')''' % locals())
    else:
        R('''cds_blind = cds''')

    # perform variance stabilization for log2 fold changes
    vsd = R('''vsd = varianceStabilizingTransformation(cds_blind)''')
    # output normalized counts (in order)
    # gzfile does not work with rpy 2.4.2 in python namespace
    # using R.gzfile, so do it in R-space

    R('''t = counts(cds, normalized=TRUE);
    write.table(t[order(rownames(t)),],
    file=gzfile('%(outfile_prefix)scounts.tsv.gz'),
    row.names=TRUE,
    col.names=NA,
    quote=FALSE,
    sep='\t') ''' % locals())

    # output variance stabilized counts (in order)
    R('''t = exprs(vsd);
    write.table(t[order(rownames(t)),],
    file=gzfile('%(outfile_prefix)svsd.tsv.gz'),
    row.names=TRUE,
    col.names=NA,
    quote=FALSE,
    sep='\t')
    ''' % locals())

    # plot correlation heatmap of variance stabilized data
    deseqPlotCorrelationHeatmap(
        '%scorrelation_heatmap.png' % outfile_prefix,
        vsd)

    # plot PCA
    deseqPlotPCA('%spca.png' % outfile_prefix,
                 vsd)

    # plot gene heatmap for all genes - order by average expression
    # subtract one to get numpy indices
    select = R.order(R.rowMeans(R.counts(cds)), decreasing=True)
    # the following uses R-based indexing
    deseqPlotGeneHeatmap(
        '%sgene_heatmap.png' % outfile_prefix,
        R['as.matrix'](R.exprs(vsd).rx(select)))

    # plot heatmap of top 200 expressed genes
    deseqPlotGeneHeatmap(
        '%sgene_heatmap_top200.png' % outfile_prefix,
        R['as.matrix'](R.exprs(vsd).rx(select[:200])))

    # Call diffential expression for all pairings of groups included in the
    # design

    all_results = []
    for combination in itertools.combinations(groups, 2):
        control, treatment = combination
        gfix = "%s_vs_%s_" % (control, treatment)

        outfile_groups_prefix = outfile_prefix + gfix
        E.info(("calling differential expression for "
                "control=%s vs treatment=%s") %
               (control, treatment))
        res = R('''res = nbinomTest(cds, '%s', '%s')''' % (control, treatment))

        # plot significance
        R.png('''%(outfile_groups_prefix)ssignificance.png''' % locals())
        R('''plot(
        res$baseMean,
        res$log2FoldChange,
        log="x",
        pch=20, cex=.1,
        col = ifelse( res$padj < %(fdr)s, "red", "black"))''' % locals())
        R['dev.off']()

        # plot pvalues against rowsums
        deseqPlotPvaluesAgainstRowsums(
            '%(outfile_groups_prefix)spvalue_rowsums.png' % locals())

        E.info("Generating output (%s vs %s)" % (control, treatment))

        # get variance stabilized fold changes - note the reversal of
        # treatment/control
        R('''vsd_l2f =
        (rowMeans(exprs(vsd)[,conditions(cds) == '%s', drop=FALSE])
        - rowMeans( exprs(vsd)[,conditions(cds) == '%s', drop=FALSE]))''' %
          (treatment, control))

        # plot vsd correlation, see Figure 14 in the DESeq manual
        # if you also want to colour by expression level
        R.png('''%(outfile_groups_prefix)sfold_transformation.png''' %
              locals())
        R('''plot(
        res$log2FoldChange, vsd_l2f,
        pch=20, cex=.1,
        col = ifelse( res$padj < %(fdr)s, "red", "black" ) )''' % locals())
        R['dev.off']()

        # plot heatmap of differentially expressed genes
        # plot gene heatmap for all genes - order by average expression
        padj_column = list(res.colnames).index('padj')
        select = R('''select = res['padj'] < %f''' % fdr)

        if R('''sum(select)''')[0] > 0:
            E.info('%s vs %s: plotting %i genes in heatmap' %
                   (treatment, control, len(select)))
            data = R.exprs(vsd).rx(select)

            if not isinstance(data, rpy2.robjects.vectors.FloatVector):
                order = R.order(R.rowMeans(data), decreasing=True)
                deseqPlotGeneHeatmap(
                    '%sgene_heatmap.png' % outfile_groups_prefix,
                    R['as.matrix'](data[order]),
                    Colv=False,
                    Rowv=True)
            else:
                E.warn('can not plot differentially expressed genes')
        else:
            E.warn('no differentially expressed genes at fdr %f' % fdr)

        # Plot pvalue histogram
        R.png('''%(outfile_groups_prefix)spvalue_histogram.png''' % locals())
        R('''pvalues = res$pval''')
        R('''hist(pvalues, breaks=50, col='skyblue' )''')
        R['dev.off']()

        # Plot diagnostic plots for FDR
        if has_replicates:
            R('''orderInPlot = order(pvalues)''')
            R('''showInPlot = (pvalues[orderInPlot] < 0.08)''')
            # Jethro - previously plotting x =
            # pvalues[orderInPlot][showInPlot]
            # pvalues[orderInPlot][showInPlot] contains all NA values
            # from pvalues which(showInPlot) doesn't... removing NA
            # values
            R('''true.pvalues <- pvalues[orderInPlot][showInPlot]''')
            R('''true.pvalues <- true.pvalues[is.finite(true.pvalues)]''')
            if R('''sum(showInPlot)''')[0] > 0:
                R.png('''%(outfile_groups_prefix)sfdr.png''' % locals())
                # failure when no replicates:
                # rpy2.rinterface.RRuntimeError:
                # Error in plot.window(...) : need finite 'xlim' values
                R('''plot( seq( along=which(showInPlot)),
                           true.pvalues,
                           pch='.',
                           xlab=expression(rank(p[i])),
                           ylab=expression(p[i]))''')
                R('''abline(a = 0, b = %(fdr)f / length(pvalues), col = "red")
                ''' % locals())
                R['dev.off']()
            else:
                E.warn('no p-values < 0.08')

        # Add log2 fold with variance stabilized l2fold value
        R('''res$transformed_log2FoldChange = vsd_l2f''')

        # Parse results and parse to file
        results, counts = deseqParseResults(control,
                                            treatment,
                                            fdr=fdr)

        all_results += results

        E.info(counts)
        outf = IOTools.openFile(
            "%(outfile_groups_prefix)ssummary.tsv" % locals(), "w")
        outf.write("category\tcounts\n%s\n" % counts.asTable())
        outf.close()

    writeExpressionResults(outfile, all_results)


def runDESeq2(outfile,
              outfile_prefix="deseq2.",
              fdr=0.1,
              ref_group=None,
              model=None,
              contrasts=None
              ):
    """
    Run DESeq2 on counts table.

    If no model is passed, then defaults to the group column in design file

    Does not make use of group tag data bc function doesn't accomodate
    multi-factor designs

    To Do: Parse results into standard output format.
    Fix fact that plotMA is hardcoded.
    """

    # load libraries
    R('''suppressMessages(library('DESeq2'))''')

    # Create metadata... this will eventually be a pandas dataframe
    if isinstance(R('''factors'''), rpy2.robjects.vectors.DataFrame):
        E.info("DESeq2: Merging additional factors in design file to"
               "create metadata table")
        R('''mdata <- cbind(groups, factors)''')
        mdata = tuple(R('''names(mdata)'''))
    else:
        R('''mdata <- data.frame(group=groups)''')
        mdata = "group"
    E.info("DESeq2 colData headers are: %s" % mdata)

    R('''design="~ group"''')

    # Check for model and that model terms are in metadata table
    if model:
        assert contrasts, "Must specifiy contrasts is model design provided"
        terms = set([x for x in re.split("\W", model) if x != ''])
        assert terms.issubset(mdata), \
            "DESeq2: design formula has terms not present in colData"
    else:
        if mdata != "group":
            E.warn("DESeq2 model specified, with no metadata in design file")
        terms = ["group", ]
        model = "~ group"
    E.info("DESeq2 design formula is: %s" % model)

    # Create DESeqDataSet, using countsTable, mdata, model
    R('''dds <- DESeqDataSetFromMatrix(countData=countsTable,
                                       colData=mdata,
                                       design=%(model)s)''' % locals())
    # WARNING: This is not done automatically... I don't know why?
    R('''colnames(dds) <- colnames(countsTable)''')
    E.info("Combined colData, design formula and counts table to create"
           " DESeqDataSet instance")

    # Rlog transform
    R('''rld <- rlog(dds)''')

    # Plot PCA of rlog transformed count data for top 500
    for factor in terms:
        outf = outfile_prefix + factor + "_PCAplot500.tiff"
        E.info("Creating PCA plot for factor: %s" % outf)
        R('''x <- plotPCA(rld, intgroup="%(factor)s")''' % locals())
        # R('''saveRDS(x, '%(outf)s')''' % locals())
        R('''tiff("%(outf)s")''' % locals())
        R('''plot(x)''')
        R('''dev.off()''')

    # Extract rlog transformed count data...
    rlog_out = IOTools.snip(outfile, ".tsv.gz") + "_rlog.tsv.gz"
    # R('''saveRDS(rld, "%s")''' % rlog_out)
    R('''write.table(assay(rld), file='%s', sep="\t", quote=FALSE)''' % rlog_out)

    # Run DESeq2
    R('''dds <- DESeq(dds)''')
    E.info("Completed DESeq2 differential expression analysis")

    # Extract contrasts...
    if contrasts:
        contrasts = (x.split(":") for x in contrasts.split(","))
    else:
        # created by loadTagData...
        groups = R('''levels(groups)''')
        contrasts = (("group",) + x for x in itertools.combinations(groups, 2))

    df_final = pandas.DataFrame()
    for combination in contrasts:
        variable, control, treatment = combination

        # Fetch results
        gfix = "%s_%s_vs_%s" % (variable, control, treatment)
        outfile_groups_prefix = outfile_prefix + gfix + "_MAplot.png"
        R('''res <- results(dds, contrast=c("%(variable)s",
                                            "%(treatment)s",
                                            "%(control)s"))''' % locals())
        E.info("Extracting contrast for levels %s (treatment) vs %s (control)"
               " for factor %s" % (treatment, control, variable))

        # plot MA plot
        R('''png("%(outfile_groups_prefix)s")''' % locals())
        R('''plotMA(res, alpha=%f)''' % fdr)
        R('''dev.off()''')
        E.info("Plotted MA plot for levels %s (treatment) vs %s (control)"
               " for factor %s" % (treatment, control, variable))

        # write data to outfile
        res_df = R('''res_df <- as.data.frame(res)''')
        df_out = com.load_data("res_df")
        df_out["treatment"] = [treatment, ]*len(df_out.index)
        df_out["control"] = [control, ]*len(df_out.index)
        df_out["variable"] = [variable, ]*len(df_out.index)
        df_out.to_csv(IOTools.openFile(outfile_prefix + gfix + ".tsv.gz", "w"),
                      sep="\t",
                      index_label="gene_id")
        E.info("Extracted results table for contrast  '%s' (treatment) vs '%s'"
               " (control) for factor '%s'" % (treatment, control, variable))

        # append to final dataframe
        df_out.reset_index(inplace=True)
        df_out.rename(columns={"index": "gene_id"}, inplace=True)
        df_final = df_final.append(df_out, ignore_index=True)

    # write final dataframe
    df_final.to_csv(IOTools.openFile(outfile, "w"), sep="\t", index=False)

Design = collections.namedtuple("Design", ("include", "group", "pair"))


def readDesignFile(design_file):
    '''reads a design file.'''

    design = collections.OrderedDict()
    with IOTools.openFile(design_file) as inf:
        for line in inf:
            if line.startswith("track"):
                continue
            track, include, group, pair = line.split("\t")[:4]
            if track in design:
                raise ValueError("duplicate track '%s'" % track)
            design[track] = Design._make((int(include), group, pair))
    return design


def plotTagStats(infile, design_file, outfile_prefix):
    '''provide summary plots for tag data.'''

    loadTagData(infile, design_file)

    nobservations, nsamples = filterTagData()

    if nobservations == 0:
        E.warn("no observations - no output")
        return

    if nsamples == 0:
        E.warn("no samples remain after filtering - no output")
        return

    groups, pairs, has_replicates, has_pairs = groupTagData()

    # import rpy2.robjects.lib.ggplot2 as ggplot2

    R('''suppressMessages(library('ggplot2'))''')
    R('''suppressMessages(library('reshape'))''')

    R('''d = melt( log10(countsTable + 1), variable_name = 'sample' )''')

    # Note that ggsave does not work if there is
    # X display.
    R.png(outfile_prefix + ".densities.png")
    R('''gp = ggplot(d)''')
    R('''pp = gp + geom_density(aes(x=value, group=sample,
    color=sample, fill=sample), alpha=I(1/3))''')
    R('''plot(pp)''')
    R['dev.off']()

    R.png(outfile_prefix + ".boxplots.png")
    R('''gp = ggplot(d)''')
    R('''pp = gp +
    geom_boxplot(aes(x=sample,y=value,color=sample,fill=sample),
    size=0.3,
    alpha=I(1/3)) +
    theme(axis.text.x = element_text( angle=90, hjust=1, size=8 ) )''')
    R('''plot(pp)''')
    R['dev.off']()


def plotDETagStats(infile, outfile_prefix,
                   additional_file=None,
                   join_columns=None,
                   additional_columns=None):
    '''provide summary plots for tag data.

    Stratify boxplots and densities according to differential
    expression calls.

    The input file is the output of any of the DE
    tools, see GeneExpressionResults for column names.

    Additional file will be joined with infile and any additional
    columns will be output as well.
    '''

    table = pandas.read_csv(IOTools.openFile(infile),
                            sep="\t")

    if additional_file is not None:
        additional_table = pandas.read_csv(
            IOTools.openFile(additional_file),
            sep="\t")
        table = pandas.merge(table,
                             additional_table,
                             on=join_columns,
                             how="left",
                             sort=False)

    # remove index. If it is numbered starting from 1, there is a bug
    # in ggplot, see https://github.com/yhat/ggplot/pull/384
    table.reset_index(inplace=True)

    # add log-transformed count data
    table['log10_treatment_mean'] = numpy.log10(table['treatment_mean'] + 1)
    table['log10_control_mean'] = numpy.log10(table['control_mean'] + 1)
    table['dmr'] = numpy.array(["insignicant"] * len(table))
    table.loc[
        (table["l2fold"] > 0) & (table["significant"] == 1), "dmr"] = "up"
    table.loc[
        (table["l2fold"] < 0) & (table["significant"] == 1), "dmr"] = "down"

    def _dplot(table, outfile, column):

        plot = ggplot.ggplot(
            ggplot.aes(column,
                       colour='dmr',
                       fill='dmr'),
            data=table) + \
            ggplot.geom_density(alpha=0.5)

        try:
            ggplot.ggsave(filename=outfile, plot=plot)
        except Exception, msg:
            E.warn("no plot for %s: %s" % (column, msg))

    def _bplot(table, outfile, column):

        plot = ggplot.ggplot(
            ggplot.aes(x=column, y='dmr'),
            data=table) + \
            ggplot.geom_boxplot()

        try:
            ggplot.ggsave(filename=outfile, plot=plot)
        except ValueError, msg:
            # boxplot fails if all values are the same
            # see https://github.com/yhat/ggplot/issues/393
            E.warn(msg)

    _dplot(table,
           outfile_prefix + ".densities_tags_control.png",
           "log10_control_mean")
    _dplot(table,
           outfile_prefix + ".densities_tags_treatment.png",
           "log10_treatment_mean")
    _bplot(table,
           outfile_prefix + ".boxplot_tags_control.png",
           "log10_control_mean")
    _bplot(table,
           outfile_prefix + ".boxplot_tags_treatment.png",
           "log10_treatment_mean")

    if additional_columns:
        for column in additional_columns:
            _dplot(table,
                   outfile_prefix + ".densities_%s.png" % column,
                   column)
            _bplot(table,
                   outfile_prefix + ".boxplot_%s.png" % column,
                   column)
    return


def runMockAnalysis(outfile,
                    outfile_prefix,
                    ref_group=None,
                    pseudo_counts=0):
    '''run a mock analysis on a count table.

    compute fold enrichment values, but do not normalize or
    perform any test.
    '''

    groups, pairs, has_replicates, has_pairs = groupTagData(ref_group)

    all_results = []
    for combination in itertools.combinations(groups, 2):
        control, treatment = combination

        R('''control_counts = rowSums( countsTable[groups == '%s'] )''' %
          control)
        R('''treatment_counts = rowSums( countsTable[groups == '%s'] )''' %
          treatment)

        # add pseudocounts to enable analysis of regions
        # that are absent/present
        if pseudo_counts:
            R('''control_counts = control_counts + %f''' % pseudo_counts)
            R('''treatment_counts = treatment_counts + %f''' % pseudo_counts)

        R('''fc = treatment_counts / control_counts''')

        results = []

        for identifier, treatment_count, control_count, foldchange in \
                zip(R('''rownames( countsTable)'''),
                    R('''treatment_counts'''),
                    R('''control_counts'''),
                    R('''fc''')):
            try:
                log2fold = math.log(foldchange)
            except ValueError:
                log2fold = "Inf"

            results.append(GeneExpressionResult._make((
                identifier,
                treatment,
                treatment_count,
                0,
                control,
                control_count,
                0,
                1,
                1,
                log2fold,
                foldchange,
                log2fold,
                "0",
                "OK")))

        all_results.extend(results)

    writeExpressionResults(outfile, all_results)


def outputTagSummary(filename_tags,
                     outfile,
                     output_filename_pattern,
                     filename_design=None):
    '''output summary values for a count table.'''

    E.info("loading tag data from %s" % filename_tags)

    if filename_design is not None:
        # load all tag data
        loadTagData(filename_tags, filename_design)

        # filter
        nobservations, nsamples = filterTagData()

    else:
        # read complete table
        R('''countsTable = read.delim('%(filename_tags)s',
        header = TRUE,
        row.names = 1,
        stringsAsFactors = TRUE,
        comment.char = '#')''' % locals())

        nobservations, nsamples = tuple(R('''dim(countsTable)'''))
        E.info("read data: %i observations for %i samples" %
               (nobservations, nsamples))
        # remove samples without data
        R('''max_counts = apply(countsTable,2,max)''')

        filter_min_counts_per_sample = 1

        empty_samples = tuple(
            R('''max_counts < %i''' % filter_min_counts_per_sample))
        sample_names = R('''colnames(countsTable)''')
        nempty_samples = sum(empty_samples)

        if nempty_samples:
            E.warn("%i empty samples are being removed: %s" %
                   (nempty_samples,
                    ",".join([sample_names[x]
                              for x, y in enumerate(empty_samples) if y])))
            R('''countsTable <- countsTable[, max_counts >= %i]''' %
              filter_min_counts_per_sample)
            nobservations, nsamples = tuple(R('''dim(countsTable)'''))

        R('''groups = factor(colnames( countsTable ))''')
        E.debug("sample names: %s" % R('''colnames(countsTable)'''))

    nrows, ncolumns = tuple(R('''dim(countsTable)'''))

    outfile.write("metric\tvalue\tpercent\n")
    outfile.write("number of observations\t%i\t100\n" % nobservations)
    outfile.write("number of samples\t%i\t100\n" % nsamples)

    # Count windows with no data
    R('''max_counts = apply(countsTable,1,max)''')

    # output distribution of maximum number of counts per window
    outfilename = output_filename_pattern + "max_counts.tsv.gz"
    E.info("outputting maximum counts per window to %s" % outfilename)
    R('''write.table(table(max_counts),
    file='%(outfilename)s',
    sep="\t",
    row.names=FALSE,
    quote=FALSE)''' %
      locals())

    # removing empty rows
    E.info("removing rows with no counts in any sample")
    R('''countsTable = countsTable[max_counts>0,]''')

    for x in range(0, 20):
        nempty = tuple(R('''sum(max_counts <= %i)''' % x))[0]
        outfile.write("max per row<=%i\t%i\t%f\n" %
                      (x, nempty, 100.0 * nempty / nrows))

    E.info("removed %i empty rows" % tuple(R('''sum(max_counts == 0)''')))
    observations, samples = tuple(R('''dim(countsTable)'''))
    E.info("trimmed data: %i observations for %i samples" %
           (observations, samples))

    # build correlation
    R('''correlations = cor(countsTable)''')
    outfilename = output_filename_pattern + "correlation.tsv"
    E.info("outputting sample correlations table to %s" % outfilename)
    R('''write.table(correlations, file='%(outfilename)s',
    sep="\t",
    row.names=TRUE,
    col.names=NA,
    quote=FALSE)''' % locals())

    # output scatter plots
    outfilename = output_filename_pattern + "scatter.png"
    E.info("outputting scatter plots to %s" % outfilename)
    R.png(outfilename, width=960, height=960)
    plotPairs()
    R['dev.off']()

    # output heatmap based on correlations
    outfilename = output_filename_pattern + "heatmap.svg"
    E.info("outputting correlation heatmap to %s" % outfilename)
    R.svg(outfilename)
    plotCorrelationHeatmap(method="correlation")
    R['dev.off']()

    # output PCA
    outfilename = output_filename_pattern + "pca.svg"
    E.info("outputting PCA plot to %s" % outfilename)
    R.svg(outfilename)
    plotPCA(groups=False)
    R['dev.off']()

    # output an MDS plot
    R('''suppressMessages(library('limma'))''')
    outfilename = output_filename_pattern + "mds.svg"
    E.info("outputting mds plot to %s" % outfilename)
    R.svg(outfilename)
    try:
        R('''plotMDS(countsTable)''')
    except rpy2.rinterface.RRuntimeError, msg:
        E.warn("can not plot mds: %s" % msg)
    R['dev.off']()


def dumpTagData(filename_tags, filename_design, outfile):
    '''output filtered tag table.'''

    if outfile == sys.stdout:
        outfilename = ""
    else:
        outfilename = outfile.name

    # load all tag data
    loadTagData(filename_tags, filename_design)

    # filter
    nobservations, nsamples = filterTagData()

    # output
    R('''write.table( countsTable,
                      file='%(outfilename)s',
                      sep='\t',
                      quote=FALSE)''' % locals())


def runTTest(outfile,
             outfile_prefix,
             fdr,
             ref_group=None):
    '''apply a ttest on the data.

    For the T-test it is best to use FPKM values as
    this method does not perform any library normalization.
    '''
    groups, pairs, has_replicates, has_pairs = groupTagData(ref_group)

    results = []
    for combination in itertools.combinations(groups, 2):
        control, treatment = combination
        r = R('''r = apply(countsTable, 1,
        function(x) { t.test(
        x[groups == '%(treatment)s'],
        x[groups == '%(control)s']) } )
        ''' % locals())

        for test_id, ttest in zip(r.names, r):
            # TS, swapped order below as assignment was incorrect
            treatment_mean, control_mean = tuple(ttest.rx2('estimate'))
            fold_change = treatment_mean / control_mean
            pvalue = tuple(ttest.rx2('p.value'))[0]
            significant = (0, 1)[pvalue < fdr]
            results.append(GeneExpressionResult._make((test_id,
                                                       treatment,
                                                       treatment_mean,
                                                       0,
                                                       control,
                                                       control_mean,
                                                       0,
                                                       pvalue,
                                                       pvalue,
                                                       numpy.log2(fold_change),
                                                       fold_change,
                                                       numpy.log2(fold_change),
                                                       significant,
                                                       "OK")))

    writeExpressionResults(outfile, results)


#####################################################################
# Pandas-based functions and matplotlib-based plotting functions ####
#####################################################################

def loadTagDataPandas(tags_filename, design_filename):
    '''load tag data for deseq/edger analysis.

    *Infile* is a tab-separated file with counts.

    *design_file* is a tab-separated file with the
    experimental design with four columns::

      track   include group   pair
      CW-CD14-R1      0       CD14    1
      CW-CD14-R2      0       CD14    1
      CW-CD14-R3      1       CD14    1
      CW-CD4-R1       1       CD4     1
      FM-CD14-R1      1       CD14    2
      FM-CD4-R2       0       CD4     2
      FM-CD4-R3       0       CD4     2
      FM-CD4-R4       0       CD4     2

    track
        name of track - should correspond to column header in *infile*
    include
        flag to indicate whether or not to include this data
    group
        group indicator - experimental group
    pair
        pair that sample belongs to (for paired tests)

    This method creates various R objects:

    countsTable : data frame with counts.
    groups : vector with groups
    pairs  : vector with pairs

    '''

    E.info("loading tag data from %s" % tags_filename)

    inf = IOTools.openFile(tags_filename)
    counts_table = pandas.read_csv(inf,
                                   sep="\t",
                                   index_col=0,
                                   comment="#")
    inf.close()

    E.info("read data: %i observations for %i samples" %
           counts_table.shape)

    E.debug("sample names: %s" % list(counts_table.columns))

    inf = IOTools.openFile(design_filename)
    design_table = pandas.read_csv(inf, sep="\t", index_col=0)
    inf.close()

    E.debug("design names: %s" % list(design_table.index))

    missing = set(counts_table.columns).difference(design_table.index)

    if missing:
        E.warn("missing samples from design file are ignored: %s" % missing)

    # remove unnecessary samples
    design_table = design_table[design_table["include"] != 0]
    E.debug("included samples: %s" % list(design_table.index))

    counts_table = counts_table[list(design_table.index)]
    E.info("filtered data: %i observations for %i samples" %
           counts_table.shape)

    return counts_table, design_table


def filterTagDataPandas(counts_table,
                        design_table,
                        filter_min_counts_per_row=1,
                        filter_min_counts_per_sample=10,
                        filter_percentile_rowsums=0):
    '''filter tag data.

    * remove rows with at least x number of counts

    * remove samples with a maximum of *min_sample_counts*

    * remove the lowest percentile of rows in the table, sorted
       by total tags per row
    '''

    # Remove windows with no data
    max_counts_per_row = counts_table.max(1)
    counts_table = counts_table[
        max_counts_per_row >= filter_min_counts_per_row]
    observations, samples = counts_table.shape
    E.info("trimmed data: %i observations for %i samples" %
           (observations, samples))

    # remove samples without data
    max_counts_per_sample = counts_table.max()

    empty_samples = max_counts_per_sample < filter_min_counts_per_sample
    sample_names = counts_table.columns
    nempty_samples = sum(empty_samples)

    if nempty_samples:
        E.warn("%i empty samples are being removed: %s" %
               (nempty_samples,
                ",".join([sample_names[x] for x, y in

                          enumerate(empty_samples) if y])))
        raise NotImplementedError("removing empty samples needs to be done")
        # R('''countsTable <- countsTable[, max_counts >= %i]''' % filter_min_counts_per_sample)
        # R('''groups <- groups[max_counts >= %i]''' % filter_min_counts_per_sample)
        # R('''pairs <- pairs[max_counts >= %i]''' % filter_min_counts_per_sample)
        # observations, samples = tuple( R('''dim(countsTable)'''))

    # percentile filtering
    if filter_percentile_rowsums > 0:
        percentile = float(filter_percentile_rowsums) / 100.0
        sum_counts = counts_table.sum(1)
        take = sum_counts > sum_counts.quantile(percentile)
        E.info("percentile filtering at level %f: keep=%i, discard=%i" %
               (filter_percentile_rowsums,
                sum(take),
                len(take) - sum(take)))
        counts_table = counts_table[take]

    return counts_table


def identifyVariablesPandas(design_table):

    # design table should have been processed by loadTagDataPandas already
    # just in case, re-filter for not included samples here

    design_table = design_table[design_table["include"] != 0]
    conds = design_table['group'].tolist()
    pairs = design_table['pair'].tolist()

    # TS, adapted from JJ code for DESeq2 design tables:
    # if additional columns present, pass to 'factors'
    if len(design_table.columns) > 3:
        factors = design_table.iloc[:, 3:]
    else:
        factors = None

    return conds, pairs, factors


def checkTagGroupsPandas(design_table, ref_group=None):
    '''compute groups and pairs from tag data table.'''

    conds, pairs, factors = identifyVariablesPandas(design_table)
    groups = list(set(conds))

    # Relevel the groups so that the reference comes first
    # how to do this in python?
    # if ref_group is not None:
    #    R('''groups <- relevel(groups, ref = "%s")''' % ref_group)

    # check this works, will need to make factors from normal df
    # TS adapted from JJ code for DESeq2 -
    # check whether there are additional factors in design file...
    if factors:
        E.warn("There are additional factors in design file that are ignored"
               " by groupTagData: ", factors)
    else:
        pass

    # Test if replicates exist - at least one group must have multiple samples
    max_per_group = max([conds.count(x) for x in groups])

    has_replicates = max_per_group >= 2

    # Test if pairs exist:
    npairs = len(set(pairs))
    has_pairs = npairs == 2

    # ..if so, at least two samples are required per pair
    if has_pairs:
        min_per_pair = min([pairs.count(x) for x in set(pairs)])
        has_pairs = min_per_pair >= 2

    return groups, pairs, conds, factors, has_replicates, has_pairs


ResultColumns = ["test_id", "treatment_name", "treatment_mean",
                 "treatment_std", "control_name", "control_mean",
                 "control_std", "p_value", "p_value_adj", "l2fold", "fold",
                 "transformed_l2fold", "significant", "status"]

ResultColumns_dtype = {"test_id": object, "treatment_name": object,
                       "treatment_mean": float, "treatment_std":
                       float, "control_name": object, "control_mean":
                       float, "control_std": float, "p_value": float,
                       "p_value_adj": float, "l2fold": float, "fold":
                       float, "transformed_l2fold": float,
                       "significant": int, "status": object}


def makeEmptyDataFrameDict():
    return {key: [] for key in ResultColumns}


def runTTestPandas(counts_table,
                   design_table,
                   outfile,
                   outfile_prefix,
                   fdr,
                   ref_group=None):
    '''apply a ttest on the data.

    For the T-test it is best to use FPKM values as
    this method does not perform any library normalization.
    Alternatively, perform normalisation on counts table using Counts.py
    '''
    stats = importr('stats')

    (groups, pairs, conds, factors, has_replicates,
     has_pairs) = checkTagGroupsPandas(design_table, ref_group)

    df_dict = makeEmptyDataFrameDict()

    for combination in itertools.combinations(groups, 2):
        # as each combination may have different numbers of samples in control
        # and treatment, calculations have to be performed on a per
        # combination basis

        control, treatment = combination
        n_rows = counts_table.shape[0]
        df_dict["control_name"].extend((control,)*n_rows)
        df_dict["treatment_name"].extend((treatment,)*n_rows)
        df_dict["test_id"].extend(counts_table.index.tolist())

        # subset counts table for each combination
        c_keep = [x == control for x in conds]
        control_counts = counts_table.iloc[:, c_keep]
        t_keep = [x == treatment for x in conds]
        treatment_counts = counts_table.iloc[:, t_keep]

        c_mean = control_counts.mean(axis=1)
        df_dict["control_mean"].extend(c_mean)
        df_dict["control_std"].extend(control_counts.std(axis=1))

        t_mean = treatment_counts.mean(axis=1)
        df_dict["treatment_mean"].extend(t_mean)
        df_dict["treatment_std"].extend(treatment_counts.std(axis=1))

        t, prob = ttest_ind(control_counts, treatment_counts, axis=1)
        df_dict["p_value"].extend(prob)
        # what about zero values?!
        df_dict["fold"].extend(t_mean / c_mean)

    df_dict["p_value_adj"].extend(
        list(stats.p_adjust(FloatVector(df_dict["p_value"]), method='BH')))

    df_dict["significant"].extend(
        [int(x < fdr) for x in df_dict["p_value_adj"]])

    df_dict["l2fold"].extend(list(numpy.log2(df_dict["fold"])))
    # note: the transformed log2 fold change is not transformed!
    df_dict["transformed_l2fold"].extend(list(numpy.log2(df_dict["fold"])))

    # set all status values to "OK"
    df_dict["status"].extend(("OK",)*n_rows)

    results = pandas.DataFrame(df_dict)
    results.set_index("test_id", inplace=True)
    results.to_csv(outfile, sep="\t", header=True, index=True)


def plotCorrelationHeatmapMatplot(counts, outfile, method="correlation",
                                  cor_method="pearson"):
    '''plot a heatmap of correlations derived from
    countsTable.
    '''
    # to do: add other methods?

    # define outside function? - Will we reuse?
    heatmap_cdict_b_to_y = {
        'red': ((0.0, 0.4, .4), (0.01, .4, .4), (1., .95, .95)),
        'green': ((0.0, 0.4, 0.4), (0.01, .4, .4), (1., .95, .95)),
        'blue': ((0.0, .9, .9), (0.01, .9, .9), (1., 0.4, 0.4))}

    cm = matplotlib.colors.LinearSegmentedColormap(
        '', heatmap_cdict_b_to_y, 256)

    df = counts.corr(method=cor_method)
    plt.pcolor(np.array(df), cmap=cm)
    plt.colorbar()
    plt.title("%(cor_method)s correlation heatmap" % locals())
    plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
    plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
    plt.tight_layout()
    pylab.savefig(outfile)


def runEdgeRPandas(counts,
                   design_table,
                   outfile,
                   outfile_prefix="edger.",
                   fdr=0.1,
                   prefix="",
                   dispersion=None,
                   ref_group=None):
    '''run EdgeR on countsTable.

    Results are stored in *outfile* and files prefixed by *outfile_prefix*.

    The dispersion is usually measuered from replicates. If there are no
    replicates, you need to set the *dispersion* explicitely.

    See page 13 of the EdgeR user guide::

       2. Simply pick a reasonable dispersion value, based on your
       experience with similar data, and use that. Although
       subjective, this is still more defensible than assuming Poisson
       variation. Typical values are dispersion=0.4 for human data,
       dispersion=0.1 for data on genetically identical model
       organisms or dispersion=0.01 for technical replicates.

    '''

    # load library
    R('''suppressMessages(library('edgeR'))''')

    (groups, pairs, conds, factors, has_replicates,
     has_pairs) = checkTagGroupsPandas(design_table, ref_group)

    if not has_replicates and dispersion is None:
        raise ValueError("no replicates and no dispersion")

    # output heatmap plot
    plotCorrelationHeatmapMatplot(counts,
                                  '%(outfile_prefix)sheatmap.png' % locals(),
                                  cor_method="spearman")

    E.info('running EdgeR: groups=%s, pairs=%s, replicates=%s, pairs=%s' %
           (groups, pairs, has_replicates, has_pairs))

    r_counts = pandas2ri.py2ri(counts)

    passDFtoRGlobalEnvironment = R('''function(df){
    countsTable <<- df}''')
    passDFtoRGlobalEnvironment(r_counts)

    if has_pairs:
        # output difference between groups
        # TS #####
        # this is performed on non-normalised data
        # should we use Counts.py to normalise first?
        # also, this isn't edgeR specific, should this be
        # moved to a seperate summary function?
        # also move the MDS plotting?
        # #####
        first = True
        pairs_df = pandas.DataFrame()
        nrows = len(counts.index)
        n = 0
        for g1, g2 in itertools.combinations(groups, 2):
            print g1, g2
            keep_a = [x == g1 for x in conds]
            counts_a = counts.iloc[:, keep_a]
            keep_b = [x == g2 for x in conds]
            counts_b = counts.iloc[:, keep_b]
            index = range(n, n+nrows)
            n += nrows
            a = counts_a.sum(axis=1)
            b = counts_b.sum(axis=1)
            diff = a-b
            diff.sort()
            temp_df = pandas.DataFrame({"cumsum": np.cumsum(diff).tolist(),
                                        "comb": "_vs_".join([g1, g2]),
                                        "id": range(0, nrows)},
                                       index=index)
            pairs_df = pairs_df.append(temp_df)
        plot_pairs = R('''function(df, outfile){
        suppressMessages(library('ggplot2'))
        p = ggplot(df, aes(y=cumsum, x=id)) +
        geom_line(aes(col=factor(comb))) +
        scale_color_discrete(name="Comparison") +
        xlab("index") + ylab("Cumulative sum")
        ggsave("%(outfile_prefix)sbalance_groups.png", plot = p)}
        ''' % locals())

        r_pairs_df = pandas2ri.py2ri(pairs_df)
        plot_pairs(r_pairs_df)

        # output difference between pairs within groups
        first = True
        legend = []
        n = 0
        pairs_in_groups_df = pandas.DataFrame()
        for pair in set(pairs):
            for g1, g2 in itertools.combinations(groups, 2):
                print pair, g1, g2
                key = "pair-%s-%s-vs-%s" % (pair, g1, g2)
                legend.append(key)
                keep_a = [x == pair and y == g1 for x, y in zip(pairs, conds)]
                keep_b = [x == pair and y == g2 for x, y in zip(pairs, conds)]

                # check pair and group combination actually exists
                if sum(keep_a) > 0 and sum(keep_b) > 0:
                    counts_a = counts.iloc[:, keep_a]
                    counts_b = counts.iloc[:, keep_b]
                    index = range(n, n+nrows)
                    n += nrows
                    a = counts_a.sum(axis=1)
                    b = counts_b.sum(axis=1)
                    diff = a-b
                    diff.sort()
                    comparison = "pair-%s-%s-vs-%s" % (pair, g1, g2)
                    temp_df = pandas.DataFrame({"cumsum": np.cumsum(diff).tolist(),
                                                "comb": comparison,
                                                "id": range(0, nrows)},
                                               index=index)
                    pairs_in_groups_df = pairs_in_groups_df.append(temp_df)

        plot_pairs_in_groups = R('''function(df, outfile){
        suppressMessages(library('ggplot2'))
        p = ggplot(df, aes(y=cumsum, x=id)) +
        geom_line(aes(col=factor(comb))) +
        scale_color_discrete(name="Comparison") +
        xlab("index") + ylab("Cumulative sum")
        ggsave("%(outfile_prefix)sbalance_pairs.png", plot = p)}
        ''' % locals())

        r_pairs_in_groups_df = pandas2ri.py2ri(pairs_in_groups_df)
        plot_pairs_in_groups(r_pairs_in_groups_df)

    # create r objects
    r_counts = pandas2ri.py2ri(counts)
    r_groups = ro.StrVector(conds)
    r_pairs = ro.StrVector(pairs)
    r_has_pairs = ro.default_py2ri(has_pairs)
    r_has_replicates = ro.default_py2ri(has_replicates)

    if dispersion is not None:
        r_dispersion = ro.default_py2ri(dispersion)
    else:
        r_dispersion = ro.default_py2ri(False)

    if ref_group is not None:
        r_ref_group = ro.default_py2ri(ref_group)
    else:
        r_ref_group = ro.default_py2ri(groups[0])

    # build DGEList object
    buildDGEList = R('''
    suppressMessages(library('edgeR'))
    function(counts, groups, ref_group){
    countsTable = DGEList(counts, group=groups)
    countsTable$samples$group <- relevel(countsTable$samples$group,
    ref = ref_group)
    countsTable = calcNormFactors(countsTable)
    return(countsTable)
    }''')

    r_countsTable = buildDGEList(r_counts, r_groups, r_ref_group)

    # output MDS plot
    # TT - should this be performed on the normalised counts table?(see above)
    R.png('''%(outfile_prefix)smds.png''' % locals())
    try:
        MDSplot = R('''function(counts){
        plotMDS(counts)}''')
        MDSplot(r_counts)
    except rpy2.rinterface.RRuntimeError, msg:
        E.warn("can not plot mds: %s" % msg)
    R['dev.off']()

    # build design matrix
    buildDesign = R('''function(countsTable, has_pairs, pairs){
    if (has_pairs==TRUE) {
        design <- model.matrix( ~pairs + countsTable$samples$group ) }
    else {
        design <- model.matrix( ~countsTable$samples$group ) }
    return(design)
    }''')

    r_design = buildDesign(r_countsTable, r_has_pairs, r_pairs)

    fitModel = R('''function(countsTable, design, has_replicates, dispersion){
    if (has_replicates == TRUE) {
        # estimate common dispersion
        countsTable = estimateGLMCommonDisp( countsTable, design )
        # estimate tagwise dispersion
        countsTable = estimateGLMTagwiseDisp( countsTable, design )
        # fitting model to each tag
        fit = glmFit( countsTable, design ) }
    else {
        # fitting model to each tag
        fit = glmFit(countsTable, design, dispersion=dispersion) }
    return(fit)
    }''')

    r_fit = fitModel(r_countsTable, r_design, r_has_replicates, r_dispersion)

    E.info("Generating output")

    # perform LR test
    lrtTest = R('''function(fit, prefix){
    lrt = glmLRT(fit)
    # save image for access to the whole of the lrt object
    save.image(paste0(prefix,"lrt.RData"))
    return(lrt)
    }''')
    r_lrt = lrtTest(r_fit, outfile_prefix)

    # return statistics table - must be a better way to do this?
    extractTable = R('''function(lrt){ return(lrt$table)}''')
    r_lrt_table = extractTable(r_lrt)

    # output cpm table
    outputCPMTable = R('''function(countsTable, outfile_prefix, lrt){
    suppressMessages(library(reshape2))
    countsTable.cpm <- cpm(countsTable, normalized.lib.sizes=TRUE)
    melted <- melt(countsTable.cpm)
    names(melted) <- c("test_id", "sample", "ncpm")

    # melt columns are factors - convert to string for sorting
    melted$test_id = levels(melted$test_id)[as.numeric(melted$test_id)]
    melted$sample = levels(melted$sample)[as.numeric(melted$sample)]

    # sort cpm table by test_id and sample
    sorted <- melted[with(melted, order(test_id, sample)),]
    gz <- gzfile(paste0(outfile_prefix,"cpm.tsv.gz"), "w" )
    write.table(sorted, file=gz, sep = "\t", row.names=FALSE, quote=FALSE)
    close(gz)}''')

    outputCPMTable(r_countsTable, outfile_prefix, r_lrt)

    # output differences between pairs
    if len(groups) == 2:
        plotSmear = R('''function(countsTable, outfile_prefix){
        png(paste0(outfile_prefix,"smaplot.png"
        plotSmear(countsTable, pair=c('%s'))
        abline(h=c(-2, 2), col='dodgerblue')
        dev.off()}''' % "','".join(groups))

    lrt_table = pandas2ri.ri2py(r_lrt_table)

    n_rows = lrt_table.shape[0]
    df_dict = makeEmptyDataFrameDict()

    # import r stats module for BH adjustment
    stats = importr('stats')

    df_dict["control_name"].extend((groups[0],)*n_rows)
    df_dict["treatment_name"].extend((groups[1],)*n_rows)
    df_dict["test_id"].extend(lrt_table.index)
    df_dict["control_mean"].extend(lrt_table['logCPM'])
    df_dict["treatment_mean"].extend(lrt_table['logCPM'])
    df_dict["control_std"].extend((0,)*n_rows)
    df_dict["treatment_std"].extend((0,)*n_rows)
    df_dict["p_value"].extend(lrt_table['PValue'])
    df_dict["p_value_adj"].extend(
        list(stats.p_adjust(FloatVector(df_dict["p_value"]), method='BH')))
    df_dict["significant"].extend(
        [int(float(x) < fdr) for x in df_dict["p_value_adj"]])
    df_dict["l2fold"].extend(list(numpy.log2(lrt_table['logFC'])))

    # TS -note: the transformed log2 fold change is not transformed!
    df_dict["transformed_l2fold"].extend(list(numpy.log2(lrt_table['logFC'])))

    # TS -note: check what happens when no fold change is available
    # TS -may need an if/else in list comprehension. Raise E.warn too?
    df_dict["fold"].extend([math.pow(2, float(x)) for x in lrt_table['logFC']])

    # set all status values to "OK"
    # TS - again, may need an if/else, check...
    df_dict["status"].extend(("OK",)*n_rows)

    results = pandas.DataFrame(df_dict)
    results.set_index("test_id", inplace=True)
    results.to_csv(outfile, sep="\t", header=True, index=True)

    counts = E.Counter()
    counts.signficant = sum(results['significant'])
    counts.insignficant = (len(results['significant']) - counts.signficant)
    counts.all_over = sum([x > 0 for x in results['l2fold']])
    counts.all_under = sum([x < 0 for x in results['l2fold']])
    counts.signficant_over = sum([results['significant'][x] == 1 and
                                 results['l2fold'][x] > 1 for
                                 x in range(0, n_rows)])
    counts.signficant_under = sum([results['significant'][x] == 1 and
                                   results['l2fold'][x] < 1 for
                                   x in range(0, n_rows)])

    outf = IOTools.openFile("%(outfile_prefix)ssummary.tsv" % locals(), "w")
    outf.write("category\tcounts\n%s\n" % counts.asTable())
    outf.close()


def outputSpikeIns(filename_tags,
                   outfile,
                   output_filename_pattern,
                   filename_design=None,
                   foldchange_max=10.0,
                   expression_max=5.0,
                   max_counts_per_bin=100,
                   expression_bin_width=0.5,
                   foldchange_bin_width=0.5,
                   iterations=1):

    E.info("loading tag data from %s" % filename_tags)

    if filename_design is not None:
        # load all tag data
        counts_table, design_table = loadTagDataPandas(
            filename_tags, filename_design)

        # filter
        counts_table = filterTagDataPandas(counts_table, design_table)

    else:
        raise NotImplementedError("reading full table not implemented")

    nobservations, nsamples = counts_table.shape

    groups = list(set(design_table["group"]))
    if len(groups) != 2:
        raise NotImplementedError(
            "spike in only implemented for one pairwise comparison")

    # select group data
    group1 = counts_table[design_table[design_table.group == groups[0]].index]
    group2 = counts_table[design_table[design_table.group == groups[1]].index]

    outfile.write("interval_id\t%s\t%s\n" %
                  ("\t".join(group1.columns), "\t".join(group2.columns)))
    outf_info = IOTools.openFile(output_filename_pattern + "info.tsv.gz", "w")
    outf_info.write("interval_id\tl10_expression\tl2fold\n")

    # important: convert to matrixes, otherwise there will be a per-name lookup
    # when l10average or l2fold are computed.
    group1 = group1.as_matrix()
    group2 = group2.as_matrix()

    # compute bins
    expression_bins = numpy.arange(0, expression_max, expression_bin_width)
    fold_change_bins = numpy.arange(-foldchange_max,
                                    foldchange_max, foldchange_bin_width)

    E.info("l10expression bins=%s" % expression_bins)
    E.info("l2fold change bins=%s" % fold_change_bins)

    # output values
    output_counts = numpy.zeros(
        (len(expression_bins) + 1, len(fold_change_bins) + 1))

    for iteration in range(iterations):

        # randomize order
        group1 = numpy.random.permutation(group1)
        group2 = numpy.random.permutation(group2)

        # compute means and foldchanges
        group1_mean = group1.mean(1)
        group2_mean = group2.mean(1)

        # compute average expression by means
        l10average = numpy.log((group1_mean + group2_mean) / 2.0)

        # compute a fold change with pseudo count of 1
        l2fold = numpy.log2((group1_mean + 1.0) / (group2_mean + 1.0))

        # digitize arrays with bins
        l10average_idx = numpy.digitize(l10average, expression_bins)
        l2fold_idx = numpy.digitize(l2fold, fold_change_bins)

        interval_id = 0
        for idx, coord in enumerate(zip(l10average_idx, l2fold_idx)):
            # assert expression_bins[coord[0]-1] <= l10average[idx] < expression_bins[coord[0]]
            # assert fold_change_bins[coord[1]-1] <= l2fold[idx] <  expression_bins[coord[1]]

            if output_counts[coord] >= max_counts_per_bin:
                continue

            output_counts[coord] += 1
            outf_info.write("spike%i\t%f\t%f\n" % (interval_id,
                                                   l10average[idx],
                                                   l2fold[idx]))

            outfile.write("spike%i\t%s\t%s\n" %
                          (interval_id,
                           "\t".join(
                               map(str, list(group1[idx]))),
                           "\t".join(map(str, list(group2[idx])))))
            interval_id += 1

    outf_info.close()

    df = pandas.DataFrame(output_counts)
    df.index = list(expression_bins) + [expression_max]
    df.columns = list(fold_change_bins) + [foldchange_max]

    df.to_csv(IOTools.openFile(output_filename_pattern + "counts.tsv.gz", "w"),
              sep="\t")

    E.info("output %i spike in intervals" % interval_id)
