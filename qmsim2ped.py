#!/usr/bin/env python

##########
#
#                               reHC-*
#  Haplotyping with Recombinations, Errors, and Missing Genotypes
#
#  Copyright (C) 2010,2011,2012  Yuri Pirola <yuri.pirola(-at-)gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#
#  This file is part of reHC-* (reHCstar).
#
#  reHC-* is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  reHC-* is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with reHC-*.  If not, see <http://www.gnu.org/licenses/>.
#
##########

##########
#
#  qmsim2ped.py
#
#  A program that converts a genotyped and structured population generated
#  by QMSim to the PED plink format.
#
##########

from __future__ import print_function

import bisect
import fnmatch
import logging
import math
import optparse
import os
import random
import re
import sys


def parse_command_line():
    usage= "usage: %prog [options]"
    parser= optparse.OptionParser(usage=usage)
    parser.add_option("-t", "--output_filename_template",
                      action="store", dest="out_templ",
                      type="string", default=
                      "gen-ped-size{size:03d}-conf{conf:01d}-length{length:04d}-"
                      "miss{miss:.3f}-recUNKNOWN-err{err:.3f}-hc{hc:01d}#txt",
                      help="The template used to compute the output filenames"
                      " (default: %default)",
                      metavar="FILENAME TEMPLATE")
    parser.add_option("-m", "--missing_genotype_probability",
                      action="store", dest="missing",
                      type="float", default=0.0,
                      help="the probability that a single-locus genotype has not been called"
                      " (default: %default)",
                      metavar="PROBABILITY [0,1]")
    parser.add_option("-e", "--error_probability",
                      action="store", dest="error",
                      type="float", default=0.0,
                      help="the probability that a single-locus genotype has been mis-called"
                      " (default: %default)",
                      metavar="PROBABILITY [0,1]")
    parser.add_option("-s", "--seed",
                      action="store", dest="seed",
                      type="int", default=122295,
                      help="the seed of the random generator"
                      " (default: %default)",
                      metavar="INT")
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      default=False,
                      help="print additional log messages")
    (options, args) = parser.parse_args()
    return options

options = parse_command_line()

missing_genotype_prob = options.missing
error_prob = options.error
seed = options.seed
homo1 = "1 1"
homo2 = "2 2"
heter = "1 2"
missing = "0 0"
hhomo1 = "1|1"
hhomo2 = "2|2"
hheter1 = "1|2"
hheter2 = "2|1"
genotype_mapping = { "0": homo1, "2": homo2, "3": heter, "4": heter, "5": missing }
haplotype_mapping = { "0": hhomo1, "2": hhomo2, "3": hheter1, "4": hheter2 }  ## No missing here!
others = { homo1: (homo2, heter),
           homo2: (homo1, heter),
           heter: (homo1, homo2) }

pheno= "phenotype"

random.seed(seed)

pedfile = { "pref": "pop1_data_",
            "skip_rows": 4 }
genfile = { "pref": "pop1_mrk_",
            "skip_rows": 4 }
mrkfile = { "pref": "lm_mrk_",
            "skip_rows": 6 }
xoverfile = { "pref": "pop1_crosso_",
              "skip_rows": 8 }

log_level = logging.DEBUG if options.verbose else logging.INFO
logging.basicConfig(level=log_level,
                    format='%(levelname)-7s [%(asctime)s]  %(message)s',
                    datefmt="%y%m%d %H%M%S")

logging.info("QMSim to PED converter -- started")
logging.info("Configuration: %s", options)
logging.info("Execution in directory: %s", os.getcwdu())
logging.info("Command line: %s", " ".join(sys.argv))

files_in_curdir = os.listdir(os.curdir)

pedfiles = sorted(fnmatch.filter(files_in_curdir, pedfile["pref"]+"*"))
genfiles = sorted(fnmatch.filter(files_in_curdir, genfile["pref"]+"*"))
mrkfiles = sorted(fnmatch.filter(files_in_curdir, mrkfile["pref"]+"*"))
xoverfiles = sorted(fnmatch.filter(files_in_curdir, xoverfile["pref"]+"*"))

assert len(pedfiles)==len(genfiles) and len(genfiles)==len(mrkfiles) and len(mrkfiles)==len(xoverfiles)



for counter, popfiles in enumerate(zip(pedfiles, genfiles, mrkfiles, xoverfiles)):
    logging.info("Deriving pedigree %d from files (%s)...",
                 counter, ", ".join(popfiles))
    pedigree = []
    logging.debug("Reading pedigree from %s...", popfiles[0])
    with open(popfiles[0], "r") as pedf:
        for rownum,row in enumerate(pedf):
            if rownum < pedfile["skip_rows"]:
                continue
            row = row.strip().split(None, 5)
            assert int(row[0]) == len(pedigree)+1
            row[4] = "1" if row[4] == "M" else "2"
            del row[5]
            del row[3]
            pedigree.append(row)
        logging.debug("Read %d individuals.", len(pedigree))

    genotypes = []
    haplotypes = []
    nloci = None
    logging.debug("Reading genotypes and haplotypes from %s...", popfiles[1])
    with open(popfiles[1], "r") as genf:
        for rownum,row in enumerate(genf):
            if rownum < genfile["skip_rows"]:
                continue
            row = row.strip().split(None, 3)
            assert int(row[0]) == len(genotypes)+1
            hapstr = row[3].split(None)
            gen = [ genotype_mapping[locus] for locus in hapstr ]
            hap = [ haplotype_mapping[locus] for locus in hapstr ]
            assert nloci == None or (nloci == len(gen) and nloci == len(hap))
            nloci = len(gen)
            genotypes.append(gen)
            haplotypes.append(hap)
        logging.debug("Read %d genotypes and haplotype pairs over %d loci.",
                      len(genotypes), nloci)

    mrk_map = []
    chrs = set()
    logging.debug("Reading marker map from %s...", popfiles[2])
    with open(popfiles[2], "r") as mrkf:
        for rownum,row in enumerate(mrkf):
            if rownum < mrkfile["skip_rows"]:
                continue
            row = row.strip().split(None, 2)
            mrk_map.append((int(row[1]), float(row[2])))
            chrs.add(int(row[1]))
        assert nloci == len(mrk_map)
        logging.debug("Read the map of %d loci over %d chromosomes (%s).",
                      len(mrk_map), len(chrs), ", ".join([ str(e) for e in chrs]))

    def mapxovers(mrk_map, xovers):
        xovers = xovers.split(None)
        outxovers = []
        for xover in xovers:
            if xover.startswith("C"):
                cmatch = re.match("C(?P<chrom>\d+){", xover)
                assert cmatch
                chrom = int(cmatch.group("chrom"))
            else:
                assert chrom
                delchrom = xover.endswith("}")
                xover = xover.rstrip("}")
                cmpos = (chrom, float(xover))
                cmidx = bisect.bisect_left(mrk_map, cmpos)
                if (cmidx==0 or cmidx==len(mrk_map) or
                    mrk_map[cmidx][0] != chrom or mrk_map[cmidx-1][0] != chrom):
                    logging.debug("Discarding crossover (%d, %f) outside the marker interval.",
                                  *cmpos)
                else:
                    outxovers.append((chrom, cmidx))
                if delchrom:
                    chrom = None
        return outxovers


    xovers = []
    logging.debug("Reading crossover events from %s...", popfiles[3])
    with open(popfiles[3], "r") as xoverf:
        for rownum,row in enumerate(xoverf):
            if rownum < xoverfile["skip_rows"]:
                continue
            row = row.strip().split(None, 1)
            if row[0] == "I":
                ind = row[1]
                ixovers = []
                xovers.append(ixovers)
            elif row[0] == "S":
                ixovers.extend([ (chrom, ind, pos, 0)
                                 for chrom,pos in mapxovers(mrk_map, row[1]) ])
            elif row[0] == "D":
                ixovers.extend([ (chrom, ind, pos, 1)
                                 for chrom,pos in mapxovers(mrk_map, row[1]) ])
            else:
                sys.exit("SOMETHING IMPOSSIBLE HAPPENED")


    orig_genotypes = genotypes

    logging.debug("Simulating missing genotypes with uniform probability %.5f...",
                  missing_genotype_prob)
    missings = [ [ 0 if (random.random() >= missing_genotype_prob) else 1
                   for j,g in enumerate(gen) ]
               for i,gen in enumerate(genotypes) ]
    genotypes = [ [ g if m==0 else missing
                    for m,g in zip(mis, gen) ]
                  for mis,gen in zip(missings, genotypes) ]



    logging.debug("Splitting genotypes into chromosomes...")
    def split_into_chromosomes(vectors, chrs, mrk_map):
        newvectors= [ [ ] for chrom in chrs ]
        for vect in vectors:
            for chrom in chrs:
                newvectors[chrom-1].append([])
            for mm,v in zip(mrk_map, vect):
                newvectors[mm[0]-1][-1].append(v)
        return newvectors

    orig_genotypes = split_into_chromosomes(orig_genotypes, chrs, mrk_map)
    genotypes = split_into_chromosomes(genotypes, chrs, mrk_map)
    haplotypes = split_into_chromosomes(haplotypes, chrs, mrk_map)
    missings = split_into_chromosomes(missings, chrs, mrk_map)
    errors = {}

    # Output
    for chrom in chrs:
        logging.info("Considering chromosome/haplotype configuration %d...", chrom)
        chromidx = chrom - 1
        first_locus = min([i for i,mm in enumerate(mrk_map) if mm[0]==chrom])
        nloci = len([0 for mm in mrk_map if mm[0]==chrom])

		# Simulate a fixed amount of errors
        all_positions = [ (i,j)
                          for i,mis in enumerate(missings[chromidx])
                          for j,m in enumerate(mis) if m==0 ]
        no_of_errors = int(math.ceil(len(all_positions)*error_prob))
        logging.debug("Simulating %d genotyping errors... (error rate=%.5f)",
                      no_of_errors, error_prob)
        error_positions = frozenset(random.sample(all_positions, no_of_errors))
        errors[chromidx] = [ [ 1 if (i,j) in error_positions else 0
                               for j,g in enumerate(gen) ]
                             for i,gen in enumerate(genotypes[chromidx]) ]
        genotypes[chromidx] = [ [ g if e==0 else random.choice(others[g])
                                  for e,g in zip(err, gen) ]
                                for err,gen in zip(errors[chromidx], genotypes[chromidx]) ]

        no_of_errors = [ sum(err) for err in errors[chromidx] ]
        logging.info("Simulated %d errors (max per individual: %d).",
                     sum(no_of_errors), max(no_of_errors))

        no_of_missings = [ sum(mis) for mis in missings[chromidx] ]
        logging.info("Simulated %d missing genotypes (max per individual: %d).",
                     sum(no_of_missings), max(no_of_missings))
        outfilename = options.out_templ.format(conf=counter, size=len(genotypes[chromidx]),
                                               length=nloci,
                                               miss=missing_genotype_prob, err=error_prob,
                                               hc=chromidx)
        outfilename = outfilename.replace(".", "_").replace("#", ".")
        logging.info("Writing genotyped pedigree to file %s...", outfilename)
        thisnloci = len(genotypes[chromidx][0])
        with open(outfilename, "w") as outf:
            # Header
            print("### Genotyped pedigree generated from QMSim files ({0})"
                  "".format(", ".join(popfiles)),
                  file=outf)
            print("### Missing genotype probability = {0:.5f}"
                  "".format(missing_genotype_prob), file=outf)
            print("### Genotyping error probability = {0:.5f}".format(error_prob), file=outf)
            print("### Seed = {0:d}".format(seed), file=outf)
            ## Errors
            print("##ERRORS_IN_A_INDIVIDUAL\t{0}\t{1}"
                  "".format(max(no_of_errors),
                            max(no_of_errors)/float(thisnloci)),
                  file=outf)
            print("##TOTAL_ERRORS\t{0}\t{1}"
                  "".format(sum(no_of_errors),
                            sum(no_of_errors)/float(thisnloci*len(genotypes[chromidx]))),
                  file=outf)
            if sum(no_of_errors)>0:
                error_list = [ "ERR {individual}\t{locus}\t({original}) -> ({miscalled})"
                               "".format(individual=pedigree[ind][0],
                                         locus=locus,
                                         original=orig_genotypes[chromidx][ind][locus],
                                         miscalled=genotypes[chromidx][ind][locus])
                               for ind,err in enumerate(errors[chromidx])
                               for locus,e in enumerate(err)
                               if e == 1 ]
                print("## ERR individual locus (original_genotype) -> (miscalled genotype)",
                      file=outf)
                print("\n".join(["# {0}".format(f) for f in error_list]),
                      file=outf)
                del error_list

            ## Recombinations
            recomb_list = [ "REC {individual} {locus} {parent}"
                            "".format(individual=xover[1],
                                      locus=xover[2]-first_locus,
                                      parent=xover[3])
                            for xoveri in xovers
                            for xover in xoveri
                            if xover[0] == int(chrom) ]
            if recomb_list:
                print("# Computed from crossover event # REC individual locus 0/1 "
                      "(0=on paternal haplotype, 1=on maternal haplotype)",
                      file=outf)
                print("\n".join(["# Computed from crossover event # {0}"
                                 "".format(f) for f in set(recomb_list)]),
                      file=outf)

            del recomb_list


            ## Haplotypes
            for ind,hap in zip(pedigree,haplotypes[chromidx]):
                print("# GENERATED_HAPLOTYPES\t0\t{indped}\tphenotype\t{haplotypes}"
                      "".format(indped="\t".join(ind),
                                haplotypes="\t".join(hap)),
                      file=outf)


            ## Genotypes
            for ind,gen in zip(pedigree,genotypes[chromidx]):
                print("0\t{indped}\tphenotype\t{genotype}"
                      "".format(indped="\t".join(ind),
                                genotype="\t".join(gen)),
                      file=outf)


logging.info("QMSim to PED converter -- terminated")
