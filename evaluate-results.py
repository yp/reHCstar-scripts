#!/usr/bin/python3


##########
#
#                               reHC-*
#  Haplotyping with Recombinations, Errors, and Missing Genotypes
#
#  Copyright (C) 2010,2011  Yuri Pirola <yuri.pirola(-at-)gmail.com>
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
#  evaluate-results.py
#
#  A program that evaluates the quality of the results computed by reHC-*
#  compared to the generated haplotype configuration.
#
##########

import os
import sys
import random
import logging
import operator
import itertools
import functools
from collections import Counter
from optparse import OptionParser

def parse_command_line():
    usage= "usage: %prog [options]"
    parser= OptionParser(usage=usage)
    parser.add_option("-o", "--original",
                      action="store", dest="original",
                      type="string", default=None,
                      help="the file containing the genotyped pedigree and the "
                      "original haplotype configuration (as produced by 'hc-generator.py')",
                      metavar="FILE")
    parser.add_option("-r", "--result",
                      action="store", dest="result",
                      type="string", default=None,
                      help="the file containing the haplotype configuration computed "
                      "by reHC-*",
                      metavar="FILE")
    parser.add_option("-f", "--full-stats",
                      action="store_true", dest="full",
                      default=False,
                      help="print statistics for each individual")
    parser.add_option("-H", "--show-header",
                      action="store_true", dest="header",
                      default=False,
                      help="print the header")
    parser.add_option("-n", "--dont-normalize-founders",
                      action="store_true", dest="no_norm_found",
                      default=False,
                      help="do not swap founders' haplotypes such that the minimum amount of "
                      "phase error is obtained")
    parser.add_option("-p", "--normalize-pseudofounders",
                      action="store_true", dest="norm_pseudo",
                      default=False,
                      help="swap haplotypes of individuals whose parents are not genotyped "
                      "such that the minimum amount of phase error is obtained")
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      default=False,
                      help="print additional log messages")
    (options, args) = parser.parse_args()

    return(options)

def list_el(v):
    return [ l for l,e in zip(itertools.count(0), v) if e==1 ]


def compute_errors(genotype, orig_hp, orig_hm, res_hp, res_hm):
    orig_gen= [ encoding[str(hp) + ' ' + str(hm)] for hp,hm in zip(orig_hp, orig_hm) ]
    res_gen= [ encoding[str(hp) + ' ' + str(hm)] for hp,hm in zip(res_hp, res_hm) ]
    gen_dif= [ 0 if (g != 0 or orig == 0 or res == 0 or orig == res) else 1
               for g,orig,res in zip(genotype, orig_gen, res_gen) ]
    pat_dif= [ 0 if orig == 0 or orig == res else 1
               for orig,res in zip(orig_hp, res_hp) ]
    mat_dif= [ 0 if orig == 0 or orig == res else 1
               for orig,res in zip(orig_hm, res_hm) ]
    mask_pat_dif= [ 0 if g == 0 or h == 0 else diff
                    for diff, g, h in zip(pat_dif, genotype, orig_hp) ]
    mask_mat_dif= [ 0 if g == 0 or h == 0 else diff
                    for diff, g, h in zip(mat_dif, genotype, orig_hm) ]
    err_gen= sum(gen_dif)
    err_pat_hap= sum(pat_dif)
    err_mat_hap= sum(mat_dif)
    err_mask_pat_hap= sum(mask_pat_dif)
    err_mask_mat_hap= sum(mask_mat_dif)

    return ( err_gen,
             err_pat_hap, err_mat_hap,
             err_mask_pat_hap, err_mask_mat_hap,
             [ list_el(l) for l in (gen_dif, pat_dif, mat_dif,
                                    mask_pat_dif, mask_mat_dif) ] )

def compute_errors_with_swap(orig_ped, genotype, orig_hp, orig_hm, res_hp, res_hm, norm_found, norm_pseudofound):
    (father, mother)= orig_ped[individual][1:3]
    if ( ( (norm_found) and (father == '0') and (mother == '0') ) or
         ( (norm_pseudofound) and
           (father != '0') and (father not in orig_ped) and
           (mother != '0') and (mother not in orig_ped) ) ):
        ris1= compute_errors(genotype, orig_hp, orig_hm, res_hp, res_hm)
        ris2= compute_errors(genotype, orig_hp, orig_hm, res_hm, res_hp)
        return ris1 if sum(ris1[0:-1]) < sum(ris2[0:-1]) else ris2
    else:
        return compute_errors(genotype, orig_hp, orig_hm, res_hp, res_hm)

def compute_recombinations(ind_id, p_id, h, ph1, ph2):
    lh= len(h)
    rec1= []
    rec2= []
    (phase1, phase2)= (0, 1)
    ph= (ph1, ph2)
    for (l, a, a1, a2) in zip(range(lh), h, ph1, ph2):
        if a == 0 or a1 == 0 or a2 == 0:
            continue
        av= (a1, a2)
        if a != av[phase1]:
            phase1= 1-phase1
            assert(a == av[phase1])
            rec1.append((ind_id, p_id, l))
        if a != av[phase2]:
            phase2= 1-phase2
            assert(a == av[phase2])
            rec2.append((ind_id, p_id, l))
    if len(rec1) <= len(rec2):
        return rec1
    else:
        return rec2



options = parse_command_line()

log_level= logging.DEBUG if options.verbose else logging.INFO

logging.basicConfig(level=log_level,
                    format='%(levelname)-6s [%(asctime)s]  %(message)s')

logging.info("EVALUATION OF (r,e)-HC-* RESULTS")

if ( not options.original or
     not os.path.isfile(options.original) or
     not os.access(options.original, os.R_OK) ):
    logging.fatal("Original file not specified or not accessible. Given '%s'",
                  options.original)
    sys.exit(1)

logging.info("Original file: '%s'", options.original)


if ( not options.result or
     not os.path.isfile(options.result) or
     not os.access(options.result, os.R_OK) ):
    logging.fatal("Result file not specified or not accessible. Given '%s'",
                  options.result)
    sys.exit(1)

logging.info("Result file:   '%s'", options.result)



# Read the original haplotype configuration
logging.info("Reading original file '%s'...", options.original)
missing_str= []
with open(options.original, 'r') as orig:
    missing_str= [ line.strip()
                    for line in orig
                    if not line.startswith('#')  ]
logging.debug("Parsing the original genotype configuration "
              "to get the genotypes...")
genotypes= {}
encoding= { '0 0':0,
            '1 1':1,
            '2 2':2,
            '1 2':3,
            '2 1':3
            }
for r in missing_str:
    split_r= r.split("\t", 6)
    genotypes[ split_r[1] ]= [ encoding[single]
                                  for single in split_r[6].split("\t") ]
del missing_str

logging.info("Re-reading original file '%s'...", options.original)
orig_ped_str= []
with open(options.original, 'r') as orig:
    orig_ped_str= [ line.replace('# GENERATED_HAPLOTYPES', '', 1).strip()
                    for line in orig
                    if line.startswith('# GENERATED_HAPLOTYPES')  ]
logging.info("Parsing the original haplotype configuration...")
orig_ped= {}
for r in orig_ped_str:
    split_r= r.split("\t", 6)
    comb_hap= split_r[6].split("\t")
    pat_hap= [ int(single.split('|')[0]) for single in comb_hap ]
    mat_hap= [ int(single.split('|')[1]) for single in comb_hap ]
    orig_gen= [ encoding[str(hp) + ' ' + str(hm)] for hp,hm in zip(pat_hap, mat_hap) ]
    orig_ped[split_r[1]]= [split_r[1], split_r[2], split_r[3], pat_hap, mat_hap, orig_gen]
del orig_ped_str;

# Read the computed haplotype configuration
logging.info("Reading result file '%s'...", options.result)
res_ped_str= []
with open(options.result, 'r') as result:
    res_ped_str= [ line.strip()
                   for line in result
                   if not line.startswith('#')  ]

logging.debug("Parsing the computed haplotype configuration...")
res_ped= {}
for r in res_ped_str:
    split_r= r.split("\t", 6)
    comb_hap= split_r[6].split("\t")
    pat_hap= [ int(single.split('|')[0]) for single in comb_hap ]
    mat_hap= [ int(single.split('|')[1]) for single in comb_hap ]
    res_gen= [ encoding[str(hp) + ' ' + str(hm)] for hp,hm in zip(pat_hap, mat_hap) ]
    res_ped[split_r[1]]= [split_r[1], split_r[2], split_r[3], pat_hap, mat_hap, res_gen]

logging.info("Checking basic consistency...")
orig_ped_ind= frozenset(orig_ped.keys())
res_ped_ind= frozenset(res_ped.keys())
genotypes_ind= frozenset(genotypes.keys())
if ( res_ped_ind != genotypes_ind or
     orig_ped_ind > genotypes_ind ):
    logging.fatal("The two pedigrees refer to different individuals.")
    sys.exit(1)

logging.info("Trimming genotypes and haplotypes...")
for ind_id in orig_ped:
    ind_orig= orig_ped[ind_id]
    ind_res= res_ped[ind_id]
    gl_orig= len(ind_orig[5])
    gl_res= len(ind_res[5])
    if gl_orig != gl_res:
        gl= min(gl_orig, gl_res)
        ind_orig[3]= ind_orig[3][0:gl]
        ind_orig[4]= ind_orig[4][0:gl]
        ind_orig[5]= ind_orig[5][0:gl]
        ind_res[3]= ind_res[3][0:gl]
        ind_res[4]= ind_res[4][0:gl]
        ind_res[5]= ind_res[5][0:gl]
        genotypes[ind_id]= genotypes[ind_id][0:gl]

logging.info("Computing recombinations...")
# Compute recombinations in the original pedigree
orig_rec= []
for ind_id in orig_ped:
    ind= orig_ped[ind_id]
    if ind[1] != '0' and (ind[1] in orig_ped):
        orig_rec.extend(compute_recombinations(ind_id, ind[1],
                                               ind[3], orig_ped[ind[1]][3], orig_ped[ind[1]][4]))
    if ind[2] != '0' and (ind[2] in orig_ped):
        orig_rec.extend(compute_recombinations(ind_id, ind[2],
                                               ind[4], orig_ped[ind[2]][3], orig_ped[ind[2]][4]))
orig_rec=frozenset(orig_rec)
# Compute recombinations in the resulting haplotype configuration
res_rec= []
for ind_id in orig_ped:
    ind= res_ped[ind_id]
    if ind[1] != '0' and ind[1] in orig_ped:
        res_rec.extend(compute_recombinations(ind_id, ind[1],
                                              ind[3], res_ped[ind[1]][3], res_ped[ind[1]][4]))
    if ind[2] != '0' and ind[2] in orig_ped:
        res_rec.extend(compute_recombinations(ind_id, ind[2],
                                              ind[4], res_ped[ind[2]][3], res_ped[ind[2]][4]))
res_rec=frozenset(res_rec)
logging.info("Computing errors...")
orig_err= []
for ind_id in orig_ped:
    ind= orig_ped[ind_id]
    orig_err.extend([ (ind_id, l)
                      for l,g1,g2 in zip(range(len(genotypes[ind_id])), genotypes[ind_id], ind[5])
                      if g1 != 0 and g1 != g2 ])
orig_err=frozenset(orig_err)
res_err= []
for ind_id in res_ped:
    ind= res_ped[ind_id]
    res_err.extend([ (ind_id, l)
                     for l,g1,g2 in zip(range(len(genotypes[ind_id])), genotypes[ind_id], ind[5])
                      if g1 != 0 and g1 != g2 ])
res_err=frozenset(res_err)


logging.info("Checking pedigree structure consistency...")
for ind_id in orig_ped:
    ind_orig= orig_ped[ind_id]
    ind_res= res_ped[ind_id]
    if ind_orig[0:3] != ind_res[0:3]:
        logging.fatal("The parents of '%s' are different in the two pedigrees. "
                      "Original: (%s) - Result: (%s).",
                      ind_id,
                      ",".join(ind_orig[0:3]),
                      ",".join(ind_res[0:3]) )
        sys.exit(1)

logging.info("Computing pedigree stats...")
n_trios= Counter( ( tuple(orig_ped[ind_id][1:3])
                    for ind_id in orig_ped
                    if orig_ped[ind_id][1:3] != ('0', '0') ) )
n_children= Counter()
n_children.update( ( orig_ped[ind_id][1]
                     for ind_id in orig_ped
                     if orig_ped[ind_id][1] != '0' ) )
n_children.update( ( orig_ped[ind_id][2]
                     for ind_id in orig_ped
                     if orig_ped[ind_id][2] != '0' ) )


logging.info("Computing differences...")

tot_gen= 0
tot_het= 0
tot_hom= 0
tot_mis= 0

tot_abs_err_gen= 0
tot_abs_err_pat_hap= 0
tot_abs_err_mat_hap= 0
tot_abs_err_mask_pat_hap= 0
tot_abs_err_mask_mat_hap= 0
tot_err_gen= 0.0
tot_err_pat_hap= 0.0
tot_err_mat_hap= 0.0
tot_err_mask_pat_hap= 0.0
tot_err_mask_mat_hap= 0.0

loci= set()
err_gen_loci= Counter()
err_pat_hap_loci= Counter()
err_mat_hap_loci= Counter()
err_mask_pat_hap_loci= Counter()
err_mask_mat_hap_loci= Counter()
loci_counters= ( err_gen_loci, err_pat_hap_loci, err_mat_hap_loci,
                 err_mask_pat_hap_loci, err_mask_mat_hap_loci )


if options.full and options.header:
    print('"input_file"', '"result_file"',
          '"individual_id"',
          '"father_id"', '"mother_id"',
          '"genotype_length"',
          '"heterozygous"', '"homozygous"', '"missing"',
          '"genotype_errors"',
          '"paternal_haplotype_errors"', '"maternal_haplotype_errors"',
          '"paternal_haplotype_errors_wo_missing"', '"maternal_haplotype_errors_wo_missing"',
          sep="\t")
descs=( "# gen_err", "# pat_hap_err", "# mat_hap_err", "# pat_hap_err_miss", "# mat_hap_err_miss" )
norm_founders= not options.no_norm_found
norm_pseudofounders= options.norm_pseudo
for individual in orig_ped:

    genotype= genotypes[individual]
    ( err_gen,
      err_pat_hap, err_mat_hap,
      err_mask_pat_hap, err_mask_mat_hap,
      err_lists ) = compute_errors_with_swap(orig_ped,
                                             genotype,
                                             orig_ped[individual][3],
                                             orig_ped[individual][4],
                                             res_ped[individual][3],
                                             res_ped[individual][4],
                                             norm_founders, norm_pseudofounders)

    n_gen= len(genotype)
    n_het= sum( (g == 3) for g in genotype )
    n_hom= sum( (g == 1 or g == 2) for g in genotype )
    n_mis= sum( (g == 0) for g in genotype )
    for counter, l in zip(loci_counters, err_lists):
        counter.update(l)
        loci= loci.union(l)

    if options.full:
        (father, mother)= orig_ped[individual][1:3]
        is_founder= (father, mother) == ('0', '0')
        is_pseudo= ( father not in orig_ped ) or ( mother not in orig_ped )
        full_sibs= n_trios[(father, mother)] - 1
        half_sibs= n_children[father] + n_children[mother] - full_sibs - 2
        print(options.original,
              options.result,
              individual,
              orig_ped[individual][1],
              orig_ped[individual][2],
              len(genotypes[individual]),
              n_het, n_hom, n_mis,
              err_gen,
              err_pat_hap, err_mat_hap,
              err_mask_pat_hap, err_mask_mat_hap,
              n_children[individual],
              full_sibs, half_sibs,
              is_founder, is_pseudo,
              sep="\t")
        err_msg= [ "{}, {}, {}".format(d, individual, ", ".join([str(i) for i in l]))
                   for d,l in zip(descs, err_lists) if len(l)>0 ]
        if len(err_msg)>0:
            print("\n".join(err_msg))

    tot_gen += n_gen
    tot_het += n_het
    tot_hom += n_hom
    tot_mis += n_mis

    tot_abs_err_gen += err_gen
    tot_abs_err_pat_hap += err_pat_hap
    tot_abs_err_mat_hap += err_mat_hap
    tot_abs_err_mask_pat_hap += err_mask_pat_hap
    tot_abs_err_mask_mat_hap += err_mask_mat_hap
    tot_err_gen += (err_gen / n_mis) if n_mis > 0 else 0
    tot_err_pat_hap += (err_pat_hap / (n_het + n_mis)) if (n_het + n_mis) > 0 else 0
    tot_err_mat_hap += (err_pat_hap / (n_het + n_mis)) if (n_het + n_mis) > 0 else 0
    tot_err_mask_pat_hap += (err_mask_pat_hap / n_het) if n_het > 0 else 0
    tot_err_mask_mat_hap += (err_mask_pat_hap / n_het) if n_het > 0 else 0

n_indiv= len(orig_ped.keys())

# tot_gen /= n_indiv
# tot_het /= n_indiv
# tot_hom /= n_indiv
# tot_mis /= n_indiv
tot_err_gen /= n_indiv
tot_err_pat_hap /= n_indiv
tot_err_mat_hap /= n_indiv
tot_err_mask_pat_hap /= n_indiv
tot_err_mask_mat_hap /= n_indiv

if not options.full:
    if options.header:
        print('"input file"',
              '"result file"',
              '"pedigree size"',
              '"genotype length"',
              '"tot heterozygous loci"',
              '"tot homozygous loci"',
              '"tot missing genotypes"',
              '"tot genotype errors"',
              '"tot paternal haplotype errors"',
              '"tot maternal haplotype errors"',
              '"tot paternal haplotype errors wo missing"',
              '"tot maternal haplotype errors wo missing"',
              '"avg genotype errors"',
              '"avg paternal haplotype errors"',
              '"avg maternal haplotype errors"',
              '"avg paternal haplotype errors wo missing"',
              '"avg maternal haplotype errors wo missing"',
              '"orig recomb"',
              '"computed recomb"',
              '"TP recomb"',
              '"FN recomb"',
              '"FP recomb"',
              '"orig errors"',
              '"computed errors"',
              '"TP errors"',
              '"FN errors"',
              '"FP errors"',
              sep="\t")

    print(options.original,
          options.result,
          n_indiv,
          tot_gen/n_indiv,
          tot_het,
          tot_hom,
          tot_mis,
          tot_abs_err_gen,
          tot_abs_err_pat_hap,
          tot_abs_err_mat_hap,
          tot_abs_err_mask_pat_hap,
          tot_abs_err_mask_mat_hap,
          tot_err_gen,
          tot_err_pat_hap,
          tot_err_mat_hap,
          tot_err_mask_pat_hap,
          tot_err_mask_mat_hap,
          len(orig_rec),
          len(res_rec),
          len(orig_rec & res_rec),
          len(orig_rec - res_rec),
          len(res_rec - orig_rec),
          len(orig_err),
          len(res_err),
          len(orig_err & res_err),
          len(orig_err - res_err),
          len(res_err - orig_err),
          sep="\t")
else:
    if options.header:
        print('"input_file"', '"result_file"',
              '"locus"',
              '"genotype_errors"',
              '"paternal_haplotype_errors"', '"maternal_haplotype_errors"',
              '"paternal_haplotype_errors_wo_missing"', '"maternal_haplotype_errors_wo_missing"',
              sep="\t")
    for locus in range(max(loci)+1):
        print(options.original,
              options.result,
              locus,
              "\t".join([ str(counter[locus]) for counter in loci_counters ]),
              sep="\t")
