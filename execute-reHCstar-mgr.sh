#!/bin/bash

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
#  execute-reHCstar-mgr.sh
#
#  A simple script that execute reHC in the 'create-read' mode on
#  every genotyped pedigree of the current directory.
#
#  Useful for 'serial' experimentations.
#
##########

LANG=C

# Check existence of required files
if [ ! -f "./config-execution-mgr.ini"  ];
then
    echo "Configuration file 'config-execution-mgr.ini' not found. Aborting..."
    exit 1
fi

# Set parameters to default values
reHCmgr_exe=""
reHCmgr_blocks=""
reHCmgr_opts=""
reHC_exe=""
reHC_sat_mode=""
error_handler=""
recomb_handler=""
pedigrees=gen-ped-*.txt
dest_dir="./rehcstar-out/"
time_limit=0   # no limit
memory_limit=0   # no limit

# Read configuration
source ./config-execution-mgr.ini

# Check required parameters...
if [ -z "${reHCmgr_exe}" ]; then echo "Parameter 'reHCmgr_exe' not defined. Aborting..."; exit 1; fi
if [ -z "${reHC_exe}" ]; then echo "Parameter 'reHC_exe' not defined. Aborting..."; exit 1; fi
if [ -z "${reHC_sat_mode}" ]; then echo "Parameter 'reHC_sat_mode' not defined. Aborting..."; exit 1; fi
if [ -z "${pedigrees}" ]; then echo "Parameter 'pedigrees' not defined. Aborting..."; exit 1; fi
time_limit=`echo ${time_limit} | sed 's/[^0-9]//g'`   # Remove non-digits
memory_limit=`echo ${memory_limit} | sed 's/[^0-9]//g'`
if [ -z "${time_limit}" ]; then echo "Parameter 'time_limit' is not correctly defined. Aborting..."; exit 1; fi
if [ -z "${memory_limit}" ]; then echo "Parameter 'memory_limit' is not correctly defined. Aborting..."; exit 1; fi

if [ ! -x "${reHCmgr_exe}" ]; then echo "Manager program '${reHCmgr_exe}' not found or not executable. Aborting..."; exit 1; fi
if [ ! -x "${reHC_exe}" ]; then echo "reHCstar program '${reHC_exe}' not found or not executable. Aborting..."; exit 1; fi

for pedigree in ${pedigrees}; do
    if [ ! -f ${pedigree} ]; then
        echo "Pedigree '${pedigree}' refers to a non-existent file. Aborting..."
        exit 1
    fi
done

dest_dir="${dest_dir%\/}"
if [ ! -e "${dest_dir}" ]; then
    mkdir -p "${dest_dir}"
fi
if [ ! -d "${dest_dir}" ]; then
    echo "Impossible to create output directory '${dest_dir}'. Aborting..."
    exit 1
fi


reHCmgr_partial="${reHCmgr_exe} ${reHCmgr_opts} ${reHCmgr_blocks} --cmd \"${reHC_exe} ${reHC_sat_mode} ${error_handler} ${reHC_io}\" --cmd-rec \"${recomb_handler}\""

echo "`date`  --  Starting experimentation"
echo "`date`  --  Starting experimentation" > execution.log
limits=""
if [ ${time_limit} -gt 0 ]; then
    echo "`date`  --  Setting time limit to ${time_limit} minutes (each)"
    echo "`date`  --  Setting time limit to ${time_limit} minutes (each)" >> execution.log
    limits="${limits} ulimit -t $(( time_limit * 60 ));"
fi
if [ ${memory_limit} -gt 0 ]; then
    echo "`date`  --  Setting memory limit to ${memory_limit} MB (each)"
    echo "`date`  --  Setting memory limit to ${memory_limit} MB (each)" >> execution.log
    limits="${limits} ulimit -v $(( memory_limit * 1024 ));"
fi
for full_pedigree in ${pedigrees}; do
    pedigree=`basename ${full_pedigree}`
    echo "`date`  --  Executing on ${pedigree}"
    echo "`date`  --  Executing on ${pedigree}" >> execution.log
    (
        if [ -f "config-execution-mgr-${pedigree}.ini" ]; then
            source "./config-execution-mgr-${pedigree}.ini"
        fi
        echo "${limits# } nice time -f \"%U %S %E %x %M %C\" -o ${dest_dir}/time-${pedigree}" \
            "${reHCmgr_partial} ${reHCmgr_local_opts}" \
            "-p ${full_pedigree} -r ${dest_dir}/hap-${pedigree} > ${dest_dir}/log-${pedigree} 2>&1" > "${dest_dir}/cmd-${pedigree}"
        source "${dest_dir}/cmd-${pedigree}"
    )
done
echo "`date`  --  Experimentation ended"
echo "`date`  --  Experimentation ended" >> execution.log
