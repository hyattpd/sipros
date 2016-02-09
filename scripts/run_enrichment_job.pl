#!/usr/bin/perl

# run_enrichment_job.pl

# This script writes a qsub script for the PBS queueing system and executes it.
# Expects there to be 101 processors available on which to run a parallel job.

# The output will be 101 .psc files in the working directory, outout.1.psc through
# output.101.psc.

# usage: run_enrichment_job.pl <sipros> <config directory> <working directory>

# sipros = the full path to the sipros binary
# config directory = the directory in which you placed the config files
# working directory = a working directory in which you have placed your FT2 files
#                     for analysis

$sipros = shift;
$config = shift;
$working_dir = shift;
if($sipros eq "") { die "usage: $0 <sipros> <config directory> <working directory>\n"; }
if($working_dir eq "") { die "usage: $0 <sipros> <config directory> <working directory>\n"; }
if($config eq "") { die "usage: $0 <sipros> <config directory> <working directory>\n"; }

$ctr = 1;

open FH, ">$working_dir/job.txt" or die "couldn't create qsub script\n";
print FH "#PBS -t 1-101\n";
print FH "$sipros -w $working_dir -c $config/config." . '$PBS_ARRAYID' . " > $working_dir/output." . '$PBS_ARRAYID' . ".psc\n";
close FH;

system "qsub $working_dir/job.txt";
