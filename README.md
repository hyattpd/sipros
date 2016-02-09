README for running SIPROS jobs

1.  Write the config files for either carbon or nitrogen to a directory.

Example:
write_config_files.pl default_sip_config.xml myfastadb.fna /home/me/config_files/ C

# Usage: write_config_files.pl <config file> <database path> <config directory> < C | N >
# config file = the config file to read in
# database path = the FASTA database you want to search, needed to write into the config file
# config directory = the directory you want to store the config files in
# C | N = carbon or nitrogen enrichment

2.  Run the enrichment job on a cluster with a PBS queue system.

Example:
run_enrichment_job.pl /home/me/sipros.exe /home/me/config_files /home/me/myrun/

# usage: run_enrichment_job.pl <polyscan> <config directory> <working directory>

# polyscan = the full path to the polyscan binary
# config directory = the directory in which you placed the config files
# working directory = a working directory in which you have placed your FT2 files
#                     for analysis

3.  Parse the raw output on the cluster.

Example:
parse_raw_output.pl -w /home/me/myrun -t 25.0 > myrawoutput.txt

# usage: parse_raw_output.pl -w <working_directory> -t <score thresholdold>
# working_directory = directory containing the 101 output files
# score thresholdold = minimum score to retain an entry, 25.0 is the default

4.  Copy the file to your local machine (or continue work on the cluster, up to you).

5.  Cluster by enrichment % and print out the three table files.

Example:
build_protein_tables.pl /home/me/zpb/finalout 10.0 < myrawoutput.txt

# usage: build_protein_tables.pl <output prefix> <threshold>   <   <input_file>
# output prefix = will write the three filenames beginning with this prefix
# threshold = clustering threshold, i.e. clusters greater than this
#             distance will not be merged
# This script reads from standard input.

6.  Filter the peptides, removing proteins that are subsets of other proteins
and proteins with fewer than some threshold # of peptides identified.

Example:
filter_peptides.pl /home/me/zpb/finalout 2 2 1 1

# usage:  filter_peptides.pl <prefix> <min # peptides> <min # scans> <do_subset> <do_grouping>
#
# prefix = The prefix for the three output files.
# min # peptides = Minimum number of peptides to retain an entry.
# min # scans = Minimum number of scans to retain an entry
# do_subset = 0 or 1, Should we delete proteins whose ids are subset of another?
# do_grouping = 0 or 1, Should we merge proteins whose ids are identical?

7.  Run ProRata for quantification information.

Example:
run_prorata.pl ProRataConfig_Template.xml DTASelect-filter_Template.txt /home/me/quant

# usage: run_prorata.pl <prorataconfig> <dtafiltertemplate> <working_dir>
# prorataconfig = ProRata config file
# dtafiltertemplate = DTA Select template to read and copy
# working_dir = The directory containing the three SIPROS protein table files.

