#!/usr/bin/perl

# write_config_files.pl

# This script takes the default config file and writes 101 config files to a specified
# directory, one for each % enrichment level of carbon or nitrogen.

# Usage: write_config_files.pl <config file> <database path> <config directory> < C | N >
# config file = the config file to read in
# database path = the FASTA database you want to search, needed to write into the config file
# config directory = the directory you want to store the config files in
# C | N = carbon or nitrogen enrichment

$config = shift;
$database = shift;
$working_dir = shift;
$enrich_target = shift;
if($config eq "") { die "usage: $0 <config file> <database path> <config directory> < C | N >\n"; }
if($database eq "") { die "usage: $0 <config file> <database path> <config directory> < C | N >\n"; }
if($working_dir eq "") { die "usage: $0 <config file> <database path> <config directory> < C | N >\n"; }
if($enrich_target eq "") { die "usage: $0 <config file> <database path> <config directory> < C | N >\n"; }

# 30 and 42 are the line numbers in the config file to modify for carbon and nitrogen %'s
if($enrich_target =~ /^[Cc]/) { $enrich_target = 30; }
elsif($enrich_target =~ /^[Nn]/) { $enrich_target = 42; }
else { die "enrich_target must be one of C or N\n"; }

$counter = 1;

open FH, $config or die "couldn't open config file\n";
@content = <FH>;
close FH;

# Loop through in 0.01 increments and set the enrichment levels in each config file
for($i = 0; $i <= 1.006; $i+= 0.01) {
  open FH, ">$working_dir/config.$counter" or die "couldn't open config.$counter for writing (check working dir: $working_dir)\n";
  for($j = 0; $j < @content; $j++) {
    if($j == 3) { print FH "        <FASTA_DATABASE>$database<\/FASTA_DATABASE>\n"; }
    elsif($j == $enrich_target && $counter != 2) { 
      if($content[$j] !~ /PERCENT/) { die "bad line in config $content[$j]\n"; }
      printf FH "\t\t<PERCENTAGE>\t%.4f,\t%.4f\t</PERCENTAGE>\n", 1.0-$i, $i;
    }
    else { print FH $content[$j]; }
  }
  close FH;
  $counter++;
}
