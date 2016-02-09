#!/usr/bin/perl

# clean_dir.pl

# usage: clean_dir.pl <directory to clean>

# This script removes all files from a directory that do not end in ".FT2".

$working_dir = shift;
opendir DH, $working_dir or die "couldn't open directory $working_dir\n";
while($line = readdir(DH)) {
  next if($line =~ /^\./);
  next if($line =~ /\.[Ff][Tt]2$/);
  unlink "$working_dir/$line";
}
closedir DH;
