#!/usr/bin/perl

# filter_peptides.pl - This script filters the three files, deleting proteins
# whose identifications are subsets of other protiens, and also deleting proteins
# with fewer than some specified number of peptide ids.

# usage:  filter_peptides.pl <prefix> <min # peptides> <min # scans> <do_subset> <do_grouping>
#
# prefix = The prefix for the three output files.
# min # peptides = Minimum number of peptides to retain an entry.
# min # scans = Minimum number of scans to retain an entry
# do_subset = 0 or 1, Should we delete proteins whose ids are subset of another?
# do_grouping = 0 or 1, Should we merge proteins whose ids are identical?

$| = 1;
$usage = "$0 <file stub> <min # peptides> <min # scans> [0|1(subset)] [0|1(group)]";

if(@ARGV < 3 || @ARGV > 5) { die $usage; }

$stub = shift;
$minpep = shift;
$minscan = shift;
$subset = shift;
if($subset eq "") { $subset = 1; }
$group = shift;
if($group eq "") { $group = 1; }

# Store the whole expanded protein input in @data
$expfile = "$stub.expanded_protein_table.txt";
open FH, $expfile or die "could not open $expfile for reading";
@data = <FH>;
close FH;

# Count proteins and peptides.  Do simple filtering of
# proteins that have less than min # of peptides or scans.
# cnt = Count of total proteins, npro = Revised count.
$state = 0;
$npro = 0; $npep = 0;
for($i = 2; $i < @data; $i++) {
  $line = $data[$i];
  if($line =~ /^pro/) { 
    $cnt++;
    @inf = split /[\n\t\r]+/, $line;
    if($inf[4] >= $minpep && $inf[5] >= $minscan) { 
      push @proteins, $line;
      $npro++;
      $state = 1;
    }
    else { $state = 0; }
  }
  if($line =~ /^pep/ && $state == 1) {
    @inf = split /[\n\t\r]+/, $line;
    $pep = "$inf[1]#$inf[2]#$inf[3]#$inf[4]";
    if($saw_pep{$pep} eq "") {
      $saw_pep{$pep} = $npep;
      $pepset[@proteins-1] .= "$npep#";
      push @peptides, $line;
      $npep++;
    }
    else {
      $pepset[@proteins-1] .= "$saw_pep{$pep}#";
    }
  }
}
print "$cnt proteins in expanded protein table\n";
print "$npro proteins after peptide count / scan filtering\n";

# Sort the peptide ids numerically (easier to compare)
for($i = 0; $i < @proteins; $i++) { 
  $line = $pepset[$i];
  @inf = split /#/, $line;
  @sinf = sort { $a <=> $b } @inf;
  $sline = join '#', @sinf;
}

# Merge groups with identical peptide ids
print "$npro proteins before grouping\n";
if($group == 1) {
  for($i = 0; $i < @proteins; $i++) {
    for($j = $i+1; $j < @proteins; $j++) {
      if($pepset[$i] eq $pepset[$j]) {
        @inf1 = split /[\r\t\n]+/, $proteins[$i];
        @inf2 = split /[\r\t\n]+/, $proteins[$j];
        $newline = "$inf1[0]\t$inf1[1],$inf2[1]\t$inf1[2]\t$inf1[3]\t$inf1[4]\t$inf1[5]\t$inf1[6]\t$inf1[7]\t$inf1[8],$inf2[8]\n";
        $proteins[$i] = $newline;
        splice @proteins, $j, 1;
        splice @pepset, $j, 1;
        $j--;
      }
    }
  }
  $npro = @proteins;
  print "$npro proteins after checking identity\n";
}

# Now look for subset using the sinf sorted lists of peptides.
# No shortcut here; have to check each peptide id.
# We make two passes, one forward, one reverse.
if($subset  == 1) {

  # Forward pass
  for($i = 0; $i < @proteins-1; $i++) {
    for($j = $i+1; $j < @proteins; $j++) {
      @inf1 = split /#/, $pepset[$i];
      @inf2 = split /#/, $pepset[$j];
      next if(@inf2 >= @inf1);
      $destroy = 1;
      for($k = 0; $k < @inf2; $k++) {
        $found = 0;
        for($l = 0; $l < @inf1; $l++) {
          if($inf1[$l] eq $inf2[$k]) { $found = 1; }
        } 
        if($found == 0) { $destroy = 0; last; }
      }
      if($destroy == 1) {
#        print "SUBSET\t$pepset[$i]\t$pepset[$j]\n";
        splice @proteins, $j, 1;
        splice @pepset, $j, 1;
        $j--;
      }
    }
  } 
  $npro = @proteins;
  print "$npro proteins after forward subset pass\n";

  # Reverse pass
  for($i = @proteins-1; $i > 0; $i--) {
    for($j = $i-1; $j >= 0; $j--) {
      @inf1 = split /#/, $pepset[$i];
      @inf2 = split /#/, $pepset[$j];
      next if(@inf2 >= @inf1);
      $destroy = 1;
      for($k = 0; $k < @inf2; $k++) {
        $found = 0;
        for($l = 0; $l < @inf1; $l++) {
          if($inf1[$l] eq $inf2[$k]) { $found = 1; }
        } 
        if($found == 0) { $destroy = 0; last; }
      }
      if($destroy == 1) {
#        print "SUBSET\t$pepset[$i]\t$pepset[$j]\n";
        splice @proteins, $j, 1;
        splice @pepset, $j, 1;
        $j++;
      }
    }
  } 
  $npro = @proteins;
  print "$npro proteins after reverse subset pass\n";
}

# Print out the revised tables.
$eptab = "$stub.filtered.expanded_protein_table.txt";
$protab = "$stub.filtered.protein_table.txt";
$peptab = "$stub.filtered.peptide_table.txt";

open FH, ">$protab" or die "couldn't open $protab for writing...\n";
open GH, ">$eptab" or die "couldn't open $eptab for writing...\n";
open WH, ">$peptab" or die "couldn't open $peptab for writing...\n";
print FH "ProteinName\tAverageEnrichmentRatio\tStandardDeviation\tPepCounts\tScanCounts\tUniqueScanCounts\tScanEnrichRatios\tDescription\n";
print WH "RawFile\tScanNumber\tSequence\tCharge\tEnrichmentRatio\tScore\tProteins\n";
print GH "pro\tProteinName\tAverageEnrichmentRatio\tStandardDeviation\tPepCounts\tScanCounts\tUniqueScanCounts\tScanEnrichRatios\tDescription\n";
print GH "pep\tRawFile\tScanNumber\tSequence\tCharge\tEnrichmentRatio\tScore\tProteins\n";

for($i = 0; $i < @proteins; $i++) {
  print GH $proteins[$i];
  $proteins[$i] =~ s/^pro\t//;
  print FH $proteins[$i];
  @inf = split /#/, $pepset[$i]; 
  foreach(@inf) { print GH $peptides[$_]; }
}
for($i = 0; $i < @peptides; $i++) {
  $peptides[$i] =~ s/^pep\t//g;
  print WH $peptides[$i];
}
close FH;
close WH;
close GH;
