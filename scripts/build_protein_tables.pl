#!/usr/bin/perl

# build_protein_tables.pl - This script takes SIPROS's parsed output and
# parses it into three files: a protein table, a peptide table, and a table
# combining both protein and peptide information.  While performing this
# parse, it also clusters peptides by enrichment % and classifies different
# clusters as separate identifications.

# usage: build_protein_tables.pl <output prefix> <threshold>   <   <input_file>
# output prefix = will write the three filenames beginning with this prefix
# threshold = clustering threshold, i.e. clusters greater than this 
#             distance will not be merged
# This script reads from standard input.

$| = 1;
if(@ARGV != 2) { die "usage: $0 <prefix> <threshold, i.e. 10.0>\n"; }
$prefix = shift;
$thresh = shift;

# Open the three output files and print header information
open FH, ">$prefix.protein_table.txt" or die "couldn't open protein table file for writing\n";
open GH, ">$prefix.expanded_protein_table.txt" or die "couldn't open expanded protein table file for writing\n";
open WH, ">$prefix.peptide_table.txt" or die "couldn't open peptide table file for writing\n";
print FH "ProteinName\tAverageEnrichmentRatio\tStandardDeviation\tPepCounts\tScanCounts\tUniqueScanCounts\tScanEnrichRatios\tDescription\n";
print WH "RawFile\tScanNumber\tSequence\tCharge\tEnrichmentRatio\tScore\tProteins\n";
print GH "pro\tProteinName\tAverageEnrichmentRatio\tStandardDeviation\tPepCounts\tScanCounts\tUniqueScanCounts\tScanEnrichRatios\tDescription\n";
print GH "pep\tRawFile\tScanNumber\tSequence\tCharge\tEnrichmentRatio\tScore\tProteins\n";

$state = 0;

# Parse the input file and build hashes from protein (prot) and peptide (pep) perspectives.
while($line = <STDIN>) {
  if($line =~ /Organized by Scan/) { $state = 1; }
  if($line =~ /Organized by Pept/) { $state = 2; }
  if($line =~ /Organized by Prot/) { $state = 3; }
  if($line =~ /Files Processed/) { $state = 4; }
  if($line =~ /Database Entries/) { $state = 5; }
  if($state == 1) {
    next if($line !~ /^\d/);
    @inf = split /[ \n\r\t]+/, $line;
    $pep = $inf[2];
    $prot = $inf[8];
    $pct = $inf[9];
    $pct =~ s/\%//g;
    $pinfo{$prot} .= "$pct#";
    $prot2scan{$prot} .= $line;
    $pep2scan{$pep} .= $line;
  }
  if($state == 2) {
    next if($line !~ /EXACT/);
    @inf = split /[ \n\r\t]+/, $line;
    $pep = $inf[0];
    $prot = $inf[8];
    $pep2prot{$pep} .= "$inf[8],"; 
  }
  if($state == 3) {
    next if($line !~ /EXACT/);
    @inf = split /[ \n\r\t]+/, $line;
    $prot = $inf[0];
    $peps{$prot} .= $line;
  }
  if($state == 4) {
    next if($line !~ /^\d/);
    @inf = split /[ \n\r\t]+/, $line;
    $raw[$inf[0]] = $inf[1];
  }
  if($state == 5) {
    next if($line !~ /^\S/);
    next if($line =~ /^Database Ent/);
    @inf = split /[\n\r\t]+/, $line;
    $prot_desc{$inf[0]} = $inf[1]; 
  }
}

foreach(keys %pep2prot) { chop $pep2prot{$_}; }

# Perform agglomerative hierarchical clustering on the enrichment %'s
# for each peptide.
@skeys = sort { $a cmp $b } (keys %pinfo);
foreach(@skeys) {
  @inf = split /#/, $pinfo{$_};
  @sinf = sort { $a <=> $b } @inf;
  @avg = (); @nele = ();

  # Initialize clusters.  Each cluster begins with 1 element.
  # nele = count of elements, clus = list of elements, avg = average
  # enrichment % of the cluster
  for($i = 0; $i < @sinf; $i++) {
    $avg[$i] = $sinf[$i];
    $nele[$i] = 1;
    $clus[$i] = $sinf[$i];
  }
  for($i = 0; $i < @sinf; $i++) {
# print "MERGE ITER $i\t$stub\t$_\t";
# for($f = 0; $f < @avg; $f++) { print "$avg[$f] $nele[$f]#"; } print "\n";
    $mindist = $thresh; $minndx = -1;
    for($j = 0; $j < @avg-1; $j++) {
      $dist = $avg[$j+1] - $avg[$j];
#      $dist = $avg[$j]*$nele[$j] + $avg[$j+1]*$nele[$j+1];
#      $dist /= ($nele[$j] + $nele[$j+1]);
      # We found a candidate for merging (minimum distance < threshold).
      if($dist < $mindist) {
        $mindist = $dist;
        $minndx = $j;
      }
    }
    last if($minndx == -1);
    
    # Merge the two clusters.
    $avg[$minndx] = $avg[$minndx]*$nele[$minndx] + $avg[$minndx+1]*$nele[$minndx+1];
    $avg[$minndx] /= ($nele[$minndx] + $nele[$minndx+1]);
    $nele[$minndx] = $nele[$minndx] + $nele[$minndx+1];
    $clus[$minndx] .= ",".$clus[$minndx+1];

    # Delete the old cluster.
    splice(@avg, $minndx+1, 1);
    splice(@nele, $minndx+1, 1);
    splice(@clus, $minndx+1, 1);
  }

  # Print out the protein/peptide information for each cluster.
  # The subroutines peptide_info and scan_count_info return this information.
  # sdev = standard deviation, peptide_counter = counter for # of peptides for a
  # protein.
  for($i = 0; $i < @avg; $i++) {
    $avg[$i] = sprintf "%.2f", $avg[$i];
    @vals = split /\,/, $clus[$i];
    $low_enrich = $vals[0];
    $high_enrich = $vals[@vals-1];
    $sdev[$i] = sprintf "%.3f", sdev(@vals);

    $peptide_counter = 0;
    %peptide_tracker = {};
    $pi = peptide_info($_, $low_enrich, $high_enrich);

    $prot_line = "$_\t$avg[$i]\t$sdev[$i]\t$peptide_counter\t";
    $prot_line .= scan_count_info($_, $low_enrich, $high_enrich);
    $prot_line .= "$clus[$i]\t$prot_desc{$_}\n";
    print FH $prot_line;
    print GH "pro\t$prot_line";
    print WH $pi;
    $pi =~ s/\n/\npep\t/g;
    $pi =~ s/pep\t$//g;
    print GH "pep\t$pi";
  }
}
close FH;
close WH;
close GH;

# Subroutines to print out information for a scan and a peptide.

sub scan_count_info() {
  my $prot = shift @_;
  my $low = shift @_;
  my $high = shift @_;
  my @lines;
  @lines = split /\n/, $prot2scan{$prot};
  my $i, my @inf;
  my $scan_count = 0;
  my $scan_from_uni = 0;
  my $uni_check = 0;
  for($i = 0; $i < @lines; $i++) {
    @inf = split /[ \n\r\t]+/, $lines[$i];
    $pct = $inf[9];
    $pct =~ s/\%//g;
    next if($pct < $low || $pct > $high);
    $scan_count++;
    $uni_check = ($pep2scan{$inf[2]} =~ tr/\n/\n/);
    if($uni_check == 1) { $scan_from_uni++; }
  }
  return "$scan_count\t$scan_from_uni\t";
}

sub peptide_info() {
  my $prot = shift @_;
  my $low = shift @_;
  my $high = shift @_;
  my @lines;
  @lines = split /\n/, $prot2scan{$prot};
  my $i, my @inf;
  my $rstring = "";
  for($i = 0; $i < @lines; $i++) {
    @inf = split /[ \n\r\t]+/, $lines[$i];
    $pct = $inf[9];
    $pct =~ s/\%//g;
    next if($pct < $low || $pct > $high);

    # count peptides
    $ptrack = "$inf[2]#$inf[7]";
    if($peptide_tracker{$ptrack} eq "") {
      $peptide_tracker{$ptrack} = 1;
      $peptide_counter++;
    }

    $rstring .= "$raw[$inf[0]]\t";
    $rstring .= "$inf[1]\t$inf[2]\t$inf[7]\t$inf[9]\t$inf[5]\t$pep2prot{$inf[2]}\n";
  }
  return $rstring;
}

sub sdev() {
  my @vals = @_;
  my $i, my $avg = 0.0, my $sdev = 0.0;
  my $nv = @vals;
  if($nv == 0) { return 0; }
  foreach(@vals) { $avg+=$_; }
  $avg /= $nv;
  foreach(@vals) { $sdev += ($_-$avg)*($_-$avg); }
  $sdev /= $nv;
  $sdev = sqrt($sdev);
  return $sdev;
}
