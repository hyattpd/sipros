#!/usr/bin/perl

# parse_raw_output.pl - This script parses the 101 raw output files
# from the working directory.  These files are of the form output.1.psc
# to output.101.psc.

# usage: parse_raw_output.pl -w <working_directory> -t <score thresholdold>
# working_directory = directory containing the 101 output files
# score thresholdold = minimum score to retain an entry, 25.0 is the default

$| = 1;

$working_dir = ".";
$threshold = 25.0;

# Parse command line arguments
for($i = 0; $i < @ARGV; $i++) {
  if($ARGV[$i] eq "-w" || $ARGV[$i] eq "-W") {
    if($i == @ARGV-1) { die "$0: -w option requires argument.\n"; }
    $working_dir = $ARGV[++$i];
  }
  elsif($ARGV[$i] eq "-t" || $ARGV[$i] eq "-T") {
    if($i == @ARGV-1) { die "$0: -t option requires argument.\n"; }
    $threshold = $ARGV[++$i];
  }
  else { die "$0: Unrecognized option $ARGV[$i]\n"; }
}

# Get .psc file list for processing
opendir(DH, $working_dir) or die "$0: couldn't open $working_dir for reading\n";
while($line = readdir(DH)) {
  if($line =~ /\.[pP][sS][cC]$/) {
    push @psc_files, "$working_dir/$line";
  }
  if($line =~ /\.[fF][tT]2$/) {
    push @ft2_files, "$working_dir/$line";
  }
}
closedir(DH);

for($i = 0; $i < @ft2_files; $i++) {
  $ft2_files[$i] = ($ft2_files[$i] =~ /\/([^\/]+)$/)[0];
}
@sft2_list = sort { $a cmp $b } @ft2_files;

# Print file names with indices
print "Files Processed by Polyscan\n\n";
print "Index\tFileName\n";
for($i = 0; $i < @sft2_list; $i++) {
  $ndx = $i+1;
  print "$ndx\t$sft2_list[$i]\n";
  $ft2_index{$sft2_list[$i]} = $ndx;
}

# Go through each psc file
$global_min = 1000.0; # used to deduce polyscan's thresholdold
for($i = 0; $i < @psc_files; $i++) {
  $psc_num = ($psc_files[$i] =~ /\.(\d+)\./)[0];
  print STDERR "Processing $psc_files[$i]...\n";
  open FH, $psc_files[$i] or die "$0: couldn't open .psc file $psc_files[$i] for reading\n";
  while($line = <FH>) {
    @inf = split /[\n\t]+/, $line;
    next if(@inf != 8);
    next if($inf[1] =~ /[^0-9]/);
    next if($inf[7] =~ /[^0-9]/);
    next if($inf[2] =~ /[^\.0-9]/);
    next if($inf[5] =~ /[^A-Z]/); 

    $inf[0] = ($inf[0] =~ /\/([^\/]+)$/)[0];
    if($inf[2] < $global_min) { $global_min = $inf[2]; }
    $ft2 = $ft2_index{$inf[0]};
    if($ft2 eq "") { die "$0: could not find ft2 file $inf[0] in the index\n"; }

    $scankey = "$ft2#$inf[1]";
    $inf[4] =~ s/[\[\]]//g;
    if($inf[4] eq $inf[5]) { $type = 0; } # Exact match
    elsif($inf[4] =~ /[^A-Z]/) { $type = 2; } # PTM
    else { 
      $type = 1;  # Polymorphism
      for($j = 0; $j < length($inf[4]); $j++) {
        if(substr($inf[4], $j, 1) ne substr($inf[5], $j, 1) &&
           substr($inf[4], $j, 1) eq "L") {
          substr($inf[4], $j, 1) = "I"; 
        }
      }
    }
    @dbinf = split /[ \n\t]+/, $inf[6];
    $locus = $dbinf[0];
    if($saw_db_header{$locus} eq "") { $saw_db_header{$locus} = $inf[6]; }

    # Case 1:  This result beats the current best score for this scan
    if(!defined($scan_best_score{$scankey}) || $inf[2] > $scan_best_score{$scankey}) {

      $scan_second_score{$scankey} = $scan_best_score{$scankey};

      $scan_best_score{$scankey} = $inf[2];
      $scan_best_charge{$scankey} = $inf[3];
      $scan_best_peptide{$scankey} = $inf[4];
      $scan_best_db_peptide{$scankey} = $inf[5];
      $scan_best_type{$scankey} = $type;
      $scan_best_locus{$scankey} = $locus;
      $scan_best_db_pos{$scankey} = $inf[7];
      $scan_best_psc{$scankey} = $psc_num;     
    }

    # Case 2:  This result ties the current best score for this scan

    elsif($inf[2] == $scan_best_score{$scankey}) {
      $scan_best_charge{$scankey} .= "#$inf[3]";
      $scan_best_peptide{$scankey} .= "#$inf[4]";
      $scan_best_db_peptide{$scankey} .= "#$inf[5]";
      $scan_best_type{$scankey} .= "#$type";
      $scan_best_locus{$scankey} .= "#$locus";
      $scan_best_db_pos{$scankey} .= "#$inf[7]";
      $scan_best_psc{$scankey} .= "#$psc_num";     
    }

    # Case 3:  This result beats the current second best score for this scan

    elsif(!defined($scan_second_score{$scankey}) || $inf[2] > $scan_second_score{$scankey}) {
      $scan_second_score{$scankey} = $inf[2];
    }

  }
  close FH;
}

$global_min = int($global_min)*1.0;
@sscans = sort { @sinfa = split /#/, $a; @sinfb = split /#/, $b; $sinfa[0] <=> $sinfb[0] or $sinfa[1] <=> $sinfb[1]; } (keys %scan_best_score);
@sheaders = sort { $a cmp $b } (keys %saw_db_header);

# Print out the database information so it won't have to be looked up later.
# Note this generates an enormous amount of output (every located protein header).
print "\n\n\nDatabase Entries\n\n";
foreach(@sheaders) { print "$_\t$saw_db_header{$_}\n"; }

# Print out the parsed output information for each scan number.
print "\n\n\nOutput Organized by Scan\n\n";
print "FileIndex\tScanNo\tPeptide\tType\tDBPeptide\tScore\tDCN\tCharge\tLocus\tEnrichPct\tBest\n";
for($i = 0; $i < @sscans; $i++) {
  $score = $scan_best_score{$sscans[$i]};
  next if($score < $threshold);
  if(!defined($scan_second_score{$sscans[$i]})) {
    # dcn = Delta CN
    $dcn = $score - $global_min;
    $dcn .= "*";
  }
  else {
    $dcn = $score - $scan_second_score{$sscans[$i]};
  }
  ($find, $snum) = split /#/, $sscans[$i];
  
  @tinf1 = split /#/, $scan_best_peptide{$sscans[$i]};
  @tinf2 = split /#/, $scan_best_type{$sscans[$i]};
  @tinf3 = split /#/, $scan_best_db_peptide{$sscans[$i]};
  @tinf4 = split /#/, $scan_best_locus{$sscans[$i]};
  @tinf5 = split /#/, $scan_best_db_pos{$sscans[$i]};
  @tinf6 = split /#/, $scan_best_psc{$sscans[$i]};
  @tinf7 = split /#/, $scan_best_charge{$sscans[$i]};

  @tind = (); @stind = ();
  for($j = 0; $j < @tinf2; $j++) { 
    $tind[$j] = $j;
  }
  @stind = sort {
    if($tinf4[$a] =~ /^Reverse/) { $r1 = 1; }
    else { $r1 = 0; }
    if($tinf4[$b] =~ /^Reverse/) { $r2 = 1; }
    else { $r2 = 0; }
    $tinf2[$a] <=> $tinf2[$b] or $r1 <=> $r2 or
    $tinf1[$a] cmp $tinf1[$b] or $tinf6[$a] <=> $tinf6[$b] or
    $tinf4[$a] cmp $tinf4[$b];
  } @tind;

  # Order of preference for best hit = EXACT match, then polymorphism, then PTM
  # (Irrelevant for this enrichment version.)
  $tstr[0] = "EXACT"; $tstr[1] = "POLY"; $tstr[2] = "PTM";
  for($j = 0; $j < @tind; $j++) {
    $n = $stind[$j];
    if($j > 0) { 
      last if($tinf2[$n] != $tinf2[0]);
      $o = $stind[$j-1]; 
      next if($tinf1[$n] eq $tinf1[$o] && $tinf4[$n] eq $tinf4[$o] && $tinf3[$n] eq $tinf3[$o]); 
    }
  
    # Protein Perspective - Compile Information
    $protein_peptide{$tinf4[$n]} .= "$tinf1[$n]#";
    $protein_type{$tinf4[$n]} .= "$tinf2[$n]#";
    $protein_dbpep{$tinf4[$n]} .= "$tinf3[$n]#";
    $protein_pos{$tinf4[$n]} .= "$tinf5[$n]#";

    $peptide_scan{$tinf1[$n]} .= "$sscans[$i]|"; 
    $peptide_locus{$tinf1[$n]} .= "$tinf4[$n]#";
    $peptide_type{$tinf1[$n]} .= "$tinf2[$n]#";
    $peptide_charge{$tinf1[$n]} .= "$tinf7[$n]#";

    # Peptide => Locus Information
    if($j == 0 || $tinf1[$n] ne $tinf1[$o]) { 
      $peptide_scan_count{$tinf1[$n]}++; 
    } 
    if($saw_pep_locus{"$tinf1[$n]#$tinf4[$n]"} != 1) {
      $peptide_locus_count{$tinf1[$n]}++;
      $saw_pep_locus{"$tinf1[$n]#$tinf4[$n]"} = 1;
    }
    print "$find\t$snum\t";
    printf "$tinf1[$n]\t$tstr[$tinf2[$n]]\t$tinf3[$n]\t%.2f\t%.2f\t$tinf7[$n]\t$ctr0\t$tinf4[$n]\t", $score, $dcn;
    $pct = $tinf6[$n] - 1;
    $pct = "$pct\%";
    print "$pct\t";
    if($j == 0) { print "*\n"; } else { print ".\n"; }
  }
}

# Output from a Protein Perspective
print "\n\n\nOutput Organized by Protein\n\n";
print "Protein\tBegin\tEnd\tPeptide\tType\tUni\tScanCnt\tLocusCnt\n";
@sloci = sort (keys %protein_peptide);

for($i = 0; $i < @sloci; $i++) {

  @pinf1 = split /#/, $protein_peptide{$sloci[$i]};
  @pinf2 = split /#/, $protein_type{$sloci[$i]};
  @pinf3 = split /#/, $protein_pos{$sloci[$i]};
  @pinf4 = split /#/, $protein_dbpep{$sloci[$i]};

  @pind = (); @spind = ();
  for($j = 0; $j < @pinf1; $j++) { 
    $pind[$j] = $j;
  }
  @spind = sort {
    $pinf3[$a] <=> $pinf3[$b] or
    length($pinf4[$a]) <=> length($pinf4[$b]) or
    $pinf1[$a] cmp $pinf1[$b] or
    $pinf2[$a] <=> $pinf2[$b];
  } @pind;

  $tstr[0] = "EXACT"; $tstr[1] = "POLY"; $tstr[2] = "PTM";
  for($j = 0; $j < @pind; $j++) {
    $n = $spind[$j];
    if($j > 0) { 
      $o = $spind[$j-1]; 
      next if($pinf1[$n] eq $pinf1[$o] && $pinf3[$n] eq $pinf3[$o]);
    }
    $end = length($pinf4[$n])+$pinf3[$n]-1;
    print "$sloci[$i]\t$pinf3[$n]\t$end\t$pinf1[$n]\t$tstr[$pinf2[$n]]\t";
    if($peptide_locus_count{$pinf1[$n]} == 1) { print "U\t"; }
    else { print "N\t"; }
    print "$peptide_locus_count{$pinf1[$n]}\t$peptide_scan_count{$pinf1[$n]}\n";
  }

}

# Print output from an individual peptide perspective.
print "\n\n\nOutput Organized by Peptide\n\n";
print "Peptide\tCharge\tUnique\tLocusCnt\tScanCnt\tType\tFileNo\tScanNo\tLocus\n";
@spep = sort (keys %peptide_type);

for($i = 0; $i < @spep; $i++) {

  @pinf1 = split /#/, $peptide_type{$spep[$i]};
  @pinf2 = split /\|/, $peptide_scan{$spep[$i]};
  @pinf3 = split /#/, $peptide_locus{$spep[$i]};
  @pinf4 = split /#/, $peptide_charge{$spep[$i]};

  @pind = (); @spind = ();
  for($j = 0; $j < @pinf1; $j++) { 
    $pind[$j] = $j;
  }
  @spind = sort {
    ($spfa, $spsa) = split /#/, $pinf2[$a];
    ($spfb, $spsb) = split /#/, $pinf2[$b];
    $pinf1[$a] <=> $pinf1[$b] or
    $spfa <=> $spfb or
    $spsa <=> $spsb or
    $pinf3[$a] cmp $pinf3[$b];
  } @pind;

  $tstr[0] = "EXACT"; $tstr[1] = "POLY"; $tstr[2] = "PTM";
  for($j = 0; $j < @pind; $j++) {
    $n = $spind[$j];
    if($j > 0) { 
      $o = $spind[$j-1]; 
      next if($pinf1[$n] eq $pinf1[$o] && $pinf2[$n] eq $pinf2[$o] && $pinf3[$n] eq $pinf3[$o]);
    }
    print "$spep[$i]\t$pinf4[$n]\t";
    if($peptide_locus_count{$pinf1[$n]} == 1) { print "U\t"; }
    else { print "N\t"; }
    print "$peptide_scan_count{$spep[$i]}\t$peptide_locus_count{$spep[$i]}\t";
    ($spf, $sps) = split /#/, $pinf2[$n];
    print "$tstr[$pinf1[$n]]\t$spf\t$sps\t$pinf3[$n]\n";
  }
}
