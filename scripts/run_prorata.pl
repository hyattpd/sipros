#!/usr/bin/perl

# run_prorata.pl - Script that runs ProRata using the expanded protein tables.

# usage: run_prorata.pl <prorataconfig> <dtafiltertemplate> <working_dir>
# prorataconfig = ProRata config file
# dtafiltertemplate = DTA Select template to read and copy
# working_dir = The directory containing the three SIPROS protein table files.

$| = 1;

if(@ARGV != 3) { die "usage: $0 <prorataconfig> <dtafiltertemplate> <working_dir>"; }
$prcfg = shift;
$dtaf = shift;
$wd = shift;

$done = 0;
opendir DH, $wd or die "couldn't open $wd";
while($line = readdir DH) {
  next if($line !~ /expanded_protein_table\.txt$/);
  if($done == 1) { die "multiple protein tables in $wd"; }
  $file = "$wd/$line";
  $done = 1;
}
closedir DH;
$config = "$wd/ProRataConfig.xml";
$dtafilter = "$wd/DTASelect-filter.txt";
$pro_out = "$wd/ProRata_Quantification_Protein.txt";
$final_output = "$wd/polyscan_prorata.output.txt";

open OH, ">$final_output" or die "couldn't open output file $final_output for writing\n";

# Open the two config files and read the data in.
open FH, $prcfg or die "couldn't open prorata config file $prcfg";
@prorata = <FH>;
close FH;
foreach(@prorata) { $_ =~ s/\r//g; }
open FH, $dtaf or die "couldn't open dta filter template file $dtaf";
@dtalines = <FH>;
close FH;
foreach(@dtalines) { $_ =~ s/\r//g; }

open FH, $file or die "couldn't open $file for reading";
while($line = <FH>) {
  if($line =~ /^pro/) {
    next if($line =~ /ProteinName/);
    @inf = split /[ \n\r\t]+/, $line;
    $enrich = sprintf "%.3f", $inf[2]/100.0;
    $comp_pro = "$inf[1]#$enrich";
    $pro2enr{$inf[1]} .= "$enrich#";
  }
  elsif($line =~ /pep/) {
    next if($line =~ /RawFile/);
    @inf = split /[ \n\r\t]+/, $line;
    $inf[1] =~ s/\.FT2$//;
    $mod = "$inf[1].$inf[2].$inf[2].$inf[4]\t$inf[3]";
    $comp2pep{$comp_pro} .= "$mod#";
  }
  else { die "ERROR: invalid line $line"; }
}
close FH;

foreach(keys %pro2enr) {
  next if($pro2enr{$_} !~ /#/);
  @ratios = split /#/, $pro2enr{$_};
  for($i = 0; $i < @ratios-1; $i++) {
    for($j = $i+1; $j < @ratios; $j++) {
      $pair1 = "$_#$ratios[$i]";
      $pair2 = "$_#$ratios[$j]";
      print STDERR "Processing pair $pair1 and $pair2 ...\n";

      # Write ProRata Xml config file
      open FH, ">$config" or die "couldn't open config file $config for writing";
      for($k = 0; $k < @prorata; $k++) {
        if($k == 33) {
          $rtmp = sprintf "%.3f", 1.0 - $ratios[$i];
          print FH "                <NATURAL>\t$rtmp,\t$ratios[$i]\t<\/NATURAL>\n";
        }
        elsif($k == 34) {
          $rtmp = sprintf "%.3f", 1.0 - $ratios[$j];
          print FH "                <ENRICHED>\t$rtmp,\t$ratios[$j]\t<\/ENRICHED>\n";
        }
        else { print FH "$prorata[$k]"; }
      }
      close FH;
     
      # Write DTA Select file
      open FH, ">$dtafilter" or die "couldn't open dta filter file $dtafilter for writing";
      for($k = 0; $k < 29; $k++) { print FH $dtalines[$k]; }
      @dinf = split /[ \n\r\t]+/, $dtalines[29];
      $dinf[0] = "$_";
      for($k = 0; $k < @dinf; $k++) {
        if($k == 0) { print FH $dinf[$k]; }
        else { print FH "\t$dinf[$k]"; }
      }
      print FH "\n";
      @mod1 = split /#/, $comp2pep{$pair1};
      for($k = 0; $k < @mod1; $k++) {
        @inf = split /\t/, $mod1[$k];
        print FH "*\t$inf[0]\t5.0\t0.5\t2000.00\t2000.00\t5000.0\t1\t500.0\t60.0\t1\tR.$inf[1].G\n";
      }
      @mod2 = split /#/, $comp2pep{$pair2};
      for($k = 0; $k < @mod2; $k++) {
        @inf = split /\t/, $mod2[$k];
        print FH "*\t$inf[0]\t5.0\t0.5\t2000.00\t2000.00\t5000.0\t1\t500.0\t60.0\t1\tR.$inf[1].G\n";
      }
      for($k = 31; $k < @dtalines; $k++) { print FH $dtalines[$k]; }

    # Sic Forma and Pratio
    chdir($wd);
    $rc = 0xffff & system "sicforma.exe";
    if($rc != 0) { die "error in sicforma\n"; }
    $rc = 0xffff & system "pratio.exe";
    if($rc != 0) { die "error in pratio\n"; }
    chdir(".");
 
    # parse output
    $state = 0;
    open FH, $pro_out or die "couldn't open $pro_out for reading";
    while($line = <FH>) {
      if($line =~ /^locus/) { $state = 1; }
      elsif($state == 1) {
        @inf = split /[ \n\t\r]+/, $line;
        print OH "$_\t$ratios[$i]\t$ratios[$j]\t$inf[1]\t$inf[2]\t$inf[3]\t$inf[4]\n";
        $state = 0;
      }
    }
    close FH;

    sleep(1);

    # Clean up files
    unlink("$wd/ProRataConfig.xml");
    unlink("$wd/DTASelect-filter.txt");
    unlink("$wd/ProRata_Quantification.qpr.xml");
    unlink("$wd/ProRata_Quantification_Peptide.txt");
    unlink("$wd/ProRata_Quantification_Protein.txt");
    opendir DH, "$wd/xic" or die "couldnt open dir $wd/xic";
    while($dline = readdir(DH)) {
      next if($dline =~ /^\./);
      unlink("$wd/xic/$dline");
    }
    closedir DH;
    rmdir("$wd/xic"); 

    }
  }
}

# unlink($config);
# unlink($dtafilter);
close OH;
