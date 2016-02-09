/********************************************************/
// ms2scan.h
//
// SIPROS written by Doug Hyatt and Chongle Pan
//
// General class for handling MS2 Scan Data 
/********************************************************/

#ifndef _MS2SCAN_H_
#define _MS2SCAN_H_

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include "scanindex.h"
#include "PeptideScorer.h"
#include "tokenvector.h"

#define MAX_LINE_LEN 10000
#define MAX_CHARGE 100
#define MAX_AMBIGUOUS_CHARGE 3

using namespace std;

// MS2Scan Data Class:  Contains all information about
// m/z's, intensities, charges, etc.

class MS2Scan
{

  double parent_mz_;
  int parent_charge_;
  vector <double> scan_mzs_;
  vector <double> scan_intensities_;
  vector <int> scan_charges_; 

public:
  MS2Scan();
  ~MS2Scan() {}

  double parent_mz() { return parent_mz_; }
  int parent_charge() { return parent_charge_; }
  vector <double> scan_mzs() { return scan_mzs_; }
  vector <double> scan_intensities() { return scan_intensities_; }
  vector <int> scan_charges() { return scan_charges_; }

  void set_parent_mz(double d) { parent_mz_ = d; }
  void set_parent_charge(int i) { parent_charge_ = i; }
  void set_scan_mzs(vector <double> vd) { scan_mzs_ = vd; }
  void set_scan_intensities(vector <double> vd) { scan_intensities_ = vd; }
  void set_scan_charges(vector <int> vi) { scan_charges_ = vi; }

  bool ReadSingleScan(ScanIndex &, int);

  bool ScorePeptidePreliminary(PeptideScorer* [MAX_CHARGE], string, double &, int &);
  bool ScorePeptideFinal(PeptideScorer* [MAX_CHARGE], string, double &, int &);

  double MzToMass();
  bool HasAmbiguousCharge() { if(parent_charge_ == 0) return true; return false; }
};

#endif 
