/********************************************************/
// scanmass.h
//
// SIPROS written by Doug Hyatt and Chongle Pan
//
// Class to maintain a list of all masses associated 
// with all scans, and functions for quick lookup of 
// all scans within a certain distance (delta) of a
// specified mass.
/********************************************************/

#ifndef SCANMASS_H_
#define SCANMASS_H_

#include <string>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include "tokenvector.h"

#define MAX_LINE_LEN 10000

using namespace std;

// The "ScanMass" class consists of a single scan number and
// its associated parent mass.

class ScanMass
{
  int file_number_;
  int scan_number_;
  double parent_mass_;

public:
  ScanMass(int fn, int sn, double pm) { file_number_ = fn; scan_number_ = sn; parent_mass_ = pm; }
  ~ScanMass() {}

  int file_number() { return file_number_; }
  int scan_number() { return scan_number_; }
  double parent_mass() { return parent_mass_; }

  void set_file_number(int i) { file_number_ = i; }
  void set_scan_number(int i) { scan_number_ = i; }
  void set_parent_mass(double d) { parent_mass_ = d; }

  void Print() { cout << file_number_+1 << "#" << scan_number_ << "\t" << parent_mass_ << endl; }
};

// ScanMassList maintains a list of all scan numbers and their
// parent masses, and methods for returning all scan numbers within
// a specified distance from a supplied parent mass.

class ScanMassList {
  vector<ScanMass *> scan_mass_list_;
  vector<string> ft2_file_;

public:

  ScanMassList(); 
  ~ScanMassList(); 

  vector <ScanMass *> scan_mass_list() { return scan_mass_list_; }
  vector <string> ft2_file() { return ft2_file_; }

  int Size() { return scan_mass_list_.size(); }

  bool ReadFt2File(vector <string>);
  void SortMasses();

  pair<int, int> GetRangeFromMass(double, double);

  ScanMass *operator [] (int i) { return scan_mass_list_[i]; }

  void Dump();

};

bool ScanMassCompare(ScanMass *, ScanMass *);

#endif /*SCANMASS_H_*/
