/********************************************************/
// scanindex.h
//
// SIPROS written by Doug Hyatt and Chongle Pan
//
// Implementation of random access to scans stored in
// the FT2 file format.  Contains two classes, a scan
// index class representing an index of a single FT2
// file, and scan index set, which implements a set of
// indices for multiple FT2 files.
/********************************************************/

#ifndef _SCAN_INDEX_H
#define _SCAN_INDEX_H

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>

#define SCAN_LINE_LEN 10000

using namespace std;

// Implementation of random access to scans contained
// in an FT2 file.  Contains a map of scan number to
// position within the file.

class ScanIndex 
{

  string scan_file_name_;
  vector <int> scan_numbers_;
  map <int, streampos> scan_number_to_file_position_;

public:

  ScanIndex() { scan_file_name_ = ""; }
  ~ScanIndex() {}

  string scan_file_name() { return scan_file_name_; }
  vector <int> scan_numbers() { return scan_numbers_; }
  map <int, streampos> scan_number_to_file_position() { return scan_number_to_file_position_; }

  void set_scan_file_name(string s) { scan_file_name_ = s; }
  void set_scan_numbers(vector <int> vi) { scan_numbers_ = vi; }
  void set_scan_number_to_file_position(map <int, streampos> mis) { scan_number_to_file_position_ = mis; }

  bool BuildIndex(const string fn = "NONE");
  bool ReadIndex(const string fn = "NONE");
  bool WriteIndex(const string fn = "NONE");

  bool GetFilePosition(int, streampos &);
  bool GetScanText(int, vector<string> &);

};

// Similar to ScanIndex, but implements random access over an
// entire set of FT2 files.  Implemented as a vector of ScanIndex
// classes with appropriate access routines.

class ScanIndexSet {

  vector <ScanIndex *> scan_indices_;

public:

  ScanIndexSet() {};
  ~ScanIndexSet() { for(size_t i = 0; i < scan_indices_.size(); i++) delete scan_indices_[i]; }

  vector <ScanIndex *> scan_indices() { return scan_indices_; }

  void set_scan_indices(vector <ScanIndex *> vsi) { scan_indices_ = vsi; }

  int FileNumber(string);

  bool BuildIndex(string);
  bool ReadIndex(string);
  bool WriteIndex(string);

  bool GetFilePosition(string, int, streampos &);
  bool GetScanText(string, int, vector<string> &);

  string GetScanFileName(int i) { return scan_indices_[i]->scan_file_name(); }

  ScanIndex * operator [] (int i) { return scan_indices_[i]; }
};

#endif /*_SCAN_INDEX_H*/
