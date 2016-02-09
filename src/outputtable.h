/********************************************************/
// outputtable.h
//
// SIPROS written by Doug Hyatt and Chongle Pan
//
// Class holding the final output of SIPROS
/********************************************************/

#ifndef _OUTPUTTABLE_H
#define _OUTPUTTABLE_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "text.h"
#include "dbpeptide.h"
#include "profrag.h"
#include "scanindex.h"

using namespace std;

// Final output class - single entry

class Output {

  int scan_file_number_;
  int scan_number_;
  double score_;

  string ms2_scan_sequence_;
  string database_sequence_;
  string database_id_;  
  int database_position_;

public:

  Output() {}
  ~Output() {}

  int scan_file_number() { return scan_file_number_; }
  int scan_number() { return scan_number_; }
  double score() { return score_; }
  string ms2_scan_sequence() { return ms2_scan_sequence_; }
  string database_sequence() { return database_sequence_; }
  string database_id() { return database_id_; }
  int database_position() { return database_position_; }

  void set_scan_file_number(int i) { scan_file_number_ = i; }
  void set_scan_number(int i) { scan_number_ = i; }
  void set_score(double d) { score_ = d; }
  void set_ms2_scan_sequence(string s) { ms2_scan_sequence_ = s; }
  void set_database_sequence(string s) { database_sequence_ = s; }
  void set_database_id(string s) { database_id_ = s; }
  void set_database_position(int i) { database_position_ = i; }
};

// Final output class - table of all entries
class OutputTable {

  vector <class Output *> entries_;

public:

  OutputTable() {}
  ~OutputTable() { Clear(); }

  vector <class Output *> entries() { return entries_; }

  void set_entries(vector <class Output *> vo) { entries_ = vo; }

  void Clear();
  void SortTable();
  void Update(class DatabasePeptideList &, double);

};

bool OutputCompare(class Output *, class Output *);

#endif
