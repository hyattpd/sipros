/********************************************************/
// denovo.h
//
// SIPROS written by Doug Hyatt and Chongle Pan
//
// Class to process DeNovo Sequence tags, search them
// against a database (allowing for polymorphisms), and
// return dbpepscans (see dbpeptide.h) for final scoring.
/********************************************************/

#ifndef _DENOVO_H
#define _DENOVO_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include "tokenvector.h"
#include "bitap.h"
#include "alphabet.h"
#include "pattern.h"
#include "profrag.h"
#include "dbpeptide.h"
#include "blosum.h"

using namespace std;

// Class for DeNovo sequencing tags

class DeNovoTag {

  string file_name_;
  int scan_number_;
  double parent_mz_;
  int charge_;
  string index_;
  double left_flanking_mass_;
  double right_flanking_mass_;
  string sequence_tag_;
  double score_;
  int set_size_; 
  ifstream infstream_;

public:

  DeNovoTag();
  ~DeNovoTag() {}

  string file_name() { return file_name_; }
  int scan_number() { return scan_number_; }
  double parent_mz() { return parent_mz_; }
  int charge() { return charge_; }
  string index() { return index_; }
  double left_flanking_mass() { return left_flanking_mass_; }
  double right_flanking_mass() { return right_flanking_mass_; }
  string sequence_tag() { return sequence_tag_; }
  double score() { return score_; }
  int set_size() { return set_size_; }

  void set_file_name(string s) { file_name_ = s; }
  void set_scan_number(int i) { scan_number_ = i; }
  void set_parent_mz(double d) { parent_mz_ = d; }
  void set_charge(int i) { charge_ = i; }
  void set_index(string s) { index_ = s; }
  void set_left_flanking_mass(double d) { left_flanking_mass_ = d; }
  void set_right_flanking_mass(double d) { right_flanking_mass_ = d; }
  void set_sequence_tag(string s) { sequence_tag_ = s; }
  void set_score(double d) { score_ = d; }
  void set_set_size(int i) { set_size_ = i; }

  bool OpenFile(string);
  int ReadNextEntry();
  void FlipTag();
  void Process(FastaDB &, ProteinFragmentRules &, Amino *, DatabasePeptideList &, ScanIndexSet &);
};

#endif
