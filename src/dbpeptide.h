/********************************************************/
// dbpeptide.h
//
// SIPROS written by Doug Hyatt and Chongle Pan
//
// Detailed classes for handling associations between a
// database of peptides and associated MS2 scan data.
/********************************************************/

#ifndef _DBPEPTIDE_H
#define _DBPEPTIDE_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "text.h"
#include "ms2scan.h"
#include "PeptideScorer.h"

using namespace std;

// The database peptide class stores information about
// an individual peptide and its associated peptide in
// a database of proteins.  There are two peptide sequences
// stored in this class:  the sequence_ variable contains
// the observed sequence (derived from scan data), and modified
// from the original database by either a polymorphism or a PTM;
// and the Text id variable, which points to the protein 
// containing the sequence from the FASTA protein database. 
// The is_ptm variable indicates whether or not this is a post-
// translational modification or a single amino acid polymorphism.

class DatabasePeptide {

  class Text *database_protein_; // Pointer to the protein sequence/id
  int database_position_;      // Position within the protein where this sequence occurs
  int length_;                 // Length of the original peptide minus additional chars from ptms/etc.
  string sequence_;            // Modified sequence: contains the ptm/polymorphism
  bool is_ptm_;                // 0 = polymorphism, 1 = ptm

public:

  DatabasePeptide(class Text *, int, int, string, bool);
  ~DatabasePeptide() {}

  class Text *database_protein() { return database_protein_; }
  int database_position() { return database_position_; }
  int length() { return length_; }
  string sequence() { return sequence_; }
  bool is_ptm() { return is_ptm_; }
 
  void set_database_protein(class Text *tp) { database_protein_ = tp; }
  void set_database_position(int i) { database_position_ = i; }
  void set_length(int l) { length_ = l; }
  void set_sequence(string s) { sequence_ = s; }
  void set_is_ptm(bool b) { is_ptm_ = b; }

  // Quick Access Functions to Text attributes
  string database_sequence() { return database_protein_->Substring(database_position_-1, length_); }
  string database_id() { return database_protein_->id(); }
};

// Each Database Peptide class contains a single peptide sequence
// and its associated database data.  However, this peptide may be
// observed over many different scans.  In order to conserve memory,
// we do not store scan number information in the Database Peptide
// class, but instead create a new structure called "_dbpepscan".  This
// structure contains information about a single scan number, and a pointer
// to the DatabasePeptide class for that particular peptide sequence.
// In this way, we do not store database information for every single scan
// number, but only a pointer to that information.

struct _dbpepscan {
  DatabasePeptide *db_peptide;
  int scan_file_number;
  int scan_number;
  int scan_charge;
  double score;
};

// The DatabasePeptideList class is a list of all the observed
// peptides and the database entries from which those peptides were
// derived.  Also contains a list of struct dbpepscans containing all
// the associated scan numbers/charges/etc.

class DatabasePeptideList {

  vector <struct _dbpepscan *> db_peptide_scans_;
  vector <DatabasePeptide *> db_peptides_;

public:

  DatabasePeptideList() {}
  ~DatabasePeptideList() { Clear(); }

  vector <struct _dbpepscan *> db_peptide_scans() { return db_peptide_scans_; }
  vector <DatabasePeptide *> db_peptides() { return db_peptides_; }

  void AddPeptide(DatabasePeptide *dp) { db_peptides_.push_back(dp); }
  void AddPeptideAndScan(DatabasePeptide *, int, int);
  int Size() { return db_peptide_scans_.size(); }
  void Process(ScanIndexSet &, double);
  void Dump(ScanIndexSet &, double);

  void ReduceToTopXHits(int);
  void Clear();

  void SortDBPepScans();

};

bool DBPepScanCompare(struct _dbpepscan *, struct _dbpepscan *);

#endif
