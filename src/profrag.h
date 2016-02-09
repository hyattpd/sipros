/********************************************************/
// profrag.h
//
// SIPROS written by Doug Hyatt and Chongle Pan
//
// Class to generate tryptic peptides according to
// provided cleavage rules.  Also generates loops for
// single amino acid polymorphisms and single post-
// translational modifications.
/********************************************************/

#ifndef _PROFRAG_H
#define _PROFRAG_H

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "alphabet.h"
#include "text.h"
#include "ptm.h"
#include "dbpeptide.h"

using namespace std;

class ProteinFragmentRules {

  string after_residues_;             // String containing valid letters to follow a cleavage site
  string before_residues_;            // String containing valid letters to precede a cleavage site
  int max_missed_cleave_;             // Maximum number of missed cleavages in the peptide
  int min_specific_cleave_;           // Minimum number of specific (KR) cleavage sites in the peptide
  int min_peptide_length_;            // Minimum acceptable peptide length
  int max_peptide_length_;            // Maximum acceptable peptide length
  bool test_start_removal_;           // If true, allows the initial methionine to be cleaved
  double min_peptide_mass_;           // Minimum acceptable peptide mass
  double max_peptide_mass_;           // Maximum acceptable peptide mass
  int max_polymorphisms_;             // Maximum # of polymorphisms/PTMs allowed 
  int max_ptms_;                      // Maximum # of PTMs allowed 
  double parent_mass_error_;          // Upper bound for acceptable parent mass delta
  double fragment_mass_error_;        // Upper bound for fragment mass delta
  double left_offset_;                // Value to add to the parent mass from left side (default 0)
  double right_offset_;               // Value to add to the parent mass from right side (default 18.01505)
  string fragmentation_method_;       // CID, ESD, etc.

public:

  ProteinFragmentRules();
  ~ProteinFragmentRules() {}

  string after_residues() { return after_residues_; }  
  string before_residues() { return before_residues_; }
  int max_missed_cleave() { return max_missed_cleave_; }      
  int min_specific_cleave() { return min_specific_cleave_; } 
  int min_peptide_length() { return min_peptide_length_; }  
  int max_peptide_length() { return max_peptide_length_; }  
  bool test_start_removal() { return test_start_removal_; } 
  double min_peptide_mass() { return min_peptide_mass_; }    
  double max_peptide_mass() { return max_peptide_mass_; }  
  int max_polymorphisms() { return max_polymorphisms_; }    
  double parent_mass_error() { return parent_mass_error_; } 
  double fragment_mass_error() { return fragment_mass_error_; }   
  double left_offset() { return left_offset_; }          
  double right_offset() { return right_offset_; }       

  void set_after_residues(string s) { after_residues_ = s; }  
  void set_before_residues(string s) { before_residues_ = s; }
  void set_max_missed_cleave(int n) { max_missed_cleave_ = n; }      
  void set_min_specific_cleave(int n) { min_specific_cleave_ = n; } 
  void set_min_peptide_length(int n) { min_peptide_length_ = n; }  
  void set_max_peptide_length(int n) { max_peptide_length_ = n; }  
  void set_test_start_removal(bool b) { test_start_removal_ = b; } 
  void set_min_peptide_mass(double d) { min_peptide_mass_ = d; }    
  void set_max_peptide_mass(double d) { max_peptide_mass_ = d; }  
  void set_max_polymorphisms(int n) { max_polymorphisms_ = n; }    
  void set_parent_mass_error(double d) { parent_mass_error_ = d; } 
  void set_fragment_mass_error(double d) { fragment_mass_error_ = d; }   
  void set_left_offset(double d) { left_offset_ = d; }          
  void set_right_offset(double d) { right_offset_ = d; }       

  bool is_cleavage_site(char, char); 
  bool is_valid_peptide(string, int, int);

  bool populate_from_xml_config();

};

// This class actually creates peptides from a protein according
// to the cleavage rules in ProteinFragmentRules.  However, it
// also generates peptides with single polymorphisms or PTMs.
class ProteinFragment {

  ProteinFragmentRules *pr;
  Text *dbid;
  string protein;

  Amino aa;
  string residues;

  PTM_List *ptmlist;

  // Cleavage map and mass map
  vector <int> cl_map;       // Cleavage site or not for each position in protein
  vector <double> mass_map;  // Masses for each position in protein
  vector <double> msum_map;  // Mass sum ending at each position in protein

  // Loop control variables for polymorphisms
  bool done_poly;          // Finished evaluating protein?
  int mutpos;              // Position of the mutation (loop #1)
  int mutval;              // Value of the mutated residue (loop #2)
  int begin;               // Begin of the peptide (loop #3)
  int end;                 // End of the peptide (loop #4)

  // Loop control variables for PTMs
  bool done_ptm;               // Finished evaluating protein?
  int ptm_pos;                 // Position of the ptm (loop #1)
  int ptm_val;                 // Index of the ptm (loop #2)
  int ptm_loop_begin;          // Begin of the peptide (loop #3)
  int ptm_loop_end;            // End of the peptide (loop #4)

  //  values 
  string peptide;          // A string of the current peptide
  double mass;             // Mass of the current peptide
  int left_cleave;         // # of cleavage sites left of the mutation
  int right_cleave;        // # of cleavage sites right of the mutation
  int mut_cleave;          // Cleavage value of the mutated residue

public:

  ProteinFragment(ProteinFragmentRules *, PTM_List *);
  ~ProteinFragment() {}

  void set_protein(Text *);
  void increment_poly();
  void increment_ptm();
  bool is_valid_poly();
  bool is_valid_ptm();
  bool is_full_unmutated();
  bool has_polymorphism();
  bool finished_poly() { return done_poly; }
  bool finished_ptm() { return done_ptm; }

  double getmass() { return mass; }
 
  DatabasePeptide *get_dbpeptide_poly();
  DatabasePeptide *get_dbpeptide_ptm();

  void dumpinfo_poly();  
  void dumpinfo_ptm();  
};

#endif
