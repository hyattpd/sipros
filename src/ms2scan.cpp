/********************************************************/
// ms2scan.h
//
// SIPROS written by Doug Hyatt and Chongle Pan
//
// General class for handling MS2 Scan Data
/********************************************************/

#include "ms2scan.h"

// Constructor

MS2Scan::MS2Scan() {
  parent_mz_ = 0.0;
  parent_charge_ = 0;
}

// Given a scan number, looks the scan number up in
// the Scan Index class, takes the returned file offset,
// seeks to that position in the FT2 file, and finally 
// reads in the information associated with that scan
// number.

bool MS2Scan::ReadSingleScan(ScanIndex & sindex, int snum) {
  vector <string> scan_data;
  istringstream input;
  double tmp_double, tmp_charge;
  scan_mzs_.clear();
  scan_intensities_.clear();
  scan_charges_.clear();
  parent_mz_ = 0.0;
  parent_charge_ = 0;
  if(sindex.GetScanText(snum, scan_data) == false) return false;
  vector <string>::iterator viter;
  for(viter = scan_data.begin(); viter != scan_data.end(); viter++) {
    TokenVector words((*viter), " \t\n\r");
    if(words.size() == 4 && words[0][0] == 'S') { 
      istringstream input(words[3]);
      input >> parent_mz_;
    }
    if(words.size() == 3 && words[0][0] == 'Z') { 
      istringstream input(words[1]);
      input >> tmp_charge;
      parent_charge_ = (int)tmp_charge;
      if(parent_charge_ > MAX_CHARGE) {
        cerr << "Warning:  charge " << parent_charge_ << " too high." << endl;
        return false;
      }
    }
    if(words.size() == 6 && words[0][0] >= '0' && words[0][0] <= '9') {
      istringstream input;
      input.str(words[0]);
      input >> tmp_double;
      scan_mzs_.push_back(tmp_double);
      input.clear();
      input.str(words[1]);
      input >> tmp_double;
      scan_intensities_.push_back(tmp_double);
      input.clear();
      input.str(words[5]);
      input >> tmp_double;
      scan_charges_.push_back((int)tmp_double);
    }
  }
  if(parent_charge_ > 0) { parent_mz_ = parent_mz_ * (double)parent_charge_ - (double)parent_charge_*1.007276466; }

/*
  cout << "parent mz is " << parent_mz_ << " and charge is " << parent_charge_ << endl;
  for(size_t i = 0; i < scan_mzs_.size(); i++) {
    cout << "mz " << scan_mzs_[i] << "\tint " << scan_intensities_[i] << "\tchg " << scan_charges_[i] << endl;
  }
*/

  return true;
}

// Preliminary Scoring Function for Scans - Does a Fast but Low Quality
// scoring of a single MS2 scan.  Due to charges of 0 in the data, we maintain
// an array of Peptide Scorer classes for all potential charges.  The scoring
// function is run on all possible charges, and the highest result is determined
// to be the "real" charge.  In this case, this determined charge is stored in
// the passed parent_charge variable (replacing the value of 0 for that scan).
// The score is also stored in a passed variable.  Function returns true or false
// upon success/failure, respectively.

bool MS2Scan::ScorePeptidePreliminary(PeptideScorer *ps[MAX_CHARGE], string peptide, double &score, int &passed_charge) {
  double tmp_double;
  if(parent_charge_ != 0) {
    score = ps[parent_charge_-1]->prelimaryScorePeptide2(peptide);
    passed_charge = parent_charge_;
  } 
  else {
    double max_score = -10000.0;
    for(int i = 0; i < MAX_AMBIGUOUS_CHARGE; i++) {
      tmp_double = ps[i]->prelimaryScorePeptide2(peptide);
      if(tmp_double > max_score) {
        score = tmp_double;
        passed_charge = i+1;
        max_score = score;
      }
    }
  }
  return true;
}

// Final Scoring Function for Scans - Does a rigorous high quality 
// scoring of a single MS2 scan.  Due to charges of 0 in the data, we maintain
// an array of Peptide Scorer classes for all potential charges.  The scoring
// function is run on all possible charges, and the highest result is determined
// to be the "real" charge.  In this case, this determined charge is stored in
// the passed parent_charge variable (replacing the value of 0 for that scan).
// The score is also stored in a passed variable.  Function returns true or false
// upon success/failure, respectively.

bool MS2Scan::ScorePeptideFinal(PeptideScorer *ps[MAX_CHARGE], string peptide, double &score, int &passed_charge) {
  double tmp_double;
  if(parent_charge_ != 0) {
    score = ps[parent_charge_-1]->primaryScorePetide(peptide);
    passed_charge = parent_charge_;
  } 
  else {
    double max_score = -10000.0;
    for(int i = 0; i < MAX_AMBIGUOUS_CHARGE; i++) {
      tmp_double = ps[i]->primaryScorePetide(peptide);
      if(tmp_double > max_score) {
        score = tmp_double;
        passed_charge = i+1;
        max_score = score;
      }
    }
  }
  return true;
}

double MS2Scan::MzToMass() {
  if(parent_charge_ == 0) return parent_mz_;
  else return (parent_mz_ + (double)parent_charge_*1.007276466)/(double)parent_charge_;
}
