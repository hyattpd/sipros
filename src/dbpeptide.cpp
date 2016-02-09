/********************************************************/
// db_peptide.h
//
// SIPROS written by Doug Hyatt and Chongle Pan
//
// Detailed classes for handling associations between a
// database of peptides and associated MS2 scan data.
/********************************************************/

#include "dbpeptide.h"

// Constructor

DatabasePeptide::DatabasePeptide(class Text *txt, int p, int l, string ps, bool isp) {
  database_protein_ = txt;
  database_position_ = p;
  length_ = l;
  sequence_ = ps;
  is_ptm_ = isp;
}

// Destructor:  Frees all memory allocated for classes and structures.

void DatabasePeptideList::Clear() {
  while(!db_peptides_.empty()) {
    delete db_peptides_.back();
    db_peptides_.pop_back();
  }
  while(!db_peptide_scans_.empty()) {
    delete db_peptide_scans_.back();
    db_peptide_scans_.pop_back();
  }
}

// This function takes scan information (file number and scan number)
// and an existing DatabasePeptide class, and makes a "dbpepscan" structure
// containing all the information.

void DatabasePeptideList::AddPeptideAndScan(DatabasePeptide *dp, int fn, int sn) {
  struct _dbpepscan *dbps_pointer = new struct _dbpepscan;
  dbps_pointer->db_peptide = dp;
  dbps_pointer->scan_file_number = fn;
  dbps_pointer->scan_number = sn;
  dbps_pointer->score = 0.0;
  db_peptide_scans_.push_back(dbps_pointer);
}

// Simple sorting function for dbpepscan structures.
void DatabasePeptideList::SortDBPepScans() {
  sort(db_peptide_scans_.begin(), db_peptide_scans_.end(), DBPepScanCompare);
}

// Comparison function for two dbpepscan structures.  Sorts first by
// scan file number, then by scan number, and finally by spectrum score.

bool DBPepScanCompare(struct _dbpepscan *s1, struct _dbpepscan *s2) {
  return ( (s1->scan_file_number < s2->scan_file_number) || 
           (s1->scan_file_number == s2->scan_file_number && s1->scan_number < s2->scan_number) ||
           (s1->scan_file_number == s2->scan_file_number && s1->scan_number == s2->scan_number && s1->score > s2->score) );
}

// Processing function for a list of database peptide scans.  Since all
// dbpepscan structures have been sorted by scan number, we only have to 
// read in the associated scan data once for each scan (whenever we see
// a new scan number).  Each dbpepscan is scored via the preliminary scoring
// function.  If it passes the preliminary scoring, a final scoring function
// is run.  

void DatabasePeptideList::Process(ScanIndexSet &sindexset, double threshold) {
  vector<struct _dbpepscan *>::iterator viter;
  struct _dbpepscan *dbps_pointer;
  PeptideScorer *ps[MAX_CHARGE];
  MS2Scan ms2_data;
  double tmp_score = 0.0;
  int tmp_charge = 0;
  string tmp_peptide;

  cout << "Beginning process... # of db_peptides " << db_peptides_.size() << " # of db_peptide_scans " << db_peptide_scans_.size() << endl; fflush(stdout);
  SortDBPepScans();
  int old_file_number = -1, old_scan_number = -1, fn = -1, sn = -1;
  for(int i = 0; i < MAX_CHARGE; i++) { ps[i] = new PeptideScorer; }
  int tmp_counter = 0;

  // Main processing loop
  for(viter = db_peptide_scans_.begin(); viter != db_peptide_scans_.end(); viter++) {

    tmp_counter++;
    if(tmp_counter % 10000 == 0) { cout << "done with " << tmp_counter << " scores..." << endl; fflush(stdout); }

    dbps_pointer = *viter;

    fn = dbps_pointer->scan_file_number;
    sn = dbps_pointer->scan_number;
    if(fn != old_file_number || sn != old_scan_number) {
      if(ms2_data.ReadSingleScan(*(sindexset[fn]), sn) == false) {
        cerr << "Failed to read MS2 scan file#: " << fn << " scan#: " << sn << endl;
        exit(12);
      }
      if(ms2_data.HasAmbiguousCharge() == false && 
         ps[ms2_data.parent_charge()-1]->setMS2(ms2_data.scan_mzs(), ms2_data.scan_intensities(), ms2_data.scan_charges(), ms2_data.parent_mz(), ms2_data.parent_charge()) == false) {
        cerr << "Failed to initialize MS2 scan " << sn << endl;
        exit(12);
      } 
      else if(ms2_data.HasAmbiguousCharge() == true) {
        for(int j = 0; j < MAX_AMBIGUOUS_CHARGE; j++) {
          if(ps[j]->setMS2(ms2_data.scan_mzs(), ms2_data.scan_intensities(), ms2_data.scan_charges(), ms2_data.MzToMass()*(double)(j+1) - ((double)(j+1))*1.007276466, j+1) == false) {
            cerr << "Failed to initialize MS2 scan " << sn << endl;
            exit(12);
          }
        }
      }
    }
    old_file_number = fn; old_scan_number = sn;
// cout << "[" << dbps_pointer->db_peptide->sequence() << "] score " << tmp_score << " " << fn << " " << sn << endl; fflush(stdout);

/* CARBON skip preliminary scoring 
    if(ms2_data.ScorePeptidePreliminary(ps, dbps_pointer->db_peptide->sequence(), tmp_score, tmp_charge) == false) {
      cerr << "Failed to score scan #: " << sn << " in file# " << fn << endl;
      exit(13);
    } 
// cout << "[" << dbps_pointer->db_peptide->sequence() << "] score " << tmp_score << " " << fn << " " << sn << endl; fflush(stdout);
    if(tmp_score < threshold) {
      dbps_pointer->score = tmp_score;
      dbps_pointer->scan_charge = tmp_charge;
      continue;
    }
CARBON */

    if(ms2_data.ScorePeptideFinal(ps, dbps_pointer->db_peptide->sequence(), tmp_score, tmp_charge) == false) {
      cerr << "Failed to score scan #: " << sn << " in file# " << fn << endl;
      exit(13);
    } 

    dbps_pointer->score = tmp_score;
    dbps_pointer->scan_charge = tmp_charge;
  }
  for(int i = 0; i < MAX_CHARGE; i++) { delete ps[i]; }
cout << "End of process... # of db_peptides " << db_peptides_.size() << " # of db_peptide_scans " << db_peptide_scans_.size() << endl; fflush(stdout);
}

// Print out info
void DatabasePeptideList::Dump(ScanIndexSet &sind, double threshold) {
  int ptm_mod = 0;
  SortDBPepScans();
  for(size_t i = 0; i < db_peptide_scans_.size(); i++) {
    if(db_peptide_scans_[i]->score < threshold) continue;
    if(db_peptide_scans_[i]->db_peptide->is_ptm() == 1) ptm_mod = 1;
    cout << sind[db_peptide_scans_[i]->scan_file_number]->scan_file_name() << "\t";
    cout << db_peptide_scans_[i]->scan_number << "\t";
    cout << db_peptide_scans_[i]->score << "\t";
    cout << db_peptide_scans_[i]->db_peptide->sequence() << "\t";
    cout << db_peptide_scans_[i]->db_peptide->database_protein()->Substring(db_peptide_scans_[i]->db_peptide->database_position()-1, 
            db_peptide_scans_[i]->db_peptide->length()) << "\t";
    cout << db_peptide_scans_[i]->db_peptide->database_id() << "\t";
    cout << db_peptide_scans_[i]->db_peptide->database_position() << "\t"; 
    cout << endl;
    fflush(stdout);
  }
}

// Function to trim to the top X hits for each scan number.  Retains
// extra hits if they are tied with the last hit.  e.g. if you specify
// top 10, and #10-20 all have the same score, it will retain 20 hits.

void DatabasePeptideList::ReduceToTopXHits(int threshold) {
  SortDBPepScans();  
  int fn = -1, sn = -1, ctr = 0;
  double old_score = -100.0;
  struct _dbpepscan *dbsp_pointer;
  for(vector<struct _dbpepscan *>::iterator viter = db_peptide_scans_.begin(); viter != db_peptide_scans_.end(); viter++) {
    dbsp_pointer = *viter;
    if(dbsp_pointer->scan_file_number != fn || dbsp_pointer->scan_number != sn) ctr = 0;
    else ctr++;
    if(ctr >= threshold && dbsp_pointer->score < old_score) {
      db_peptide_scans_.erase(viter);
      viter--;
      continue;
    }
    fn = dbsp_pointer->scan_file_number;
    sn = dbsp_pointer->scan_number;
    old_score = dbsp_pointer->score;
  }
}
