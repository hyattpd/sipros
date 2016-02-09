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

#include "scanmass.h"

// Constructor for ScanMassList

ScanMassList::ScanMassList() {
  ft2_file_.clear();
  scan_mass_list_.clear();
}

// Destructor - frees up vectors
ScanMassList::~ScanMassList() {
  for(vector<ScanMass*>::iterator vp = scan_mass_list_.begin(); vp != scan_mass_list_.end(); vp++)
    delete(*vp);
}

// Function for reading in a single FT2 file, parsing all
// the scans, and storing a simple list of scan numbers with
// their parent masses.

bool ScanMassList::ReadFt2File(vector<string> s) {
  char line[MAX_LINE_LEN];  
  ScanMass *scan_mass_pointer;
  int tmp_scan = -1, tmp_charge = -1;
  double tmp_mass_sline, tmp_mass_zline = -1.0;  // masses from the two different lines in the ft2 file (the S line and the Z line)
  for(size_t i = 0; i < s.size(); i++) { ft2_file_.push_back(s[i]); }
  for(size_t i = 0; i < s.size(); i++) {
    ifstream ft2_stream(ft2_file_[i].c_str());
    if(ft2_stream.is_open() == 0) { return false; } 
    while(!ft2_stream.eof()) {
      ft2_stream.getline(line, MAX_LINE_LEN);
      if(line[0] == 'S') {
        TokenVector words(line, " \r\t\n");
        istringstream input(words[1]);
        input >> tmp_scan;
        input.clear();
        input.str(words[3]);
        input >> tmp_mass_sline; 
      } 
      if(line[0] == 'Z') {
        TokenVector words(line, " \r\t\n");
        istringstream input(words[1]);
        input >> tmp_charge;
        input.clear();
        input.str(words[2]);
        input >> tmp_mass_zline; 
        if(tmp_charge != 0) {
          scan_mass_pointer = new ScanMass(i, tmp_scan, tmp_mass_sline*(double)tmp_charge - (double)tmp_charge*1.007276466);
          scan_mass_list_.push_back(scan_mass_pointer);
          scan_mass_pointer = new ScanMass(i, tmp_scan, tmp_mass_sline*(double)tmp_charge - (double)tmp_charge*1.007276466 - 1.007276466);
          scan_mass_list_.push_back(scan_mass_pointer);
          scan_mass_pointer = new ScanMass(i, tmp_scan, tmp_mass_sline*(double)tmp_charge - (double)tmp_charge*1.007276466 - 2*1.007276466);
          scan_mass_list_.push_back(scan_mass_pointer);
          scan_mass_pointer = new ScanMass(i, tmp_scan, tmp_mass_sline*(double)tmp_charge - (double)tmp_charge*1.007276466 - 3*1.007276466);
          scan_mass_list_.push_back(scan_mass_pointer);
          scan_mass_pointer = new ScanMass(i, tmp_scan, tmp_mass_sline*(double)tmp_charge - (double)tmp_charge*1.007276466 + 1.007276466);
          scan_mass_list_.push_back(scan_mass_pointer);
          scan_mass_pointer = new ScanMass(i, tmp_scan, tmp_mass_sline*(double)tmp_charge - (double)tmp_charge*1.007276466 + 2*1.007276466);
          scan_mass_list_.push_back(scan_mass_pointer);
          scan_mass_pointer = new ScanMass(i, tmp_scan, tmp_mass_sline*(double)tmp_charge - (double)tmp_charge*1.007276466 + 3*1.007276466);
          scan_mass_list_.push_back(scan_mass_pointer);
        }
        else {
          for(int j = 1; j <= 3; j++) {
            scan_mass_pointer = new ScanMass(i, tmp_scan, tmp_mass_sline*(double)j - (double)j*1.007276466); 
            scan_mass_list_.push_back(scan_mass_pointer);
            scan_mass_pointer = new ScanMass(i, tmp_scan, tmp_mass_sline*(double)j - (double)j*1.007276466 - 1.007276466); 
            scan_mass_list_.push_back(scan_mass_pointer);
            scan_mass_pointer = new ScanMass(i, tmp_scan, tmp_mass_sline*(double)j - (double)j*1.007276466 - 2*1.007276466); 
            scan_mass_list_.push_back(scan_mass_pointer);
            scan_mass_pointer = new ScanMass(i, tmp_scan, tmp_mass_sline*(double)j - (double)j*1.007276466 - 3*1.007276466); 
            scan_mass_list_.push_back(scan_mass_pointer);
            scan_mass_pointer = new ScanMass(i, tmp_scan, tmp_mass_sline*(double)j - (double)j*1.007276466 + 1.007276466); 
            scan_mass_list_.push_back(scan_mass_pointer);
            scan_mass_pointer = new ScanMass(i, tmp_scan, tmp_mass_sline*(double)j - (double)j*1.007276466 + 2*1.007276466); 
            scan_mass_list_.push_back(scan_mass_pointer);
            scan_mass_pointer = new ScanMass(i, tmp_scan, tmp_mass_sline*(double)j - (double)j*1.007276466 + 3*1.007276466); 
            scan_mass_list_.push_back(scan_mass_pointer);
          }
        }
      }
    }
  }
  return true;
}

// Sort the Mass List numerically by parent mass value

void ScanMassList::SortMasses() {
  sort(scan_mass_list_.begin(), scan_mass_list_.end(), ScanMassCompare);
}

// Given a specified mass and an acceptable mass error, returns
// the lower and upper boundaries within the list of the set of all
// masses that fit the search criteria.

pair<int, int> ScanMassList::GetRangeFromMass(double target, double error) {
  pair <int, int> p;
  int low = 0, high = scan_mass_list_.size()-1, mid;
  double lb, ub;  // lower and upper bounds on acceptable parent mass values
  ub = target + error;
  lb = target - error;
  while((high-low) > 1) {
    mid = (high+low)/2;
    if(scan_mass_list_[mid]->parent_mass() > target) high = mid;
    else low = mid;
  }

  // Iterate till we get to the first element > than the lower bound
  int ndx = low;
  if(scan_mass_list_[ndx]->parent_mass() >= lb) {
    while(ndx > 0 && scan_mass_list_[ndx]->parent_mass() >= lb) { ndx--; }
    ndx++;
  }
  else {
    while(ndx < (int)scan_mass_list_.size() && scan_mass_list_[ndx]->parent_mass() < lb) { ndx++; }
  }

  if(ndx == (int)scan_mass_list_.size() || scan_mass_list_[ndx]->parent_mass() > ub) { p = make_pair(-1, -1); }
  else {
    low = ndx;
    while(ndx < (int)scan_mass_list_.size() && scan_mass_list_[ndx]->parent_mass() <= ub) { ndx++; }
    high = ndx-1;
    p = make_pair(low, high);
  }

  if(p.first == -1) return p;
  if(scan_mass_list_[p.first]->parent_mass() < lb) cerr << "ERROR L " << scan_mass_list_[p.first]->parent_mass() << " " << lb << endl;
  if(scan_mass_list_[p.second]->parent_mass() > ub) cerr << "ERROR U " << scan_mass_list_[p.second]->parent_mass() << " " << ub << endl;
  return p;
}

// Print function

void ScanMassList::Dump() {
  cout << scan_mass_list_.size() << endl;
  for(vector<ScanMass *>::iterator vp = scan_mass_list_.begin(); vp != scan_mass_list_.end(); vp++) {
    (*vp)->Print();
  }
}

bool ScanMassCompare(ScanMass *s1, ScanMass *s2) {
  return (s1->parent_mass() < s2->parent_mass());
}
