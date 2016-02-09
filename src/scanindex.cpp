/********************************************************/
// scanindex.cpp
//
// SIPROS written by Doug Hyatt and Chongle Pan
//
// Implementation of random access to scans stored in
// the FT2 file format.  Contains two classes, a scan
// index class representing an index of a single FT2
// file, and scan index set, which implements a set of
// indices for multiple FT2 files.
/********************************************************/

#include "scanindex.h"
#include "tokenvector.h"

// This routine builds a scan index from an FT2 file
// and returns false if there are any errors.

bool ScanIndex::BuildIndex(const string fn) {
  if(fn != "NONE") scan_file_name_ = fn;
  ifstream scan_file_stream(scan_file_name_.c_str());
  if(scan_file_stream.is_open() == 0) return false; 
  char line[SCAN_LINE_LEN];
  int scan_number;
  streampos file_position = 0;
  while(!scan_file_stream.eof()) {
    scan_file_stream.getline(line, SCAN_LINE_LEN);
    if(line[0] == 'S') {
      TokenVector words(line, "\t\n\r ");
      istringstream string_input(words[1]);
      string_input >> scan_number;
      scan_number_to_file_position_[scan_number] = file_position;
// cout << "S\t" << scan_number << "\t" << file_position << endl;
    }
    file_position += scan_file_stream.gcount();
  }
  scan_file_stream.close();
  return true;
}

// This routine writes the current scan index to the
// supplied filename, returning false on any errors.

bool ScanIndex::WriteIndex(const string fn) {
  if(fn != "NONE") scan_file_name_ = fn;
  string index_file_name = scan_file_name_ + ".spx";
  ofstream indexFile(index_file_name.c_str());
  map<int, streampos>::iterator i;
  for(i = scan_number_to_file_position_.begin(); i != scan_number_to_file_position_.end(); i++) {
    indexFile.write((const char *)(&(i->first)), sizeof(int));
    if(indexFile.bad()) return false;
    indexFile.write((const char *)(&(i->second)), sizeof(streampos));
    if(indexFile.bad()) return false;
  }
  return true;
}

// This routine reads a scan index from a file and
// returns false on any errors.

bool ScanIndex::ReadIndex(const string fn) {
  if(fn != "NONE") scan_file_name_ = fn;
  string index_file_name = scan_file_name_ + ".file_positionx";
  ifstream index_file_stream(index_file_name.c_str());
  while(index_file_stream.good()) {
    int scan_number;
    streampos file_position;
    index_file_stream.read((char *)(&scan_number), sizeof(int));
    index_file_stream.read((char *)(&file_position), sizeof(streampos));
    if(index_file_stream.good()) scan_number_to_file_position_[scan_number] = file_position;	
  }	
  if(index_file_stream.bad() || !index_file_stream.eof()) return false;
  return true;	
}

// This routine takes a scan number and a reference to a file 
// file position, and stores the scan number's file position in
// the supplied reference.  Returns false if it can't find the scan
// number in the index.

bool ScanIndex::GetFilePosition(int scan_num, streampos & file_pos) {
  map<int, streampos>::iterator iter = scan_number_to_file_position_.find(scan_num);
  if(iter == scan_number_to_file_position_.end()) {
    cerr << "Scanparse Warning: Scan number not in index" << endl;
    return false; 
  }
  file_pos = iter->second;
  return true;
}

// This routine takes a scan number and a reference to a vector of
// strings, locates the scan number in the FT2 file, and puts the lines
// into the supplied vector.  Returns false on any error.

bool ScanIndex::GetScanText(int scan_num, vector <string> &scan_lines) {
  streampos file_pos;
  if(GetFilePosition(scan_num, file_pos) == false) return false;
  int saw_scan = 0;
  char line[SCAN_LINE_LEN];
  ifstream scan_file_stream(scan_file_name_.c_str());
  if(scan_file_stream.is_open() == 0) return false;
  scan_file_stream.seekg(file_pos);
  while(!scan_file_stream.eof()) {
    scan_file_stream.getline(line, SCAN_LINE_LEN);
    if(line[0] == 'S') {
      saw_scan++;
      if(saw_scan == 2) break;
    }
    scan_lines.push_back(line);
// cout << scan_num << "\t" << line << endl;
  }	
  scan_file_stream.close();
  return true;
}

// End of ScanIndex methods

// Scan Index Set Functions:  These do the same thing
// as their individual counterparts, just traversing a
// vector of indices to locate the correct entry.  See
// the individual functions in the scan index class for
// descriptions.

int ScanIndexSet::FileNumber(string fn) {
  for(size_t i = 0; i < scan_indices_.size(); i++) 
    if(scan_indices_[i]->scan_file_name() == fn) return i;
  return -1;
}

bool ScanIndexSet::BuildIndex(string fn) {
  int ndx = FileNumber(fn);
  if(ndx == -1) {
    ScanIndex *index_pointer = new ScanIndex;
    scan_indices_.push_back(index_pointer);
    return index_pointer->BuildIndex(fn);
  }
  return scan_indices_[ndx]->BuildIndex();
}

bool ScanIndexSet::WriteIndex(string fn) {
  int ndx = FileNumber(fn);
  if(ndx == -1) return false;
  return scan_indices_[ndx]->WriteIndex();
}

bool ScanIndexSet::ReadIndex(string fn) {
  int ndx = FileNumber(fn);
  if(ndx == -1) return false;
  return scan_indices_[ndx]->ReadIndex();
}

bool ScanIndexSet::GetFilePosition(string fn, int scan_num, streampos &file_pos) {
  int ndx = FileNumber(fn);
  if(ndx == -1) return false;
  return scan_indices_[ndx]->GetFilePosition(scan_num, file_pos);
}

bool ScanIndexSet::GetScanText(string fn, int scan_num, vector <string> &scan_lines) {
  int ndx = FileNumber(fn);
  if(ndx == -1) return false;
  return scan_indices_[ndx]->GetScanText(scan_num, scan_lines);
}
