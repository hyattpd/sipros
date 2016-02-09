/********************************************************/
// outputtable.h
//
// SIPROS written by Doug Hyatt and Chongle Pan
//
// Class holding the final output of SIPROS
/********************************************************/

#include "outputtable.h"

using namespace std;

// Erase all current output and clear the vectors

void OutputTable::Clear() {
  while(!entries_.empty()) {
    delete entries_.back();
    entries_.pop_back();
  }
}

// Takes a database peptide list and adds all results
// above a specified threshold score to the output table.

void OutputTable::Update(class DatabasePeptideList &dbpl, double thresh) {
cerr << "presort" << endl; fflush(stderr);
  dbpl.SortDBPepScans();
cerr << "postsort" << endl; fflush(stderr);
  for(size_t i = 0; i < dbpl.db_peptide_scans().size(); i++) {
cout << "raw dbpepscan " << i << endl; fflush(stdout);
    if(dbpl.db_peptide_scans()[i]->score < thresh) continue;
cout << "postloop dbpepscan " << i << endl; fflush(stdout);
    Output *tmp_output = new class Output(); 
cout << "post new " << i << endl; fflush(stdout);
    tmp_output->set_score(dbpl.db_peptide_scans()[i]->score);
cout << "post sscore " << i << endl; fflush(stdout);
    tmp_output->set_scan_file_number(dbpl.db_peptide_scans()[i]->scan_file_number);
cout << "post ssscan " << i << endl; fflush(stdout);
    tmp_output->set_scan_number(dbpl.db_peptide_scans()[i]->scan_number);
cout << "post ssscanno " << i << endl; fflush(stdout);
    tmp_output->set_ms2_scan_sequence(dbpl.db_peptide_scans()[i]->db_peptide->sequence());
cout << "post ms2set " << i << endl; fflush(stdout);
    tmp_output->set_database_sequence(dbpl.db_peptide_scans()[i]->db_peptide->database_sequence());
cout << "post dbseq " << i << endl; fflush(stdout);
    tmp_output->set_database_id(dbpl.db_peptide_scans()[i]->db_peptide->database_id());
cout << "post dbid " << i << endl; fflush(stdout);
    tmp_output->set_database_position(dbpl.db_peptide_scans()[i]->db_peptide->database_position());
cout << "post dbpos " << i << endl; fflush(stdout);
  }
}

void OutputTable::SortTable() {
  sort(entries_.begin(), entries_.end(), OutputCompare);
}

bool OutputCompare(class Output *s1, class Output *s2) {
  return ( (s1->scan_file_number() < s2->scan_file_number()) ||
           (s1->scan_file_number() == s2->scan_file_number() && s1->scan_number() < s2->scan_number()) ||
           (s1->scan_file_number() == s2->scan_file_number() && s1->scan_number() == s2->scan_number() && s1->score() > s2->score()) );
}
