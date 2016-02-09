/********************************************************/
// text.cpp
//
// SIPROS written by Doug Hyatt and Chongle Pan
//
// Text class for use in pattern matching algorithms.
// The "text" is the object being searched to see if
// "pattern" can be found within it.  This file also
// includes a FASTA database class (implemented as a
// collection of text classes).
/********************************************************/

#include "text.h"

// FASTA Database constructor.  Reads from the file 
// specified in the supplied string.  An error in reading
// the file is fatal (we exit and throw an error).

FastaDB::FastaDB(string s) {
  if(ReadDatabase(s) == -1) {
    cerr << "Failed to read FASTA database" << s << "!" << endl;
    exit(1);
  }
}

// Function to read in a FASTA database.  Stores headers
// in the id field and sequences in the word field.  Prints
// a warning if any database line contains a non-A-Z character
// but doesn't return an error.  Only returns an error if there's
// a problem actually reading the database.

int FastaDB::ReadDatabase(string s) { 
  class Text *tx;
  char line[10000];
  ifstream db_file(s.c_str());
  if(db_file.is_open() == 0) return -1;
  while(!db_file.eof()) {
    db_file.getline(line, 10000);
    string sline(line);
    if(line[0] == '>') {
      tx = new Text;
      tx->set_id(sline.substr(1, sline.find_last_not_of("\r\n")));
      tx->set_word("");
      sequences_.push_back(tx);
    }
    else {
      if(sline.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ\r\n") != string::npos) {
        cerr << "WARNING: Invalid character in database line: " << endl << sline << endl;
      }
      sequences_[sequences_.size()-1]->Append(sline.substr(0, sline.find_last_not_of("\r\n")+1));
    } 
  }
  return 0;
}

// Prints the headers/sequences contained in the database.

void FastaDB::PrintSequences(int cpl) {
  for(int i = 0; i < (int)sequences_.size(); i++) {
    cout << ">" << sequences_[i]->id() << endl;
    for(int j = 0; j < (int)(sequences_[i]->Length()); j+=cpl) {
      if(j+cpl < (int)(sequences_[i]->Length())) cout << sequences_[i]->Substring(j, cpl) << endl; 
      else cout << sequences_[i]->Substring(j, sequences_[i]->Length()-j) << endl;
    }
  }
}
