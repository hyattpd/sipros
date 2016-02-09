/********************************************************/
// text.h
//
// SIPROS written by Doug Hyatt and Chongle Pan
//
// Text class for use in pattern matching algorithms.
// The "text" is the object being searched to see if
// "pattern" can be found within it.  This file also 
// includes a FASTA database class (implemented as a 
// collection of text classes).
/********************************************************/

#ifndef _TEXT_H
#define _TEXT_H

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

// A simple text class implementation, contains an
// identifier and a string representing the actual text
// to be searched.

class Text {
  string id_;      // Text identifier for this object.
  string word_;    // The actual string.

public:

  Text() { word_ = ""; id_ = ""; }
  Text(string s) { word_ = s; id_ = ""; }
  Text(string s1, string s2) { word_ = s1; id_ = s2; }
  ~Text() { }

  string id() { return id_; }
  string word() { return word_; }

  void set_id(string s) { id_ = s; }
  void set_word(string s) { word_ = s; } 

  void Append(string s) { word_.append(s); }
  int Length() { return word_.length(); }
  string Substring(int p, int l) { return word_.substr(p, l); }
};

// FASTA Database class implemented as a collection of
// Text classes.  The ids in the texts correspond to the
// FASTA headers.

class FastaDB {

  vector <Text *> sequences_;    // Vector of sequences

public:

  FastaDB(string);
  ~FastaDB() { for(size_t i = 0; i < sequences_.size(); i++) delete sequences_[i]; }

  vector <Text *> sequences() { return sequences_; }

  void set_sequences(vector <Text *> vt) { sequences_ = vt; }

  int ReadDatabase(string);  
  void PrintSequences(int);  

  int Size() { return sequences_.size(); }
  Text* operator[] (int ndx) { return sequences_[ndx]; }
};

#endif
