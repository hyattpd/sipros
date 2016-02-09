/********************************************************/
// bitap.h
//
// SIPROS written by Doug Hyatt and Chongle Pan
//
// Implementation of Uwe Manber's bitap algorithm for
// approximate string matching.  This file contains two
// classes.  A Match class representing a single match
// between a pattern and a text, and a Bitap class which
// consists of a list of hits along with the actual bitap
// methods for producing the list.  Requires the Pattern
// and Text classes (representing the string and the
// database) to work.
/********************************************************/

#ifndef _BITAP_H
#define _BITAP_H

#include "alphabet.h"
#include "pattern.h"
#include "text.h"

using namespace std;

// The match class contains a single match between a pattern
// and text, with position and mismatch information included.

class Match {

  Pattern *pattern_pointer_; // A pointer to the pattern we're looking for
  Text *text_pointer_;       // A pointer to the text we're searching in
  int position_;             // Position within the text of this pattern
  int mismatch_count_;       // Number of mismatches in the located pattern 

public:

  Match() {}
  ~Match() {}

  Pattern *pattern_pointer() { return pattern_pointer_; }
  Text *text_pointer() { return text_pointer_; }
  int position() { return position_; }
  int mismatch_count() { return mismatch_count_; }
 
  void set_pattern_pointer(Pattern *p) { pattern_pointer_ = p; }
  void set_text_pointer(Text *t) { text_pointer_ = t; }
  void set_position(int i) { position_ = i; }
  void set_mismatch_count(int i) { mismatch_count_ = i; }
   
  string MatchString();
  string TextString() { return text_pointer_->word(); }
  string PatternString() { return pattern_pointer_->word(); }

  void Print(bool = 0);
  char Residue(int pos) { return text_pointer_->word()[pos]; } 

};

// The bitap class contains a list of matches, and the
// methods for running the bitap algorithm on a supplied
// pattern and text.

class Bitap {

  vector <Match *> matches_;  // A set of all the matches

public:

  Bitap();
  ~Bitap() { for(size_t i = 0; i < matches_.size(); i++) delete matches_[i]; }

  vector <Match *> matches() { return matches_; }

  void set_matches(vector <Match *> vm) { matches_ = vm; }

  void KMismatch(Text *, Pattern *, Alphabet &, int);
  void Reset() { matches_.clear(); }
  int Size() { return matches_.size(); }

  void Print(bool = 0);

  Match * operator [] (int ndx) { return matches_[ndx]; }
};

#endif
