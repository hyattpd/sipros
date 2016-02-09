/********************************************************/
// pattern.h
//
// SIPROS written by Doug Hyatt and Chongle Pan
//
// Pattern class for use in pattern matching algorithms
/********************************************************/

#ifndef _PATTERN_H
#define _PATTERN_H

#include "alphabet.h"

using namespace std;

// Simple pattern class consisting of an identifier,
// a string of letters, an alphabet class to represent
// the valid letters that can be used in the string, and
// a vector of bitmaps for each letter in the alphabet,
// i.e. ABC generates 100, 010, and 001.

class Pattern {

  string id_;                            // Text identifier
  string word_;                          // The actual pattern
  vector <unsigned long> letter_mask_;   // Bitmap for each letter in the alphabet
  Alphabet *pattern_alphabet_;           // The alphabet of valid characters

public:

  Pattern(string, Alphabet *);
  ~Pattern();

  string id() { return id_; }
  string word() { return word_; }
  vector <unsigned long> letter_mask() { return letter_mask_; }
  Alphabet *pattern_alphabet() { return pattern_alphabet_; }

  void set_id(string s) { id_ = s; }
  void set_word(string s) { word_ = s; }
  void set_letter_mask(vector <unsigned long> vul) { letter_mask_ = vul; }
  void set_pattern_alphabet(Alphabet *a) { pattern_alphabet_ = a; }

  int Size() { return word_.size(); }
  void SetPattern(string, Alphabet *);
  unsigned long BitMask(int i) { return letter_mask_[i]; } 

  void PrintMasks();
  void PrintBitString(unsigned long);
};


#endif
