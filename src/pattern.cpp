/********************************************************/
// pattern.cpp
//
// SIPROS written by Doug Hyatt and Chongle Pan
//
// Pattern class for use in pattern matching algorithms
/********************************************************/

#include "pattern.h"

// Constructor: creates a pattern from a supplied
// string and alphabet, including the bitmasks.

Pattern::Pattern(string s, Alphabet *a) { 
  word_ = s;
  pattern_alphabet_ = a;
  for(int i = 0; i < a->size(); i++) letter_mask_.push_back((unsigned long)(~0)); 
  for(size_t i = 0; i < s.length() && i < 32; i++)
    letter_mask_[a->LetterToNumber(s[i])] &= ~(1UL << i);
}

Pattern::~Pattern() {}

// Print out info on the bitmasks for debugging purposes.

void Pattern::PrintMasks() {
  cout << "printing pattern masks..." << endl;
  for(int i = 0; i < pattern_alphabet_->size(); i++)  {
    cout << pattern_alphabet_->NumberToLetter(i) << "\t";
    PrintBitString(letter_mask_[i]);
  }
}

void Pattern::PrintBitString(unsigned long l) {
  for(int i = 31; i >= 0; i--) {
    if( ((1UL << i) & l) == (1UL << i) ) cout << "1";
    else cout << "0";
  }
  cout << endl;
}
