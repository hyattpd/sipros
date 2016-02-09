/********************************************************/
// blosum.h
//
// SIPROS written by Doug Hyatt and Chongle Pan
//
// Simple implementation of the BLOSUM matrix from BLAST
/********************************************************/

#include "blosum.h"

// Constructor for Blosum matrix.  Sets residues to their
// appropriate column/row #.

Blosum62::Blosum62() {
  letter_to_number_['A'] = 0; letter_to_number_['R'] = 1; letter_to_number_['N'] = 2; letter_to_number_['D'] = 3; letter_to_number_['C'] = 4; 
  letter_to_number_['Q'] = 5; letter_to_number_['E'] = 6; letter_to_number_['G'] = 7; letter_to_number_['H'] = 8; letter_to_number_['I'] = 9; 
  letter_to_number_['L'] = 10; letter_to_number_['K'] = 11; letter_to_number_['M'] = 12; letter_to_number_['F'] = 13; letter_to_number_['P'] = 14; 
  letter_to_number_['S'] = 15; letter_to_number_['T'] = 16; letter_to_number_['W'] = 17; letter_to_number_['Y'] = 18; letter_to_number_['V'] = 19; 
  letter_to_number_['B'] = 20; letter_to_number_['J'] = 21; letter_to_number_['Z'] = 22; letter_to_number_['X'] = 23; letter_to_number_['*'] = 24; 
  bitscore_ = 0; positives_ = 0;
}

// Scores a pair of peptides using the BLOSUM62 matrix.

int Blosum62::ScorePair(string s1, string s2) {
  bitscore_ = 0; positives_ = 0;
  if(s1.length() != s2.length()) return -1;
  int tmp_string;
  for(int i = 0; i < (int)(s1.length()); i++) {
    tmp_string = blos62[letter_to_number_[s1[i]]][letter_to_number_[s2[i]]];
    if(tmp_string > 0) positives_++;
    bitscore_ += tmp_string;
  }
  return 0;
}
