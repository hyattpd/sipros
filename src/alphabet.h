/********************************************************/
// alphabet.h 
//
// SIPROS written by Doug Hyatt and Chongle Pan        
//
// Simple class for protein and nucleic acid alphabets.
// Handles isomorphic characters (I = L, etc.), mass
// information for each peptide, and conversion between
// letters and numbers.
/********************************************************/

#ifndef _ALPHABET_H
#define _ALPHABET_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "proNovoConfig.h"

using namespace std;

class Alphabet {

  int size_;                   // # of characters in the alphabet
  vector < char > letters_;    // Vector of characters in the alphabet
  map < char, int > index_;    // Map for converting a character in the
                               // alphabet to a numerical index.

public:

  Alphabet();
  Alphabet(string s) { CreateAlphabetFromString(s); }
  ~Alphabet() {}

  int size() { return size_; }
  vector <char> letters() { return letters_; }
  map <char, int> index() { return index_; }

  void set_size(int s) { size_ = s; }
  void set_letters(vector <char> vc) { letters_ = vc; }
  void set_index(map <char, int> mci) { index_ = mci; }
 
  void CreateAlphabetFromString(string);
  void AddSingleLetter(char, int); 
  void SetEquivalentCharacters(char, char);

  int LetterToNumber(char c) { return index_[c]; }
  char NumberToLetter(int i) { return letters_[i]; }

  void PrintInfo();
};

// Amino Acid Alphabet:  This class is an instance of the
// alphabet class specific to proteins.  It sets I,J,L all
// equal to each other and is also capable of returning the
// monoisotopic mass for any letter (residue) in the alphabet.

class Amino : public Alphabet {

  map<char, double> monoisotopic_mass_;  // Map that returns a mass for a letter

public:
  Amino();
  ~Amino() {}

  map <char, double> monoisotopic_mass() { return monoisotopic_mass_; }

  void set_monoisotopic_mass(map <char, double> mcd) { monoisotopic_mass_ = mcd; }

  void EquateLeucineWithIsoleucine();
  double MonoisotopicMass(char c) { return monoisotopic_mass_[c]; }
};

// DNA Alphabet: Simple DNA alphabet (instance of alphabet
// class.

class DNAAlphabet : public Alphabet {

public:
  DNAAlphabet();
  ~DNAAlphabet() {}

};

#endif
