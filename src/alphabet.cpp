/********************************************************/
// alphabet.cpp
//
// SIPROS written by Doug Hyatt and Chongle Pan
//
// Simple class for protein and nucleic acid alphabets.
// Handles isomorphic characters (I = L, etc.), mass
// information for each peptide, and conversion between
// letters and numbers.
/********************************************************/

#include "alphabet.h"

// Constructor: adds all ASCII chars to the alphabet
// as equivalent (index 0).

Alphabet::Alphabet() {
  for(int i = 0; i < 256; i++) AddSingleLetter(i, 0);
  size_ = 1;
}

// Creates an alphabet out of the characters contained
// in the supplied string, i.e. 'actg' for a nucleic acid
// alphabet, etc.  All other ASCII characters do remain in
// the alphabet with index position 0.  This is so any parser 
// is capable of handling any character it encounters, regardless
// of if it is included in the specified alphabet.

void Alphabet::CreateAlphabetFromString(string s) {
  for(int i = 0; i < 256; i++) AddSingleLetter(i, 0);
  for(size_t i = 0; i < s.length(); i++) AddSingleLetter(s[i], i+1);
  size_ = 0;
  for(int i = 0; i < 256; i++) if(index_[i] > size_) size_ = index_[i];
  size_ += 1;
}

// Adds a single letters_ to the alphabet, removing it from
// index 0 and assigning it the next available positive integer
// index.

void Alphabet::AddSingleLetter(char c, int i) {
  vector<char>::iterator viter;
  map<char, int>::iterator miter = index_.find(c);
  if(miter != index_.end()) {
    index_.erase(miter);
    for(viter = letters_.begin(); viter != letters_.end(); viter++) {
      if(*viter == c) letters_.erase(viter);
    }
  }
  letters_.push_back(c);
  index_[c] = i;
}

// Erases the second supplied character from the alphabet and
// adds it to the alphabet again with the same index as the first
// supplied character.  End result is that the two characters
// share the same index and are therefore equivalent.

void Alphabet::SetEquivalentCharacters(char c1, char c2) {
  if(index_[c2] >= 0 && index_[c2] < size_) AddSingleLetter(c1, index_[c2]);
}

// Prints the basic information included in the alphabet.  Useful
// for debugging purposes.

void Alphabet::PrintInfo() {
  for(int i = 0; i < (int)letters_.size(); i++) cout << letters_[i] << endl; 
  for(map<char,int>::iterator miter = index_.begin(); miter != index_.end(); miter++) cout << miter->first << "\t" << miter->second << endl;  
}

// Constructor for the amino acid alphabet.  Adds all 20 residues
// plus BJZX to the alphabet.  Adds monoisotopic mass information
// for the 20 residues plus J.

Amino::Amino() { 
  string s = "ACDEFGHIKLMNPQRSTVWYJ";
  string tmp_c;
  CreateAlphabetFromString(s); 
  for(size_t i = 0; i < s.size(); i++) {
    tmp_c = s[i];
    monoisotopic_mass_[s[i]] = ProNovoConfig::getResidueMass(tmp_c);
  }
/*
  monoisotopic_mass_['A'] = 71.037114;
  monoisotopic_mass_['C'] = 103.00919;
  monoisotopic_mass_['D'] = 115.02694;
  monoisotopic_mass_['E'] = 129.04259;
  monoisotopic_mass_['F'] = 147.06841;
  monoisotopic_mass_['G'] = 57.021464;
  monoisotopic_mass_['H'] = 137.08591;
  monoisotopic_mass_['I'] = 113.08406;
  monoisotopic_mass_['J'] = 113.08406;
  monoisotopic_mass_['K'] = 128.09496;
  monoisotopic_mass_['L'] = 113.08406;
  monoisotopic_mass_['M'] = 131.04048;
  monoisotopic_mass_['N'] = 114.04293;
  monoisotopic_mass_['P'] = 97.052764;
  monoisotopic_mass_['Q'] = 128.05858;
  monoisotopic_mass_['R'] = 156.10111;
  monoisotopic_mass_['S'] = 87.032029;
  monoisotopic_mass_['T'] = 101.04768;
  monoisotopic_mass_['V'] = 99.068414;
  monoisotopic_mass_['W'] = 186.07931;
  monoisotopic_mass_['Y'] = 163.06333;
*/
}

// Sets L and J to be equivalent to I.

void Amino::EquateLeucineWithIsoleucine() {
  SetEquivalentCharacters('L', 'I');
  SetEquivalentCharacters('J', 'I');
}

// Constructor for simple DNA alphabet.  Sets
// upper and lower case to be equivalent.  Doesn't
// handle ambiguous characters.

DNAAlphabet::DNAAlphabet() {
  CreateAlphabetFromString("actg");
  SetEquivalentCharacters('A', 'a');
  SetEquivalentCharacters('C', 'c');
  SetEquivalentCharacters('G', 'g');
  SetEquivalentCharacters('T', 't');
}

