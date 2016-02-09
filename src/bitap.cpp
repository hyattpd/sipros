/********************************************************/
// bitap.cpp
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

#include "bitap.h"

Bitap::Bitap() {}

// K Mismatch Bitap Algorithm lifted directly from the
// Wikipedia page on Bitap and modified slightly.
// Locates all instances of a pattern p in a text t
// containing at most k mismatches.

void Bitap::KMismatch(Text *t, Pattern *p, Alphabet &a, int k) {
  int i, j, seen;
  unsigned long tmp, old_R;
  unsigned long *R = new unsigned long [k+1];
  Match *hsp;

  for(i = 0; i <= k; i++) R[i] = ~1;

  for(i = 0; i < (int)(t->word().size()); i++) {
    old_R = R[0];
    R[0] |= p->BitMask(a.LetterToNumber(t->word()[i]));
    R[0] <<= 1;
    for(j = 1; j <= k; j++) {
      tmp = R[j];
      R[j] = (old_R & (R[j] | p->BitMask(a.LetterToNumber(t->word()[i])))) << 1;
      old_R = tmp;
    }

    seen = 0;
    for(j = 0; j <= k; j++) {
      if((seen == 0) && (0 == (R[j] & (1UL << p->Size())))) {
        hsp = new Match;
        hsp->set_mismatch_count(j);
        hsp->set_position(i + 2 - p->Size());
        hsp->set_pattern_pointer(p);
        hsp->set_text_pointer(t);
        matches_.push_back(hsp);
        seen = 1;
      }
    }
  } 

  delete R;
}

// Prints out each match in the bitap class.
// The flag just determines if we print and endline
// after each hit or not (0 = yes, 1 = no).

void Bitap::Print(bool flag) {
  for(size_t i = 0; i < matches_.size(); i++) matches_[i]->Print(flag);
}

// This returns the match string from within the text.
// i.e. complete with any mismatches if there were any.

string Match::MatchString() {
  return text_pointer_->word().substr(position_-1, pattern_pointer_->Size());
}

// Prints information on the match.

void Match::Print(bool flag) {
  cout << pattern_pointer_->id() << "\t";
  cout << pattern_pointer_->word() << "\t";
  if(text_pointer_->id().length() > 40) cout << ">" << text_pointer_->id().substr(0, 37) << "..." << "\t";
  else cout << ">" << text_pointer_->id() << "\t";
  cout << MatchString() << "\t" << position_ << "\t" << mismatch_count_;
  if(flag == 0) cout << endl;
}
