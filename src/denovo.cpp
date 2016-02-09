/********************************************************/
// denovo.h
//
// SIPROS written by Doug Hyatt and Chongle Pan
//
// Class to process DeNovo Sequence tags, search them
// against a database (allowing for polymorphisms), and
// return dbpeptidescans (see dbpeptide.h) for final scoring.
/********************************************************/

#include "denovo.h"

DeNovoTag::DeNovoTag() {
  file_name_ = ""; scan_number_ = 0; parent_mz_ = 0.0;
  charge_ = 0; index_ = ""; left_flanking_mass_ = 0.0; 
  right_flanking_mass_ = 0.0; sequence_tag_ = ""; score_ = 0.0;
  set_size_ = 0;
}

// Opens a denovo tag file from Vonode.  Returns true if
// file is opened successfully, false otherwise.
bool DeNovoTag::OpenFile(string fn) {
  infstream_.open(fn.c_str(), ifstream::in);
  if(infstream_.is_open() == 0) return false;
  return true;
}

// Reads the next entry from a Vonode file and stores
// all the information in the DeNovo class variables.

int DeNovoTag::ReadNextEntry() {
  string line;
  if(infstream_.is_open() == 0) return -1;
  if(infstream_.eof()) return 1;
  if(!getline(infstream_, line)) return 1;
  TokenVector words(line, "\r\t\n");
//  if(words.size() != 10) return -1;  temporary for parsing other files
  if(words.size() < 10) return -1;
  vector<string>::iterator iter;
  iter = words.begin();
  file_name_ = *iter;
  file_name_ += ".FT2";
  istringstream input;
  input.str(*(++iter));
  input >> scan_number_;
  input.clear();
  input.str(*(++iter));
  input >> parent_mz_;
  input.clear();
  input.str(*(++iter));
  input >> charge_;
  input.clear();
  index_ = *(++iter);
  input.str(*(++iter));
  input >> left_flanking_mass_;
  input.clear();
  sequence_tag_ = *(++iter);
  input.str(*(++iter));
  input >> right_flanking_mass_;
  input.clear();
  input.str(*(++iter));
  input >> score_;
  input.clear();
  input.str(*(++iter));
  input >> set_size_;
  input.clear();
/*
cout << "Filename: " << file_name_ << endl;
cout << "ScanNum: " << scan_number_ << endl;
cout << "ParentMZ: " << parent_mz_ << endl;
cout << "ChargeState: " << charge_ << endl;
cout << "Index: " << index_ << endl;
cout << "LeftMass: " << left_flanking_mass_ << endl;
cout << "RightMass: " << right_flanking_mass_ << endl;
cout << "SeqTag: " << sequence_tag_ << endl;
cout << "Score: " << score_ << endl;
cout << "SetSize: " << set_size_ << endl;
*/
  return 0;
}

// Flips a Vonode tag sequence and its flanking masses.

void DeNovoTag::FlipTag() {
  double tmp;
  tmp = right_flanking_mass_;
  right_flanking_mass_ = left_flanking_mass_;
  left_flanking_mass_ = tmp;
  char tmp_char;
  for(int i = 0; i < (int)(sequence_tag_.length()/2); i++) {
    tmp_char = sequence_tag_[sequence_tag_.length() - i - 1];
    sequence_tag_[sequence_tag_.length() - i - 1] = sequence_tag_[i];
    sequence_tag_[i] = tmp_char;
  }
}

// Takes a FASTA database and processes the current Vonode tag
// by searching it against all sequences in the database according
// to the cleavage/polymorphism/ptm rules specified in the Protein
// Fragment Rules class.

void DeNovoTag::Process(FastaDB& fd, ProteinFragmentRules& pr, Amino *a, DatabasePeptideList &dbpl, ScanIndexSet &sindexset) {
  int i, j, k, l, m, n, o, cl, mismatch = 0, lctr, rctr;
  int lmismatch[100], rmismatch[100], lclv, rclv;
  int cl_lim = 2 - pr.min_specific_cleave(), missed_cleavage;
  string ltmp[100], rtmp[100], peptide;
  string res = "ACDEFGHIKLMNPQRSTVWY";
  Blosum62 blos62_matrix;
  DatabasePeptide *dbp;
  if(sequence_tag_ == "SequenceTag") return;
  Pattern denovo_pat(sequence_tag_, a);
  denovo_pat.set_id(index_); 
  Bitap b;
  double lmass, rmass;
  mismatch = pr.max_polymorphisms();
  if(denovo_pat.Size() < 3) return;
  else if(denovo_pat.Size() == 3 && mismatch > 1) mismatch = 1;
  else mismatch = 2;
  for(i = 0; i < fd.Size(); i++) { b.KMismatch(fd[i], &denovo_pat, *a, mismatch); }
  for(i = 0; i < b.Size(); i++) {
    mismatch = b[i]->mismatch_count();
    lmass = 0.0;
    for(j = b[i]->position() - 1; j >= 0; j--) {
      lctr = 0; 
      if(j != b[i]->position() - 1) lmass += a->MonoisotopicMass(b[i]->Residue(j));
      if(lmass > left_flanking_mass_ + 200.0) break;
      if(lmass < left_flanking_mass_ - 200.0) continue;
      if(fabs(lmass+pr.left_offset() - left_flanking_mass_) <= pr.fragment_mass_error()) { 
        ltmp[lctr] = b[i]->TextString().substr(j, b[i]->position() - 1 - j);
        lmismatch[lctr] = 0;
        lctr++;
      }
      if(mismatch < pr.max_polymorphisms()) {
        for(k = j; k < b[i]->position() - 1; k++) {
          for(cl = 0; cl < 20; cl++) {
            if(res[cl] == b[i]->TextString()[k]) continue;
            if(fabs(lmass+pr.left_offset() - a->MonoisotopicMass(b[i]->TextString()[k]) + a->MonoisotopicMass(res[cl]) - left_flanking_mass_) < pr.fragment_mass_error()) {
              lmismatch[lctr] = 1; 
              ltmp[lctr] = b[i]->TextString().substr(j, b[i]->position() - 1 - j);
              ltmp[lctr][k-j] = res[cl];
              lctr++;
            }
          }
        } 
      }
      for(m = 0; m < lctr; m++) {
        peptide = ltmp[m] + b[i]->PatternString();
        if(j != 0 && (j != 1 || pr.test_start_removal() == false) && pr.is_cleavage_site(b[i]->TextString()[j-1], peptide[0]) == 0) lclv = 1;
        else lclv = 0;
        if(lclv > cl_lim) continue;
        mismatch = b[i]->mismatch_count() + lmismatch[m];
        rmass = 0.0; 
        for(k = b[i]->position() - 2 + b[i]->MatchString().length(); k < (int)b[i]->TextString().length(); k++) {
          rctr = 0;
          if(k != (int)(b[i]->position() - 2 + b[i]->MatchString().length())) rmass += a->MonoisotopicMass(b[i]->Residue(k));
          if(rmass > right_flanking_mass_ + 200.0) break;
          if(rmass < right_flanking_mass_ - 200.0) continue;
          if(fabs(rmass+pr.right_offset() - right_flanking_mass_) <= pr.fragment_mass_error()) {
            rtmp[rctr] = b[i]->TextString().substr(b[i]->position()-1+b[i]->MatchString().length(), k-b[i]->position()-b[i]->MatchString().length()+2);
            rmismatch[rctr] = 0;
            rctr++;
          }
          if(mismatch < pr.max_polymorphisms()) {
            for(l = k; l > (int)(b[i]->position()-2+b[i]->MatchString().length()); l--) {
              for(cl = 0; cl < 20; cl++) {
                if(res[cl] == b[i]->TextString()[l]) continue;
                if(fabs(rmass+pr.right_offset() - a->MonoisotopicMass(b[i]->TextString()[l]) + a->MonoisotopicMass(res[cl]) - right_flanking_mass_) < pr.fragment_mass_error()) {
                  rmismatch[rctr] = 1; 
                  rtmp[rctr] = b[i]->TextString().substr(b[i]->position()-1+b[i]->MatchString().length(), k-b[i]->position()-b[i]->MatchString().length()+2);
                  rtmp[rctr][l-b[i]->position()-b[i]->MatchString().length()+1] = res[cl];
                  rctr++;
                }
              }
            } 
          }
          if(rctr > 0) {
            for(n = 0; n < rctr; n++) {
              peptide += rtmp[n];
              if(k != (int)(b[i]->TextString().length()-1) && pr.is_cleavage_site(peptide[peptide.length()-1], b[i]->TextString()[k+1]) == 0) rclv = 1;
              else rclv = 0;
              if(lclv+rclv > cl_lim) continue;
              // Check for missed cleavage sites
              missed_cleavage = 0;
              for(o = 0; o < (int)(peptide.length()-1); o++)
                if(pr.is_cleavage_site(peptide[o], peptide[o+1]) == 1) missed_cleavage++;
              if(missed_cleavage > pr.max_missed_cleave()) continue;
              cout << "DB: " << b[i]->TextString().substr(j, k-j+1) << "\t" << "Tag: " << ltmp[m] << b[i]->PatternString() << rtmp[n] << "\t" << lmismatch[m] << " " << b[i]->mismatch_count() << " " << rmismatch[n] << "\t";
              string dbpseq = ltmp[m] + b[i]->PatternString() + rtmp[n];
              if(sindexset.FileNumber(file_name_) == -1 && sindexset.BuildIndex(file_name_) == false) {
                cerr << "Failed to build index file " << file_name_ << endl;  
                exit(11);
              }
              dbp = new DatabasePeptide(b[i]->text_pointer(), b[i]->position() - ltmp[m].length(), dbpseq.length(), dbpseq, 0);
              dbpl.AddPeptide(dbp);
              dbpl.AddPeptideAndScan(dbp, sindexset.FileNumber(file_name_), scan_number_);
              blos62_matrix.ScorePair(denovo_pat.word(), b[i]->MatchString());
              b[i]->Print(1); 
              cout << "\t" << blos62_matrix.bitscore() << "\t" << blos62_matrix.positives() << endl;
            }
          }
        }
      }
    }
  }
}
