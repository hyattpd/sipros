/********************************************************/
// profrag.cpp
//
// SIPROS written by Doug Hyatt and Chongle Pan
//
// Class to generate tryptic peptides according to
// provided cleavage rules.  Also generates loops for
// single amino acid polymorphisms and single post-
// translational modifications.
/********************************************************/

#include "profrag.h"

// Constructor - Initializes to Reasonable Tryptic Cleavage Rules
// Eventually write constructor for a configuration file read
ProteinFragmentRules::ProteinFragmentRules() {
  after_residues_ = "KR";
  before_residues_ = "ACDEFGHIJKLMNPQRSTVWXY[]";
  max_missed_cleave_ = 2;
  min_specific_cleave_ = 2;
  min_peptide_length_ = 6;
  max_peptide_length_ = 60;
  test_start_removal_ = true;
  min_peptide_mass_ = 600.0;
  max_peptide_mass_ = 6000.0;
  max_polymorphisms_ = 2;
  parent_mass_error_ = 0.04;
  fragment_mass_error_ = 0.04;
  left_offset_ = 0.0;
  right_offset_ = 18.01505;
}

// Just returns true or false if the cleavage site has valid
// characters before and after the site.  After typically must
// be a K or an R.  True = a specific cleavage site.
bool ProteinFragmentRules::is_cleavage_site(char c1, char c2) {
  if(after_residues_.find(c1, 0) != string::npos &&
     before_residues_.find(c2, 0) != string::npos)
    return true;
  else return false;
}

// Constructor: sets protein, protein fragment rules, and the ptm list.
ProteinFragment::ProteinFragment(ProteinFragmentRules *pfr, PTM_List *pl) {
  protein = "";
  pr = pfr;
  ptmlist = pl;
}

void ProteinFragment::set_protein(class Text *pro) {
  dbid = pro;
  protein = pro->word();
  protein_with_termini = "["+protein+"]";
cout << "1: " << protein << endl;
cout << "2: " << protein_with_termini << endl;
  residues = "ACDEFGHIKLMNPQRSTVWY*";
  cl_map.clear();
  for(int i = 0; i < (int)protein.length()-1; i++) {
    if(pr->is_cleavage_site(protein[i], protein[i+1]) == 1)
      cl_map.push_back(1);
    else cl_map.push_back(0);
  } 
  cl_map.push_back(1);
  mass_map.clear();
  msum_map.clear();
  double msum = 0.0;
  for(int i = 0; i < (int)protein.length(); i++) {
    msum += aa.MonoisotopicMass(protein[i]);
    msum_map.push_back(msum);
    mass_map.push_back(aa.MonoisotopicMass(protein[i]));
  } 
  mutpos = -1;
  mutval = -1;
  begin = -1;
  end = -1;
  peptide = "";
  mass = 0.0;
  left_cleave = 0;
  right_cleave = 0;
  ptm_pos = -1;
  ptm_val = -1;
  ptm_loop_begin = -1;
  ptm_loop_end = -1;
  done_poly = false;
  done_ptm = false;
// for(int i = 0; i < (int)protein.length(); i++) { cout << "[" << i << " " << cl_map[i] << "]"; } cout << endl;
}

void ProteinFragment::increment_poly() {
  int clcheck = 0, hit_end = 0;
  if(end != -1) {
    end++;
    if(end < (int)protein.length()) { 
      mass = 18.01505 + msum_map[end];
      if(begin > 0) mass -= msum_map[begin-1];
      mass -= mass_map[mutpos];
      mass += aa.MonoisotopicMass(residues[mutval]);
      if(end > mutpos && (end == (int)(protein.length()-1) || cl_map[end] == true)) right_cleave++;
      else if(left_cleave+mut_cleave+right_cleave-1 >= pr->max_missed_cleave()) hit_end = 1;
      if(hit_end == 0) {
        clcheck = (left_cleave + mut_cleave + right_cleave-1);
        if(clcheck <= pr->max_missed_cleave()) return;
      }
    }
    end = -1;
  }
  if(begin != -1) {
    right_cleave = 0;
    begin--;
    end = mutpos;
    if(begin >= 0) {  
      mass = 18.01505 + msum_map[end];
      if(begin > 0) mass -= msum_map[begin-1];
      mass -= mass_map[mutpos];
      mass += aa.MonoisotopicMass(residues[mutval]);
      if(begin < mutpos && cl_map[begin] == true) left_cleave++;
      clcheck = (left_cleave + mut_cleave + right_cleave-1);
      if(clcheck <= pr->max_missed_cleave()) return;
    }
    end = -1;
    begin = -1;
  }
  if(mutval != -1) {
    mutval++;
    begin = mutpos;
    end = mutpos; 
    left_cleave = 0;
    right_cleave = 0;
    if(mutval < 20) { 
      mass = 18.01505 + msum_map[end];
      if(begin > 0) mass -= msum_map[begin-1];
      mass -= mass_map[mutpos];
      mass += aa.MonoisotopicMass(residues[mutval]);
      if(mutpos != (int)(protein.length())-1) {
        mut_cleave = pr->is_cleavage_site(residues[mutval], protein[mutpos+1]);
        mut_cleave -= cl_map[mutpos];
      }
      else mut_cleave = 1;
      if(begin < mutpos && cl_map[begin] == true) left_cleave++;
      if(end > mutpos && (end == (int)(protein.length()-1) || cl_map[end] == true)) right_cleave++;
      clcheck = (left_cleave + mut_cleave + right_cleave-1);
      if(clcheck <= pr->max_missed_cleave()) return;
    }
    mutval = -1;
    begin = -1;
    end = -1;
  }
  if(mutpos != -1) {
    mutpos++;
    begin = mutpos;
    end = mutpos;
    left_cleave = 0;
    right_cleave = 0;
    mutval = 0;
    if(mutpos < (int)protein.length()) {
      mass = 18.01505 + msum_map[end];
      if(begin > 0) mass -= msum_map[begin-1];
      mass -= mass_map[mutpos];
      mass += aa.MonoisotopicMass(residues[mutval]);
      if(mutpos != (int)(protein.length())-1) {
        mut_cleave = pr->is_cleavage_site(residues[mutval], protein[mutpos+1]);
        mut_cleave -= cl_map[mutpos];
      }
      else mut_cleave = 1;
      if(begin < mutpos && cl_map[begin] == true) left_cleave++;
      if(end > mutpos && (end == (int)(protein.length()-1) || cl_map[end] == true)) right_cleave++;
      clcheck = (left_cleave + mut_cleave + right_cleave-1);
      if(clcheck <= pr->max_missed_cleave()) return;
    }
    done_poly = true;
    mutpos = -1;
    mutval = -1;
    begin = -1;
    end = -1;
  }
  mutpos = 0;
  begin = mutpos;
  end = mutpos;
  mutval = 0;
  mass = 18.01505 + msum_map[end];
  if(begin > 0) mass -= msum_map[begin-1];
  mass -= mass_map[mutpos];
  mass += aa.MonoisotopicMass(residues[mutval]);
  if(mutpos != (int)(protein.length())-1) {
    mut_cleave = pr->is_cleavage_site(residues[mutval], protein[mutpos+1]);
    mut_cleave -= cl_map[mutpos];
  }
  else mut_cleave = 1;
  if(begin < mutpos && cl_map[begin] == true) left_cleave++;
  if(end > mutpos && (end == (int)(protein.length()-1) || cl_map[end] == true)) right_cleave++;
  clcheck = (left_cleave + mut_cleave + right_cleave-1);
}

void ProteinFragment::increment_ptm() {
  int clcheck = 0, hit_ptm_loop_end = 0;
  if(ptm_loop_end != -1) {
    ptm_loop_end++;
    if(ptm_loop_end < (int)protein_with_termini.length()) { 
      mass = 18.01505 + msum_map[ptm_loop_end];
      if(ptm_loop_begin > 0) mass -= msum_map[ptm_loop_begin-1];
      mass += ptmlist->mass_shift(ptm_val);
      if(ptm_loop_end > ptm_pos && (ptm_loop_end == (int)(protein_with_termini.length()-1) || cl_map[ptm_loop_end] == true)) right_cleave++;
      else if(left_cleave+right_cleave-1 >= pr->max_missed_cleave()) hit_ptm_loop_end = 1;
      if(hit_ptm_loop_end == 0) {
        clcheck = (left_cleave + right_cleave-1);
        if(clcheck <= pr->max_missed_cleave()) return;
      }
    }
    ptm_loop_end = -1;
  }
  if(ptm_loop_begin != -1) {
    right_cleave = 0;
    ptm_loop_begin--;
    ptm_loop_end = ptm_pos;
    if(ptm_loop_begin >= 0) {  
      mass = 18.01505 + msum_map[ptm_loop_end];
      if(ptm_loop_begin > 0) mass -= msum_map[ptm_loop_begin-1];
      mass += ptmlist->mass_shift(ptm_val);
      if(ptm_loop_begin < ptm_pos && cl_map[ptm_loop_begin] == true) left_cleave++;
      clcheck = (left_cleave + right_cleave-1);
      if(clcheck <= pr->max_missed_cleave()) return;
    }
    ptm_loop_end = -1;
    ptm_loop_begin = -1;
  }
  if(ptm_val != -1) {
    ptm_val++;
    ptm_loop_begin = ptm_pos;
    ptm_loop_end = ptm_pos; 
    left_cleave = cl_map[ptm_pos];
    right_cleave = 0;
    while(ptmlist->residue(ptm_val) != protein_with_termini[ptm_pos] && ptm_val < ptmlist->size())
      ptm_val++;
    if(ptm_val < ptmlist->size()) { 
      mass = 18.01505 + msum_map[ptm_loop_end];
      if(ptm_loop_begin > 0) mass -= msum_map[ptm_loop_begin-1];
      mass += ptmlist->mass_shift(ptm_val);
      if(ptm_loop_begin < ptm_pos && cl_map[ptm_loop_begin] == true) left_cleave++;
      if(ptm_loop_end > ptm_pos && (ptm_loop_end == (int)(protein_with_termini.length()-1) || cl_map[ptm_loop_end] == true)) right_cleave++;
      clcheck = (left_cleave + right_cleave-1);
      if(clcheck <= pr->max_missed_cleave()) return;
    }
    ptm_val = -1;
    ptm_loop_begin = -1;
    ptm_loop_end = -1;
  }
  if(ptm_pos != -1) {
    ptm_pos++;
    ptm_loop_begin = ptm_pos;
    ptm_loop_end = ptm_pos;
    left_cleave = cl_map[ptm_pos];
    right_cleave = 0;
    ptm_val = 0;
    if(ptm_pos < (int)protein_with_termini.length()) {
      mass = 18.01505 + msum_map[ptm_loop_end];
      if(ptm_loop_begin > 0) mass -= msum_map[ptm_loop_begin-1];
      mass += ptmlist->mass_shift(ptm_val);
      if(ptm_loop_begin < ptm_pos && cl_map[ptm_loop_begin] == true) left_cleave++;
      if(ptm_loop_end > ptm_pos && (ptm_loop_end == (int)(protein_with_termini.length()-1) || cl_map[ptm_loop_end] == true)) right_cleave++;
      clcheck = (left_cleave + right_cleave-1);
      if(clcheck <= pr->max_missed_cleave()) return;
    }
    done_ptm = true;
    ptm_pos = -1;
    ptm_val = -1;
    ptm_loop_begin = -1;
    ptm_loop_end = -1;
  }
  ptm_pos = 0;
  ptm_loop_begin = ptm_pos;
  ptm_loop_end = ptm_pos;
  ptm_val = 0;
  mass = 18.01505 + msum_map[ptm_loop_end];
  if(ptm_loop_begin > 0) mass -= msum_map[ptm_loop_begin-1];
  mass += ptmlist->mass_shift(ptm_val);
  if(ptm_loop_begin < ptm_pos && cl_map[ptm_loop_begin] == true) left_cleave++;
  if(ptm_loop_end > ptm_pos && (ptm_loop_end == (int)(protein_with_termini.length()-1) || cl_map[ptm_loop_end] == true)) right_cleave++;
  clcheck = (left_cleave + right_cleave-1);
}

DatabasePeptide * ProteinFragment::get_dbpeptide_poly() {
  peptide = protein.substr(begin, end-begin+1);
  peptide[mutpos-begin] = residues[mutval];
  return (new DatabasePeptide(dbid, begin+1, peptide, 0));
}

DatabasePeptide * ProteinFragment::get_dbpeptide_ptm() {
  peptide = protein.substr(ptm_loop_begin, ptm_loop_end-ptm_loop_begin+1);
  peptide.insert(ptm_pos - ptm_loop_begin + 1, 1, ptmlist->symbol(ptm_val));
  return (new DatabasePeptide(dbid, ptm_loop_begin+1, peptide, 1));
}

bool ProteinFragment::is_valid_poly() {
  if(end - begin + 1 < pr->min_peptide_length()) return false;
  if(begin != 0 && ((mutpos == begin && pr->is_cleavage_site(protein[mutpos-1], residues[mutval]) == false)
     || (mutpos != begin && cl_map[begin-1] == false))) return false;
  if(end != (int)protein.length()-1 && ((mutpos == end && pr->is_cleavage_site(residues[mutval], 
     protein[mutpos+1]) == false) || (mutpos != end && cl_map[end] == false))) return false;
  if(residues[mutval] == protein[mutpos] && mutpos != begin) return false;
  if(residues[mutval] == 'I' && protein[mutpos] == 'L') return false;
  if(residues[mutval] == 'L' && protein[mutpos] == 'I') return false;
  return true;
}

bool ProteinFragment::is_valid_ptm() {
  if(ptmlist->residue(ptm_val) != protein[ptm_pos]) return false;
  if(ptm_loop_end - ptm_loop_begin + 1 < pr->min_peptide_length()) return false;
  if(ptm_loop_begin != 0 && ((ptm_pos == ptm_loop_begin && pr->is_cleavage_site(protein[ptm_pos-1], protein[ptm_pos]) == false)
     || (ptm_pos != ptm_loop_begin && cl_map[ptm_loop_begin-1] == false))) return false;
  if(ptm_loop_end != (int)protein.length()-1 && ((ptm_pos == ptm_loop_end && pr->is_cleavage_site(protein[ptm_pos], 
     protein[ptm_pos+1]) == false) || (ptm_pos != ptm_loop_end && cl_map[ptm_loop_end] == false))) return false;
  return true;
}

bool ProteinFragment::is_full_unmutated() {
  if(residues[mutval] != protein[mutpos]) return false;
  if(begin != 0) return false;
  if(end != (int)protein.length()-1) return false;
  return true;
}

void ProteinFragment::dumpinfo_poly() {
  string s = protein.substr(begin, end-begin+1);
  s[mutpos-begin] = residues[mutval];
  cout << "DUMP\t" << begin << "\t" << end << "\t" << mutpos << "\t" << residues[mutval] << "\t" << mass << "\t" << s << "\t" << left_cleave << "\t" << mut_cleave << "\t" << right_cleave << "\t" << is_valid_poly() << endl;
}

void ProteinFragment::dumpinfo_ptm() {
  string s = protein.substr(ptm_loop_begin, ptm_loop_end-ptm_loop_begin+1);
  s.insert(ptm_pos - ptm_loop_begin + 1, 1, ptmlist->symbol(ptm_val));
  cout << "DUMP\t" << ptm_loop_begin << "\t" << ptm_loop_end << "\t" << ptm_pos << "\t" << ptmlist->symbol(ptm_val) << "\t" << mass << "\t" << s << "\t" << left_cleave << "\t" << right_cleave << "\t" << is_valid_ptm() << endl;
}
