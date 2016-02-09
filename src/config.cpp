/********************************************************/
// config.cpp
//
// SIPROS written by Doug Hyatt and Chongle Pan
//
// simple class to hold configuration information
/********************************************************/

#include "config.h"

// Constructor - Initializes to Reasonable Tryptic Cleavage Rules
// Eventually write constructor for a configuration file read
Config::Config(class PTM_List *cptm, class ProteinFragmentRules *cpfr) {
  database_path_ = "";
  working_directory_ = "";
  ptm_config_ = cptm;
  pfr_config_ = cpfr;
}

char Config::get_separator() {
  #if _WIN32
  return '\\';
  #else
  return '/';
  #endif
}

bool Config::read_config_file(string filename) {
  char line[MAX_LINE_LEN];
  ifstream config_stream(filename.c_str());
  if(config_stream.is_open() == 0) { return false; }
  while(!config_stream.eof()) {
    config_stream.getline(line, MAX_LINE_LEN);
    string sline(line);
    TokenVector words(sline, ",\r\t\n");
    if(words.size() < 2) continue;
    if(words[0] == "ptm") {
      if(words.size() != 4) { 
        cerr << "PTM Lines must contain 3 additional fields: ptm,residue,symbol,mass_offset" << endl; 
        return false; 
      }
      istringstream input(words[3]);
      double mass_shift;
      input >> mass_shift;
      cerr << "Adding PTM        :\t" << words[1][0] << "," << words[2][0] << "," << mass_shift << endl;
      if(ptm_config_->add_ptm(words[1][0], words[2][0], mass_shift) == false) {
        cerr << "Invalid PTM line, format should be: ptm,residue,symbol,mass_offset" << endl; 
        return false; 
      } 
    }
    if(words[0] == "database") { 
      database_path_ = words[1];
      cerr << "Database          :\t" << database_path_ << endl;
    }
    if(words[0] == "directory") { 
      working_directory_ = words[1];
      if(working_directory_[working_directory_.size()-1] != get_separator())
        working_directory_ = working_directory_ + get_separator();
      cerr << "Working Directory :\t" << working_directory_ << endl;
    }
    if(words[0] == "after_residues") { 
      pfr_config_->set_after_residues(words[0]);
      cerr << "After Residues    :\t" << words[0] << endl;
    }
    if(words[0] == "before_residues") { 
      pfr_config_->set_before_residues(words[0]);
      cerr << "Before Residues   :\t" << words[0] << endl;
    }
    if(words[0] == "max_missed_cleave") { 
      istringstream input(words[1]);
      int tmp_int;
      input >> tmp_int;
      pfr_config_->set_max_missed_cleave(tmp_int);
      cerr << "Max Missed Cleave :\t" << tmp_int << endl;
    }
    if(words[0] == "min_specific_cleave") { 
      istringstream input(words[1]);
      int tmp_int;
      input >> tmp_int;
      pfr_config_->set_min_specific_cleave(tmp_int);
      cerr << "Min Spec Cleave   :\t" << tmp_int << endl;
    }
    if(words[0] == "min_peptide_length") { 
      istringstream input(words[1]);
      int tmp_int;
      input >> tmp_int;
      pfr_config_->set_min_peptide_length(tmp_int);
      cerr << "Min Peptide Length:\t" << tmp_int << endl;
    }
    if(words[0] == "max_peptide_length") { 
      istringstream input(words[1]);
      int tmp_int;
      input >> tmp_int;
      pfr_config_->set_max_peptide_length(tmp_int);
      cerr << "Max Peptide Length:\t" << tmp_int << endl;
    }
    if(words[0] == "test_start_removal") { 
      istringstream input(words[1]);
      bool tmp_bool;
      input >> tmp_bool;
      pfr_config_->set_test_start_removal(tmp_bool);
      cerr << "Test Start Removal:\t" << tmp_bool << endl;
    }
    if(words[0] == "min_peptide_mass") { 
      istringstream input(words[1]);
      double tmp_dub;
      input >> tmp_dub;
      pfr_config_->set_min_peptide_mass(tmp_dub);
      cerr << "Min Peptide Mass  :\t" << tmp_dub << endl;
    }
    if(words[0] == "max_peptide_mass") { 
      istringstream input(words[1]);
      double tmp_dub;
      input >> tmp_dub;
      pfr_config_->set_max_peptide_mass(tmp_dub);
      cerr << "Max Peptide Mass  :\t" << tmp_dub << endl;
    }
    if(words[0] == "max_polymorphisms") { 
      istringstream input(words[1]);
      int tmp_int;
      input >> tmp_int;
      pfr_config_->set_max_polymorphisms(tmp_int);
      cerr << "Max Polymorphisms :\t" << tmp_int << endl;
    }
    if(words[0] == "parent_mass_error") { 
      istringstream input(words[1]);
      double tmp_dub;
      input >> tmp_dub;
      pfr_config_->set_parent_mass_error(tmp_dub);
      cerr << "Parent Mass Error :\t" << tmp_dub << endl;
    }
    if(words[0] == "fragment_mass_error") { 
      istringstream input(words[1]);
      double tmp_dub;
      input >> tmp_dub;
      pfr_config_->set_fragment_mass_error(tmp_dub);
      cerr << "Fragment Mass Err :\t" << tmp_dub << endl;
    }
  }
  return true;
}
