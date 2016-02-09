/********************************************************/
// config.h
//
// SIPROS written by Doug Hyatt and Chongle Pan
//
// Simple class to hold config file information
/********************************************************/

#ifndef _CONFIG_H
#define _CONFIG_H

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "tokenvector.h"
#include "profrag.h"
#include "ptm.h"

#define MAX_LINE_LEN 10000

using namespace std;

class Config {

  string database_path_;              // Path to the database we will be searching
  string working_directory_;          // Working directory containing FT2 files

  class PTM_List *ptm_config_;               // Config for PTMs
  class ProteinFragmentRules *pfr_config_;   // Config for Cleavage Rules

public:

  Config(class PTM_List *, class ProteinFragmentRules *);
  ~Config() {}

  string database_path() { return database_path_; }  
  string working_directory() { return working_directory_; }
  class PTM_List *ptm_config() { return ptm_config_; }  
  class ProteinFragmentRules *pfr_config() { return pfr_config_; }  

  void set_database_path(string s) { database_path_ = s; }  
  void set_working_directory(string s) { working_directory_ = s; }
  void set_ptm_config(class PTM_List *p) { ptm_config_ = p; }  
  void set_pfr_config(class ProteinFragmentRules *p) { pfr_config_ = p; }  

  bool read_config_file(string);

  char get_separator();
};

#endif
