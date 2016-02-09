#include "ptm.h"

bool PTM_List::add_ptm(char res, char sym, double mass) {
  char check[10];
  check[0] = res; check[1] = '\0';
  string scheck(check);
  if(scheck.find_first_not_of("ACDEFGHIJKLMNPQRSTVWY[]") != string::npos) return false;
  _residue.push_back(res);
  _symbol.push_back(sym);
  _mass_shift.push_back(mass);
  return true;
}

void PTM_List::dump() {
  for(int i = 0; i < size(); i++) { cout << i << "\t" << _residue[i] << "\t" << _symbol[i] << "\t" << _mass_shift[i] << endl; }
}

bool PTM_List::populate_from_xml_config() {
  map<char, string> mPTMinfo;
  ProNovoConfig::getPTMinfo(mPTMinfo);
  map<char, string>::iterator iter;
  string tmp_sym, tmp_res_list;
  double tmp_mass;

  for(iter = mPTMinfo.begin(); iter != mPTMinfo.end(); iter++) {
    tmp_sym = iter->first;
    tmp_mass = ProNovoConfig::getResidueMass(tmp_sym);
    tmp_res_list = iter->second;
    for(size_t i = 0; i < tmp_res_list.size(); i++) {
      if(add_ptm(tmp_res_list[i], iter->first, tmp_mass) == false)
        return false;
    }
  }
  return true;
}
