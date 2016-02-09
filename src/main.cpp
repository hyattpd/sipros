/********************************************************/
// main.cpp
//
// SIPROS written by Doug Hyatt and Chongle Pan
//
// SIPROS Main Program:  Reads in a configuration file
// and performs polymorphism and PTM searches using the
// specified parameters.
/********************************************************/

#include "config.h"
#include "directoryStructure.h"
#include "denovo.h"
#include "scanmass.h"
#include "ms2scan.h"
#include "scanindex.h"
#include "dbpeptide.h"
#include "outputtable.h"

using namespace::std;

int main(int argc, char **argv) {


  // Grab command line arguments
  vector<string> vsArguments;
  while(argc--) vsArguments.push_back(*argv++);

  // Parse the arguments
  string sConfigFilename = "";
  string sWorkingDirectory = "";
  unsigned int i  = 0;
  for(i = 1; i < vsArguments.size()-1; i++) {
    if(vsArguments[i] == "-w") { sWorkingDirectory = vsArguments[++i]; }
    else if (vsArguments[i] == "-c") { sConfigFilename = vsArguments[++i]; }
    else if (vsArguments[i] == "-h" || vsArguments[i] == "--help") {
      cout << "Usage: -w WorkingDirectory -c ConfigurationFile" << endl;
      cout << "Default: WorkingDirectory is the current directory; ";
      cout << "ConfigurationFile is SIPROSConfig.xml in the working directory" << endl;
      exit(0);
    }
    else { // unknown arguments
      cerr << "Unknown option " << vsArguments[i] << endl << endl; 
      exit(1);
    }
  }

  // If no -w option, the default working directory is current directory
  if(sWorkingDirectory == "") sWorkingDirectory = ".";

  // If no -c option, the default configuration file is SIPROSConfig.xml
  if(sConfigFilename == "")
    sConfigFilename = sWorkingDirectory + ProNovoConfig::getSeparator() + "SIPROSConfig.xml";

  // Load configuration file.
  if(!ProNovoConfig::setFilename(sConfigFilename)) {
    cerr << "Could not load config file " << sConfigFilename << endl << endl;
    exit(2);
  }

  bool rc;
  cerr << "SIPROS v1.0" << endl;

  // Get ptms from the config file
  PTM_List ptms;
  if(!ptms.populate_from_xml_config()) {
    cerr << "Error in parsing PTM rules from config " << endl;
    exit(3);
  }

  // Get protein cleavage rules from the config file
  ProteinFragmentRules prules;
  if(!prules.populate_from_xml_config()) {
    cerr << "Error in parsing protein cleavage rules from config " << endl;
    exit(4);
  }

  // Set working Directory and Database Path
  Config cfg(&ptms, &prules);
  cfg.set_working_directory(sWorkingDirectory + ProNovoConfig::getSeparator());
  cfg.set_database_path(ProNovoConfig::getFASTAfilename());

  // Get all FT2 files in working directory.  If the vector is empty,
  // exit with an error.
  vector<string> ft2files;
  DirectoryStructure working_dir(cfg.working_directory());
  working_dir.setPattern(".ft2");
  working_dir.getFiles(ft2files);
  working_dir.setPattern(".FT2");
  working_dir.getFiles(ft2files);
  if(ft2files.size() == 0) {
    cerr << "ERROR: no .ft2 files found in specified directory " << cfg.working_directory() << endl;
    exit(2);
  }
  for(size_t i = 0; i < ft2files.size(); i++) {
    ft2files[i] = cfg.working_directory() + ft2files[i];
    cerr << "FT2 File          :\t" << ft2files[i] << endl;
  }

  cerr << endl << endl;

  // Build an index of all scans and parent masses from the FT2 files.
  ScanMassList sm;
  rc = sm.ReadFt2File(ft2files);
  if(rc == false) { cerr << "error reading scanlist" << endl; }
  sm.SortMasses();
//  sm.Dump();
  cerr << "Read in scanmasslist... " << sm.Size() << " scans." << endl; fflush(stderr);

  // Build an index of all scans and file offset locations of those scans
  // from the FT2 files.
  ScanIndexSet sis;
  for(unsigned int i = 0; i < ft2files.size(); i++) {
    rc = sis.BuildIndex(ft2files[i]);
    if(rc == false) { cerr << "error building scanindex " << i << endl; fflush(stderr); }
  }

  // Read the FASTA database from the specified location into memory.
  FastaDB fd(cfg.database_path());
  cerr << "Read in database... " << fd.Size() << " sequences " << endl; fflush(stderr);
  ProteinFragmentRules pr;
  ProteinFragment pf(&pr, &ptms);
  DatabasePeptide *dbpep;
  DatabasePeptideList dbpl;
  OutputTable final_output;

  // Main Loop
  for(int i = 0; i < fd.Size(); i++) {
    cerr << "Processing Protein " << i+1 << " of " << fd.Size() << endl; fflush(stderr);
    pf.set_protein(fd[i]);

  // For each protein, generate all polymorphic peptides for that protein and
  // check them against our list of scans/parent masses. 
    while(pf.finished_poly() == false) {
      pf.increment_poly();
// pf.dumpinfo_poly();
      if(pf.is_valid_poly() == false) continue;
// CARBON
if(pf.has_polymorphism() == true) continue;
//      if(pf.is_full_unmutated() == false) continue;
      dbpep = pf.get_dbpeptide_poly();
// cout << "[" << dbpep->pseq() << "]" << endl;
      pair<int, int> massrange = sm.GetRangeFromMass(pf.getmass(), pr.parent_mass_error());
      if(massrange.first == -1 || massrange.second == -1) {
        delete dbpep;
        continue;
      }
      dbpl.AddPeptide(dbpep);
      for(int j = massrange.first; j <= massrange.second; j++) {
// cout << "VALIDSCAN\t" <<  (i+1) << "\t" << sm[j]->getfilenum() << "\t" << sm[j]->getscannum() << "\t" << pf.getmass() << endl;
        dbpl.AddPeptideAndScan(dbpep, sm[j]->file_number(), sm[j]->scan_number());
      }
//    cout << "finished peptide new count is " << dbpl.Size() << " and bounds in masslist were " << massrange.first << " " << massrange.second << endl; 
    }

    // For each protein, generate all peptides with PTMs and check them against our
    // list of scans/parent masses.
/*
    while(pf.finished_ptm() == false) {
      pf.increment_ptm();
// pf.dumpinfo_ptm(); fflush(stdout);
      if(pf.is_valid_ptm() == false) continue;
//      if(pf.is_full_unmutated() == false) continue;
      dbpep = pf.get_dbpeptide_ptm();
// cout << "[" << dbpep->sequence() << "]" << endl; fflush(stdout);
      pair<int, int> massrange = sm.GetRangeFromMass(pf.getmass(), pr.parent_mass_error());
      if(massrange.first == -1 || massrange.second == -1) {
        delete dbpep;
        continue;
      }
      dbpl.AddPeptide(dbpep);
      for(int j = massrange.first; j <= massrange.second; j++) {
// cout << "VALIDSCAN\t" <<  (i+1) << "\t" << sm[j]->file_number() << "\t" << sm[j]->scan_number() << "\t" << pf.getmass() << endl; fflush(stdout);
        dbpl.AddPeptideAndScan(dbpep, sm[j]->file_number(), sm[j]->scan_number());
      }
//    cout << "finished peptide new count is " << dbpl.Size() << " and bounds in masslist were " << massrange.first << " " << massrange.second << endl;  fflush(stdout);
    }
*/

    // When we hit 6,000,000 entries, we score all the potential scans and output the ones 
    // that pass both thresholds (preliminary and final).  Those that fail the preliminary
    // scoring function never make it to the final scoring function.
    if(dbpl.Size() >= 6000000) {
cerr << "before proc" << endl; fflush(stderr);
      dbpl.Process(sis, 15.5);
cerr << "after proc" << endl; fflush(stderr);
//      final_output.Update(dbpl, 15.0);
cerr << "after output update" << endl; fflush(stderr);
      dbpl.Dump(sis, 10.0);
// ship off to final output
cerr << "clearing..." << endl; fflush(stderr);
      dbpl.Clear();
cerr << "cleared!..." << endl; fflush(stderr);
    }
    cerr << "Done with Protein " << i+1 << " of " << fd.Size() << "size of dbpl is " << dbpl.Size() << endl; fflush(stderr);
  }

  // Output the final piece (if less < 6M).
  dbpl.Process(sis, 15.5);
  dbpl.Dump(sis, 10.0);
//  final_output.Update(dbpl, 15.0);
  dbpl.Clear();
}
