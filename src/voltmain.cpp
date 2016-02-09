//============================================================================
// Name        : Volt.cpp
// Author      : CP
// Version     :
// Copyright   :
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "PeptideScorer.h"
#include "proNovoConfig.h"

using namespace std;

int main( int argc, char * argv[] )
{

	// get all arguments into vsArguments
	vector<string> vsArguments;
	while(argc--)
		vsArguments.push_back(*argv++);

	// read the options -w working directory and -c configuration file
	string sConfigFilename = "";
	string sWorkingDirectory = "";
	unsigned int i  = 0;
	while( i < vsArguments.size())
	{
		if(vsArguments[i] == "-w")
		{
			i = i + 1;
			if( i >= vsArguments.size())
				break;
			sWorkingDirectory = vsArguments[i];
		}
		else if (vsArguments[i] == "-c")
		{
			i = i + 1;
			if( i >= vsArguments.size())
				break;
			sConfigFilename = vsArguments[i];
		}
		else if (vsArguments[i] == "--help" )
		{
			cout << "Usage: -w WorkingDirectory -c ConfigurationFile" << endl;
			cout << "Default: WorkingDirectory is the current directory; ConfigurationFile is ProNovoConfig.xml in the working directory" << endl;
			return 0;
		}
		else
		{
			// unknown arguments
		}
		i = i + 1;
	}

	// if no -w option, the default working directory is current directory
	if( sWorkingDirectory == "" )
		sWorkingDirectory = ".";
	// if no -c option, the default configuration file is ProNovoConfig.xml
	if( sConfigFilename == "" )
		sConfigFilename = sWorkingDirectory + ProNovoConfig::getSeparator() + "ProNovoConfig.xml";

	// load configuration file.
	if( !ProNovoConfig::setFilename( sConfigFilename ) )
		return 0;


	// test ProNovoConfig
	vector<string> vsSingleResidueNames = ProNovoConfig::vsSingleResidueNames;
	vector<double> vdSingleResidueMasses = ProNovoConfig::vdSingleResidueMasses;

	for (unsigned int n = 0; n < vsSingleResidueNames.size(); ++n)
	{
		cout << "Residue " << vsSingleResidueNames[n] << "\t" << vdSingleResidueMasses[n] << endl;
	}
	cout << "TerminusMassN	= " << ProNovoConfig::getTerminusMassN() << endl;
	cout << "TerminusMassC	= " << ProNovoConfig::getTerminusMassC() << endl;
	cout << "FASTA DB =	" << ProNovoConfig::getFASTAfilename() << endl;

	cout << ProNovoConfig::getFragmentationMethod() << endl;
	// retrieve <MAX_PTM_COUNT> and <MAX_POLYMORPHISM_COUNT>
	cout << ProNovoConfig::getMaxPTMcount() << endl;
	cout << ProNovoConfig::getMaxPolymorphismCount() << endl;

	// retrieve <PEPTIDE_MASS_RANGE>
	cout << ProNovoConfig::getMaxPeptideMass() << endl;
	cout << ProNovoConfig::getMinPeptideMass() << endl;

    vector<int> viParentMassWindows = ProNovoConfig::getParentMassWindows();
    for (unsigned int a = 0; a < viParentMassWindows.size(); ++a) {
		cout << "Window = " << viParentMassWindows[a] << endl;
	}

	// retrieve <CLEAVAGE_RULES>
	cout << ProNovoConfig::getCleavageAfterResidues() << endl;
	cout << ProNovoConfig::getCleavageBeforeResidues() << endl;
	cout << ProNovoConfig::getMaxMissedCleavages() << endl;
	cout << ProNovoConfig::getMinSpecificCleavages() << endl;
	cout << boolalpha << ProNovoConfig::getTestStartRemoval() << endl;

	map<char, string> mPTMinfo;
	ProNovoConfig::getPTMinfo(mPTMinfo);
    map<char, string>::iterator iter;


    for(iter = mPTMinfo.begin(); iter != mPTMinfo.end(); iter++)
    {
      cout << iter->first << " -> " << iter->second << endl;
    }

/*
	// test peptideScore
	PeptideScorer mainPeptideScorer;
	vector<double> vdMZ;
	vector<double> vdIntensity;
	vector<int> viChargeState;


	ifstream in("E:\\Research\\DeNovo\\workspace\\Volt\\MS838_Cycle8.txt");
//	ifstream in("D:\\Sipros\\Volt\\MS838_Cycle8.txt");
	string sLine;
	double dMZ;
	double dInt;
	double dResolution;
	double dNoise;
	double dSignal;
	int iChargeState;

	while(getline(in, sLine))
	{
		istringstream issStream(sLine);
		issStream >> dMZ >> dInt >> dResolution >> dNoise >> dSignal  >> iChargeState;
		vdMZ.push_back(dMZ);
		vdIntensity.push_back(dInt);
		viChargeState.push_back(iChargeState);
	}

	mainPeptideScorer.setMS2(vdMZ, vdIntensity, viChargeState,
			525.619265, 3);


	vector<string> vsPeptide;
	vector<double> vdScore;
	vsPeptide.push_back("[YRPPAESAASGITVR]");
//	vsPeptide.push_back("GENSTKSWAAGIEAR");
//	vsPeptide.push_back("VYPAGVSVVRRSQR");
//	vsPeptide.push_back("TAEIEKKAAEKKAR");
//	vsPeptide.push_back("KTKSSTAAQPMAAER");
//	vsPeptide.push_back("ATLAEAADISSENQR");
//	vsPeptide.push_back("QLAFVLEASQLAER");
//	vsPeptide.push_back("SITSKGAGSPGVIAVTK");
//	vsPeptide.push_back("GFISPMPFPSSSVVR");


	cout << "Peptide number = " << vsPeptide.size() << endl;
	for( unsigned int i = 0; i < vsPeptide.size(); ++i )
	{

		cout << vsPeptide[i]
//		<< "\t\t" << mainPeptideScorer.prelimaryScorePeptide( vsPeptide[i])
		<< "\t\t" << mainPeptideScorer.prelimaryScorePeptide2( vsPeptide[i])
//		<< "\t\t" << mainPeptideScorer.prelimaryScorePeptide3( vsPeptide[i])
		<< "\t\t" << mainPeptideScorer.primaryScorePetide( vsPeptide[i])
//		<< "\t\t" << mainPeptideScorer.primaryScorePetide2( vsPeptide[i])
//		<< "\t\t" << mainPeptideScorer.primaryScorePetide3( vsPeptide[i])
//		<< "\t\t" << mainPeptideScorer.primaryScorePetide4( vsPeptide[i])
		<< endl;
	}
	time_t time0, time1;
	double dTimeElapsed;
	time (& time0);


	for (unsigned int k = 0; k < 400; ++k)
	{
		for( unsigned int i = 1; i < vsPeptide.size(); ++i )
		{
//			mainPeptideScorer.prelimaryScorePeptide( vsPeptide[i]);
//			mainPeptideScorer.prelimaryScorePeptide2( vsPeptide[i]);
			mainPeptideScorer.prelimaryScorePeptide3( vsPeptide[i]);
//			mainPeptideScorer.primaryScorePetide( vsPeptide[i]);
//			mainPeptideScorer.primaryScorePetide2( vsPeptide[i]);
//			mainPeptideScorer.primaryScorePetide3( vsPeptide[i]);
//			mainPeptideScorer.primaryScorePetide4( vsPeptide[i]);
		}
	}


	time (& time1);
	dTimeElapsed = difftime (time1,time0);
	cout << "Time elapsed = " << difftime (time1,time0)
	<< " sec" << endl << endl;


	// test ProNovoConfig::configIsotopologue.computeProductIon
	vector< vector<double> > vvdYionMass;
	vector< vector<double> > vvdYionProb;
	vector< vector<double> > vvdBionMass;
	vector< vector<double> > vvdBionProb;
	ProNovoConfig::configIsotopologue.computeProductIon(
			"YRPPAESAASGITVR",
			vvdYionMass, vvdYionProb,
			vvdBionMass, vvdBionProb);
	unsigned int n;
	unsigned int m;
	cout << "Y ion" << endl;
	for (n = 0; n < vvdYionMass.size(); ++n)
	{
		cout << n+1 << '\t';
		for (m = 0; m < vvdYionMass[n].size(); ++m)
		{
			cout << vvdYionMass[n][m] << "-" << vvdYionProb[n][m] << '\t';
		}
		cout << endl;
	}

	cout << "B ion" << endl;
	for (n = 0; n < vvdBionMass.size(); ++n)
	{
		cout << n+1 << '\t';
		for (m = 0; m < vvdBionMass[n].size(); ++m)
		{
			cout << vvdBionMass[n][m] << "-" << vvdBionProb[n][m] << '\t';
		}
		cout << endl;
	}



	vector<double> vdMass;
	vector<double> vdProb;
	double dMono = 0;
	ProNovoConfig::configIsotopologue.getPolyAveragine(3073.59, dMono, vdMass, vdProb);
	cout << "dMono = " << dMono << endl;
	for (unsigned int i = 0; i < vdMass.size(); ++i)
	{
		cout << vdMass[i] << "	=	"<< vdProb[i] << endl;
	}

*/

	cout << "Done!" << endl;

	return 0;
}
