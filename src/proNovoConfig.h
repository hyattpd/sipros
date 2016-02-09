#ifndef PRONOVOCONFIG_H_
#define PRONOVOCONFIG_H_

#include "tinyxml.h"
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "isotopologue.h"

using namespace std;

class Isotopologue;

class ProNovoConfig
{
	public:
		/*
		 * Sets up sessionwide configuration for ProRata
		 * the configurations are loaded in to memory as static variables
		 */

		static bool setFilename( const string & sConfigFileName );

		static bool setWorkingDirectory( const string & sDirectoryName );

		static string getWorkingDirectory()
		{ return sWorkingDirectory; }

		static char getSeparator();

		/*
		 * get the version number of ProNovo
		 */

		static string getProNovoVersion()
		{ return "1.0"; }

		// retrieve <INPUT_MS_FILE_TYPE>
		static string getMSfileType()
		{ return sMSFileType; }

		// retrieve <OUTPUT_ID_FILE_TYPE>
		static string getIDfileType()
		{ return sIDFileType; }

		// retrieve <FASTA_DATABASE>
		static string getFASTAfilename()
		{ return sFASTAFilename; }

		// retrieve <FRAGMENTATION_METHOD>
		static string getFragmentationMethod()
		{ return sFragmentationMethod; }

		// retrieve <MASS_ACCURACY>	<PARENT_ION>
		static double getMassAccuracyParentIon()
		{ return dMassAccuracyParentIon; }

		// retrieve <MASS_ACCURACY>	<FRAGMENT_IONS>
		static double getMassAccuracyFragmentIon()
		{ return dMassAccuracyFragmentIon; }

		static vector<int> getParentMassWindows()
		{ return viParentMassWindows;}

		// retrieve <MAX_PTM_COUNT> and <MAX_POLYMORPHISM_COUNT>
		static int getMaxPTMcount()
		{ return iMaxPTMcount;}
		static int getMaxPolymorphismCount()
		{ return iMaxPolymorphismCount;}

		// retrieve <PEPTIDE_MASS_RANGE>
		static double getMaxPeptideMass()
		{ return dMaxPeptideMass;}
		static double getMinPeptideMass()
		{ return dMinPeptideMass;}

		// retrieve <CLEAVAGE_RULES>
		static string getCleavageAfterResidues()
		{ return sCleavageAfterResidues;}
		static string getCleavageBeforeResidues()
		{ return sCleavageBeforeResidues;}
		static int getMaxMissedCleavages()
		{ return iMaxMissedCleavages;}
		static int getMinSpecificCleavages()
		{ return iMinSpecificCleavages;}
		static bool getTestStartRemoval()
		{ return bTestStartRemoval;}

		static bool getPTMinfo(map<char, string> & mPTMinfo);

		// retrieve <ATOM_ISOTOPIC_COMPOSITION>
		// the input character is the atom name CHONPS
		static bool getAtomIsotopicComposition(
				char cAtom,
				vector<double> & vdAtomicMass,
				vector<double> & vdComposition);

		// retrieve <RESIDUE_ATOMIC_COMPOSITION>
		static bool getResidueAtomicComposition( string & sAtomicCompositionTable );

		static Isotopologue configIsotopologue;
		static vector<string> vsSingleResidueNames;
		static vector<double> vdSingleResidueMasses;

		static vector<string> vsDoubleResidueNames;
		static vector<double> vdDoubleResidueMasses;

		static double getResidueMass( string sResidue );

		static double getTerminusMassN()
		{ return dTerminusMassN; }
		static double getTerminusMassC()
		{ return dTerminusMassC; }

		static double getProtonMass()
		{ return 1.0072765; }

		static double getNeutronMass()
		{ return 1.008665; }

		static double getMassH2O()
		{ return dMassH2O; }

		static double getMassNH3()
		{ return dMassNH3; }

		static double dnorm(double mean, double sd, double x)
		{
			double SQRT2PI  = 2.506628;
			double a = (x - mean) / sd;
  			return exp(- a * a / 2) / ( SQRT2PI * sd );
		}

		static double pnorm( double dMean, double dStandardDeviation, double dRandomVariable )
		{
			double dZScore = ( dRandomVariable - dMean ) / dStandardDeviation ;
			double dProbability = 0.5 * erfc( -dZScore / sqrt( 2.0 ) );
			return dProbability;
		}

		static double scoreError( double dMassError)
		{

		//	pnorm function
			return ( 1.0 - pnorm( 0, (getMassAccuracyFragmentIon() / 2), fabs(dMassError) ) ) * 2.0;

		//	dnorm function
		//	return  ( dnorm( 0, (getMassAccuracyFragmentIon() / 2.0), fabs(dMassError) ) ) /
		//			( dnorm( 0, (getMassAccuracyFragmentIon() / 2.0), 0 ) )	;

		//  sigmoid function
		//	return ( 1/(1+exp(dMassError*600-3)));

		}

	protected:
		ProNovoConfig();

	private:

		static ProNovoConfig* ProNovoConfigSingleton;

		// the filename of the configuration file
		static string sFilename;

		// the working directory
		static string sWorkingDirectory;

		// load all static parameters into the memory
		// when openning the file
		bool getParameters( TiXmlDocument & );

		// retrieve the value as string from an element
		// the element is located by giving its hierarchial path
		static string getValue( TiXmlDocument &, const vector<string>& );

		static void replaceDelimitor( string & sLine, char cOldDelimitor, char cNewDelimitor );

		// variables from the PEPTIDE_IDENTIFICATION element
		static string sMSFileType;
		static string sIDFileType;
		static string sFASTAFilename;
		static string sFragmentationMethod;

		static int iMaxPTMcount;
		static int iMaxPolymorphismCount;

		static double dMaxPeptideMass;
		static double dMinPeptideMass;

		static string sCleavageAfterResidues;
		static string sCleavageBeforeResidues;
		static int iMaxMissedCleavages;
		static int iMinSpecificCleavages;
		static bool bTestStartRemoval;

		static double dMassAccuracyParentIon;
		static double dMassAccuracyFragmentIon;
		static vector<int> viParentMassWindows;

		static double dTerminusMassN;
		static double dTerminusMassC;

		static double dMassH2O;
		static double dMassNH3;


};

#endif /*PRONOVOCONFIG_H_*/
