
#include "proNovoConfig.h"

string ProNovoConfig::sFilename = "ProNovoConfig.xml";

#if _WIN32
	string ProNovoConfig::sWorkingDirectory = ".\\";
#else
	string ProNovoConfig::sWorkingDirectory = ".";
#endif

// variables from the PEPTIDE_IDENTIFICATION element
string ProNovoConfig::sMSFileType = "FT2";
string ProNovoConfig::sIDFileType = "SQT";
string ProNovoConfig::sFASTAFilename = "";
string ProNovoConfig::sFragmentationMethod = "CID";

int ProNovoConfig::iMaxPTMcount = 0;
int ProNovoConfig::iMaxPolymorphismCount = 0;

double ProNovoConfig::dMaxPeptideMass = 0;
double ProNovoConfig::dMinPeptideMass = 10000;

string ProNovoConfig::sCleavageAfterResidues = "KR";
string ProNovoConfig::sCleavageBeforeResidues = "ACDEFGHIJKLMNPQRSTVWXY";
int ProNovoConfig::iMaxMissedCleavages = 2;
int ProNovoConfig::iMinSpecificCleavages = 2;
bool ProNovoConfig::bTestStartRemoval = false;

double ProNovoConfig::dMassAccuracyParentIon = 0.05;
double ProNovoConfig::dMassAccuracyFragmentIon = 0.05;
vector<int> ProNovoConfig::viParentMassWindows;

ProNovoConfig* ProNovoConfig::ProNovoConfigSingleton = 0;

vector<string> ProNovoConfig::vsSingleResidueNames;
vector<double> ProNovoConfig::vdSingleResidueMasses;

vector<string> ProNovoConfig::vsDoubleResidueNames;
vector<double> ProNovoConfig::vdDoubleResidueMasses;

double ProNovoConfig::dTerminusMassN = 1.0072765;
double ProNovoConfig::dTerminusMassC = 17.0265;

double ProNovoConfig::dMassH2O = 18.0106;
double ProNovoConfig::dMassNH3 = 17.0265;

Isotopologue ProNovoConfig::configIsotopologue;

ProNovoConfig::ProNovoConfig()
{}

bool ProNovoConfig::setFilename( const string & sConfigFileName )
{
	if ( ProNovoConfigSingleton == 0 )
	{
		ProNovoConfigSingleton = new ProNovoConfig;
	}

	sFilename = sConfigFileName;

	// Creat a TinyXML document for ProNovoConfig.XML
	TiXmlDocument txdConfigFile;

	// Try loading the file.
	if ( ! ( txdConfigFile.LoadFile( sFilename.c_str() ) ) )
	{
		cout << "ERROR! Loading Configuration file" << endl;
		return false;
	}

	if( !ProNovoConfigSingleton->getParameters( txdConfigFile ) )
	{
		return false;
	}
	// If everything goes fine return 0.
	return true;

}

bool ProNovoConfig::setWorkingDirectory( const string & sDirectoryName )
{
	if( sDirectoryName[ sDirectoryName.size() - 1 ] == ProNovoConfig::getSeparator() )
	{
		sWorkingDirectory = sDirectoryName;
	}
	else
	{
		sWorkingDirectory = sDirectoryName + ProNovoConfig::getSeparator();
	}

	return true;
}

char ProNovoConfig::getSeparator()
{
#if _WIN32
		return '\\' ;
#else
		return '/' ;
#endif
}

bool ProNovoConfig::getAtomIsotopicComposition( char cAtom,
		vector<double> & vdAtomicMass,
		vector<double> & vdComposition)
{

	// clear the input vectors
	vdAtomicMass.clear();
	vdComposition.clear();

	// Creat a TinyXML document for ProNovoConfig.XML
	TiXmlDocument txdConfigFile;

	// Try loading the file.
	if ( ! ( txdConfigFile.LoadFile( sFilename.c_str() ) ) )
	{
		cout << "ERROR! Loading Configuration file" << endl;
		return false;
	}

	string sData;
	istringstream issStream;
	double dValue;
	string sAtom = "X";
	sAtom[0] = cAtom;

	// creat the path to the MASS_DA element of the input sAtom
	vector<string> vsTagList;
	vsTagList.push_back( "PRONOVO_CONFIG" );
	vsTagList.push_back( "PEPTIDE_IDENTIFICATION" );
	vsTagList.push_back( "ATOM_ISOTOPIC_COMPOSITION" );
	vsTagList.push_back( sAtom );
	vsTagList.push_back( "MASS_DA" );

	// get the text inside and extract the value
	sData = getValue( txdConfigFile, vsTagList );
	replaceDelimitor( sData, ',', '\t' );
	// clear end of file state
	issStream.clear();
	// re-set the string associated with issStream
	issStream.str( sData );
	while( !( issStream.eof() ) )
	{
		issStream >> dValue;
		vdAtomicMass.push_back( dValue );
	}

	// move to the PERCENTAGE element
	vsTagList.pop_back();
	vsTagList.push_back( "PERCENTAGE" );
	sData = getValue( txdConfigFile, vsTagList );
	replaceDelimitor( sData, ',', '\t' );
	issStream.clear();
	issStream.str( sData );
	while( !( issStream.eof() ) )
	{
		issStream >> dValue;
		vdComposition.push_back( dValue );
	}

	return true;
}

bool ProNovoConfig::getResidueAtomicComposition(string & sAtomicCompositionTable)
{
	sAtomicCompositionTable = "";

	// Creat a TinyXML document for ProNovoConfig.XML
	TiXmlDocument txdConfigFile;

	// Try loading the file.
	if ( ! ( txdConfigFile.LoadFile( sFilename.c_str() ) ) )
	{
		cout << "ERROR! Loading Configuration file" << endl;
		return false;
	}

	/*
	 * move the node pointer to RESIDUE_ATOMIC_COMPOSITION
	 * if failed, return an empty map
	 */
	TiXmlNode * txnTemp = NULL;

	txnTemp = txdConfigFile.FirstChild( "PRONOVO_CONFIG" );
	if ( ! txnTemp )
		return false;

	txnTemp = txnTemp->FirstChild( "PEPTIDE_IDENTIFICATION" );
	if ( ! txnTemp )
		return false;

	txnTemp = txnTemp->FirstChild( "RESIDUE_ATOMIC_COMPOSITION" );
	if ( ! txnTemp )
		return false;

	// a node pointer to a atomic composition element
	TiXmlNode * txnTable = NULL;


	// a text pointer for retrieving the text inside ISOTOPOLOGUE element
	TiXmlText * txsText = NULL;

		/*
		 * loop thru all the text nodes inside ISOTOPOLOGUE;
		 * the text nodes can be separated by the comment nodes or other
		 * cast the node to a text node, only if it is of type TEXT, which equals 4
		 * then concatenate all the text
		 */
		for( txnTable = txnTemp->FirstChild("R"); txnTable; txnTable = txnTable->NextSibling("R") )
		{
			txsText =  txnTable->FirstChild()->ToText();
			sAtomicCompositionTable.append( txsText->Value() );
			sAtomicCompositionTable.append( "\n" );

		}
		replaceDelimitor( sAtomicCompositionTable, ',', '\t' );


	return true;

}

bool ProNovoConfig::getPTMinfo(map<char, string> & mPTMinfo)
{
	mPTMinfo.clear();

	// Creat a TinyXML document for ProNovoConfig.XML
	TiXmlDocument txdConfigFile;

	// Try loading the file.
	if ( ! ( txdConfigFile.LoadFile( sFilename.c_str() ) ) )
	{
		cout << "ERROR! Loading Configuration file" << endl;
		return false;
	}

	/*
	 * move the node pointer to RESIDUE_ATOMIC_COMPOSITION
	 * if failed, return an empty map
	 */
	TiXmlNode * txnTemp = NULL;

	txnTemp = txdConfigFile.FirstChild( "PRONOVO_CONFIG" );
	if ( ! txnTemp )
		return false;

	txnTemp = txnTemp->FirstChild( "PEPTIDE_IDENTIFICATION" );
	if ( ! txnTemp )
		return false;

	txnTemp = txnTemp->FirstChild( "PTM" );
	if ( ! txnTemp )
		return false;

	// a node pointer to a atomic composition element
	TiXmlNode * txnTable = NULL;


	// a text pointer for retrieving the text inside ISOTOPOLOGUE element
	TiXmlText * txsText = NULL;
		/*
		 * loop thru all the text nodes inside ISOTOPOLOGUE;
		 * the text nodes can be separated by the comment nodes or other
		 * cast the node to a text node, only if it is of type TEXT, which equals 4
		 * then concatenate all the text
		 */
		for( txnTable = txnTemp->FirstChild("M"); txnTable; txnTable = txnTable->NextSibling("M") )
		{
			txsText =  txnTable->FirstChild("SYMBOL")->FirstChild()->ToText();

			string sSymbol = txsText->Value();

			if(sSymbol.length() != 1)
			{
				if(isalpha(sSymbol[0]))
				{
					cout << "ERROR: An invalid symbol for PTM " << sSymbol << endl;
					continue;
				}
			}

			txsText =  txnTable->FirstChild("RESIDUES")->FirstChild()->ToText();
			string sResidues = txsText->Value();

		//	cout << sSymbol << " -- " << sResidues << endl;
			mPTMinfo[sSymbol[0]] = sResidues;

		}

	return true;

}

bool ProNovoConfig::getParameters( TiXmlDocument & txdConfigFile )
{
	// strings used to specify the path
	string sMainTag = "PRONOVO_CONFIG";
	string sModuleTag;

	// push back the element name in the hierarchical order
	// the top level goes first and the leaf node goes last
	vector<string> vsTagList;

	string sTemp;
	istringstream issStream;

	// Extract the elements inside <PEPTIDE_IDENTIFICATION>
	sModuleTag = "PEPTIDE_IDENTIFICATION";

	vsTagList.clear();
	vsTagList.push_back( sMainTag );
	vsTagList.push_back( sModuleTag );
	vsTagList.push_back( "INPUT_MS_FILE_TYPE" );
	sMSFileType = getValue( txdConfigFile, vsTagList );

	vsTagList[2] = "OUTPUT_ID_FILE_TYPE";
	sIDFileType = getValue( txdConfigFile, vsTagList );

	vsTagList[2] = "FASTA_DATABASE";
	sFASTAFilename = getValue( txdConfigFile, vsTagList );

	vsTagList[2] = "FRAGMENTATION_METHOD";
	sFragmentationMethod = getValue( txdConfigFile, vsTagList );

	vsTagList[2] = "MAX_PTM_COUNT";
	sTemp = getValue( txdConfigFile, vsTagList );
	issStream.clear();
	issStream.str( sTemp );
	issStream >> iMaxPTMcount;

	vsTagList[2] = "MAX_POLYMORPHISM_COUNT";
	sTemp = getValue( txdConfigFile, vsTagList );
	issStream.clear();
	issStream.str( sTemp );
	issStream >> iMaxPolymorphismCount;

	vsTagList[2] = "MASS_ACCURACY";
	vsTagList.push_back( "PARENT_ION" );
	sTemp = getValue( txdConfigFile, vsTagList );
	issStream.clear();
	issStream.str( sTemp );
	issStream >> dMassAccuracyParentIon;

	if( dMassAccuracyParentIon > 0.1)
	{
		cout << "Please provide high-resolution data with a parent ion mass accuracy better than 0.05!" << endl;
		cout << "Job aborted." << endl;
		return false;
	}

	vsTagList[3] = "FRAGMENT_IONS";
	sTemp = getValue( txdConfigFile, vsTagList ) ;
	issStream.clear();
	issStream.str( sTemp );
	issStream >> dMassAccuracyFragmentIon;
	if( dMassAccuracyFragmentIon > 0.1)
	{
		cout << "Please provide high-resolution data with a product ion mass accuracy better than 0.05!" << endl;
		cout << "Job aborted." << endl;
		return false;
	}

	vsTagList[3] = "PARENT_MASS_WINDOWS";
	sTemp = getValue( txdConfigFile, vsTagList ) ;
	issStream.clear();
	issStream.str( sTemp );
	string sField;
	viParentMassWindows.clear();
	while(getline(issStream, sField, ','))
	{
		istringstream issField(sField);
		int iWindow;
		issField >> iWindow;
		viParentMassWindows.push_back(iWindow);
	}

	// read PEPTIDE_MASS_RANGE
	vsTagList[2] = "PEPTIDE_MASS_RANGE";
	vsTagList[3] = "MIN";
	sTemp = getValue( txdConfigFile, vsTagList );
	issStream.clear();
	issStream.str( sTemp );
	issStream >> dMinPeptideMass;

	vsTagList[3] = "MAX";
	sTemp = getValue( txdConfigFile, vsTagList ) ;
	issStream.clear();
	issStream.str( sTemp );
	issStream >> dMaxPeptideMass;

	// read <CLEAVAGE_RULES>
	vsTagList[2] = "CLEAVAGE_RULES";
	vsTagList[3] = "AFTER_RESIDUES";
	sCleavageAfterResidues = getValue( txdConfigFile, vsTagList );

	vsTagList[3] = "BEFORE_RESIDUES";
	sCleavageBeforeResidues = getValue( txdConfigFile, vsTagList );

	vsTagList[3] = "MAX_MISSED_CLEAVAGE";
	sTemp = getValue( txdConfigFile, vsTagList ) ;
	issStream.clear();
	issStream.str( sTemp );
	issStream >> iMaxMissedCleavages;

	vsTagList[3] = "MIN_SPECIFIC_CLEAVAGE";
	sTemp = getValue( txdConfigFile, vsTagList ) ;
	issStream.clear();
	issStream.str( sTemp );
	issStream >> iMinSpecificCleavages;

	vsTagList[3] = "TEST_START_REMOVAL";
	sTemp = getValue( txdConfigFile, vsTagList );
	if(sTemp == "TRUE" || sTemp == "True" || sTemp == "true" ||sTemp == "T" )
		bTestStartRemoval = true;
	else
		bTestStartRemoval = false;

	// setup configIsotopologue, this must be done after the other parameters have been configured.
	string sResidueAtomicComposition;
	getResidueAtomicComposition(sResidueAtomicComposition);
	configIsotopologue.setupIsotopologue(sResidueAtomicComposition);
	configIsotopologue.getSingleResidueMostAbundantMasses(vsSingleResidueNames, vdSingleResidueMasses, dTerminusMassN, dTerminusMassC);
//	configIsotopologue.getDoubleResidueMonoisotopicMasses(vsDoubleResidueNames, vdDoubleResidueMasses);
	dMassH2O = configIsotopologue.computeMassH2O();
	dMassNH3 = configIsotopologue.computeMassNH3();

	return true;
}

double ProNovoConfig::getResidueMass( string sResidue )
{
	unsigned int i;
	double dResidueMass = 0.0;
	if( sResidue == "|||")
	{
		dResidueMass = 0.0;
		return dResidueMass;
	}
	for (i = 0; i < vsSingleResidueNames.size(); ++i)
	{
		if( vsSingleResidueNames[i] == sResidue )
		{
			dResidueMass = vdSingleResidueMasses[i];
			return dResidueMass;
		}
	}
	string sResidueReverse = sResidue;
	reverse( sResidueReverse.begin(), sResidueReverse.end());
	for (i = 0; i < vsDoubleResidueNames.size(); ++i)
	{
		if( vsDoubleResidueNames[i] == sResidue || vsDoubleResidueNames[i] == sResidueReverse )
		{
			dResidueMass = vdDoubleResidueMasses[i];
			return dResidueMass;
		}
	}
	cout << "ERROR: cannot find residue " << sResidue << endl;
	return dResidueMass;
}

/*
 * An utility method to safely extract information from an XML tag.
 * the text inside an XML element is extracted and pasted together if separated
 * by comments or elements.
 * the element is reached by giving a vector of the element name in the order
 * of their hierarchy. Arbitrary number of level can be reached
 */

string ProNovoConfig::getValue( TiXmlDocument &txdDoc,
		const vector<string> &vsTagList )
{
	// Check to see if the provided XML node is valid.
	// If yes, extract the value from the node return it.
	// If no, return emply string.

	// creat a pointer to a node
	TiXmlNode * txnTemp = NULL;

	// check if the tree path is not empty
	if ( vsTagList.size() < 1 )
		return string("");

	// iterator for the input vsTagList
	vector<string>::const_iterator itrTagListItr;
	itrTagListItr = vsTagList.begin();

	// move the pointer to txddoc's first child with the specified tag name
	txnTemp = txdDoc.FirstChild( (*itrTagListItr ).c_str() );

	// check if this element exists
	if ( ! txnTemp )
	{
		cout << "ERROR: TAG\"" << (*itrTagListItr) <<
			"\" not found in the configuration file." << endl;
		return string("");
	}

	itrTagListItr++;

	// move the pointer down the hierarchial tree of elements
	for( ; itrTagListItr != vsTagList.end(); itrTagListItr++ )
	{

		txnTemp = txnTemp->FirstChild( (*itrTagListItr ).c_str() );

		if ( ! txnTemp )
		{
			cout << "ERROR: TAG\"" << (*itrTagListItr) <<
				"\" not found in the configuration file." << endl;
			return string("");
		}

	}

	// move the iterator back to point it to the last element name
	itrTagListItr--;

	/*
	 * inside the pointed element, there could be a mixture of
	 * text nodes, comment nodes and element nodes
	 * loop thru every nodes and for each text nodes, retrieve their text
	 * concatenate the text together and return them
	 */

	TiXmlText *txs;
	string sTemp = "";

	// point txnTemp to the child nodes and loop thru every child node
	for( txnTemp = txnTemp->FirstChild(); txnTemp; txnTemp = txnTemp->NextSibling() )
	{
		// if this node is pointing to a node of type TEXT, which equals 4 in enum NodeType
		if( txnTemp->Type() == 4 )
		{
			// cast txnTemp to a text node
			txs = txnTemp->ToText();
			// get txnTemp's value and then append it to sTemp
			if( txs )
				sTemp.append( txs->Value() );
		}
	}

	return sTemp;
}

void ProNovoConfig::replaceDelimitor( string & sLine, char cOldDelimitor, char cNewDelimitor )
{
	int iLength = sLine.length();
	for( int i = 0; i < iLength; ++i )
	{
		if( sLine[i] == cOldDelimitor )
			sLine[i] = cNewDelimitor;
	}
	return;
}

