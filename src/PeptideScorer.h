#ifndef PEPTIDESCORER_H_
#define PEPTIDESCORER_H_

#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include "proNovoConfig.h"
#include "TandemMassSpectrum.h"
#define BIN_RES 1000

using namespace std;

class ProductIon
{
	char cIonType;
	int iIonNumber;
	int iCharge;
	int iMostAbundantPeakIndex;
	double dMostAbundantMass;
	double dMostAbundantMZ;
	double dMZError;
	double dMassError;
	double dScoreWeight;
	bool bComplementaryFragmentObserved;
public:
	ProductIon();
	~ProductIon();
	// set the basic info
	void setProductIon(char cIonTypeInput, int iIonNumberInput, int iChargeInput);
	// set info for the found
	void setObservedInfo(double dMZErrorInput, double dWeightInput,
			double dMostAbundantMZInput, int iMostAbundantPeakIndexInput);
	void setComplementaryFragmentObserved(bool bComplementaryFragmentObservedInput);
	char getIonType() {return cIonType;};
	int getIonNumber() {return iIonNumber;};
	int getCharge() {return iCharge;};
	double getMZError() {return dMZError;};
	double getMassError() {return dMassError;};
	double getScoreWeight() {return dScoreWeight;};
	double getMostAbundantMass() {return dMostAbundantMass;};
	double getMostAbundantMZ() {return dMostAbundantMZ;};
	int getMostAbundantPeakIndex() {return iMostAbundantPeakIndex;};
	bool getComplementaryFragmentObserved() {return bComplementaryFragmentObserved;};
};


class LessProductIon
{
	string sKey;
public:
	LessProductIon(string sKeyInput){sKey = sKeyInput;};
	bool operator()(ProductIon ion1, ProductIon ion2 ) const;
};

class PeptideScorer
{
	/*
	 * configurations
	 */
	double dMassTolerance;
	map<char, double> mapResidueMass;
	double dNtermMass;
	double dCtermMass;
	double dProtonMass;

	/*
	 * the current MS2 to be used for scoring peptides
	 */
	vector<double> vdMZ;
	vector<double> vdIntensity;
	vector<double> vdNormalizedIntensity;
	vector<int> viZinput; 						// charge state from input
	vector<int> vbPeakPresenceBins;
	double dParentMZ;
	int iParentChargeState;


	/*
	 * MS2 data and functions for preliminary scoring function
	 */
	// monoisotopic mass and intensity
	vector<double> vdMonoMass;
	vector<double> vdMonoInt;
	vector<int> vdMonoZ;


	// vdYMass = masses for y1, y2, y3 ... in this order
	// vdBMass = masses for b1, b2, b3 ... in this order
	bool computeYBionMass(const string & sPeptide,
			vector<double> & vdYmass, vector<double> & vdBmass,
			double & dPeptideMass);

	// binarySearch returns true if found at least one value
	// vdList have to be in ascending order
	// indices for the found values in vdList are returned in viIndex4Found
	bool binarySearch(const double & dTarget,
			const vector<double> & vdList, const double & dTolerance,
			vector<int> & viIndex4Found );

	bool searchMZ(const double & dTargetMZ, int & iIndex4Found );

	/*
	 * MS2 data and functions for primary scoring function
	 */
	// get the index for the maximum value in vdData
	int getMaxValueIndex(const vector<double> & vdData);
	// get the index for the maximum value in the values within viIndexRange of vdData
	// viIndexRange is a list of indices for vdData
	int getMaxValueIndex(const vector<double> & vdData, const vector<int> & viIndexRange);


	bool findProductIon(const vector<double> & vdIonMass,
						 const vector<double> & vdIonProb,
						 const int & iCharge,
						 double & dScoreWeight,
						 double & dAverageMZError,
						 double & dMostAbundantObservedMZ,
						 int & iMostAbundantPeakIndex );

	// identical to findProductIon except searching for NH3 loss and water loss
	bool findProductIon2(const vector<double> & vdIonMass,
						 const vector<double> & vdIonProb,
						 const int & iCharge,
						 double & dScoreWeight,
						 double & dAverageMZError,
						 double & dMostAbundantObservedMZ,
						 int & iMostAbundantPeakIndex );


public:
	PeptideScorer();
	~PeptideScorer();

	// set current MS2 scan: vdMZ, vdIntensity, dParentMZ, iChargeState
	// calling this function clears the data for previous MS2 scan
	// return false, if there is a problem with the mass spectrum
	bool setMS2(
			const vector<double> & vdMZinput,
			const vector<double> & vdIntensityInput,
			const vector<int> & viZinputTemp,
			const double & dParentMZinput,
			const int & iParentChargeStateInput);

	// PTM must be a non-alphabetic symbol
	// preliminary scoring function for a input peptide
	// return the preliminary score

	// time = 13
	double prelimaryScorePeptide( const string & sPeptide);
	// based on prelimaryScorePeptide, add searches for undeisotoped +1 ions and add intensity scoring
	// time = 16
	double prelimaryScorePeptide2( const string & sPeptide);
	// searches for +1 ions and +2 ions and their water and NH3 loss peaks, add intensity scoring
	// no searches in MonoMass list
	// time = 11
	double prelimaryScorePeptide3( const string & sPeptide);


	// primary scoring function for a input peptide
	// return the primary score
	double primaryScorePetide( const string & sPeptide);
	// identical to primaryScorePetide, except adding intensity scoring
	double primaryScorePetide2( const string & sPeptide);
	// identical to primaryScorePetide2, except using polyaveragine and thus being 5 times faster
	double primaryScorePetide3( const string & sPeptide);
	// identical to primaryScorePetide3, except searching for NH3 and water loss
	double primaryScorePetide4( const string & sPeptide);


};

#endif /* PEPTIDESCORER_H_ */
