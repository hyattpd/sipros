#ifndef MASSSPECTRALPEAKCLUSTER_H_
#define MASSSPECTRALPEAKCLUSTER_H_

#include <string>
#include <vector>
#include <math.h>
#include "MassSpectralPeak.h"
#include "proNovoConfig.h"

using namespace std;

class MassSpectralPeakCluster
{
	vector<MassSpectralPeak> vPeakList;
	int iChargeState;
	bool testPeakValidity( MassSpectralPeak & peak );
	bool testPeakCompatibility( MassSpectralPeak & peak1, MassSpectralPeak & peak2 );

public:
	MassSpectralPeakCluster();
	~MassSpectralPeakCluster();

	bool addPeak( MassSpectralPeak peakInput );

	void insertPeak( const MassSpectralPeak & peakInput ){vPeakList.push_back(peakInput);};

	void insertPeak( const vector<MassSpectralPeak> & vPeakListInput);

	void getPeakList(vector<MassSpectralPeak> & vPeakListOutput);

	void sortPeaks(); // by MZ

	MassSpectralPeak getMonoisotopicPeak();
	int getPeakCount();
	int getChargeState() {return iChargeState;};
	double computeAverageNeutralMonoisotopicMass();
	double computeTotalSignalIntensity();
	double computeAggregateSignalToNoiseRatio();
	double computeAverageResolution();
	void print();

};

class LessMassSpectralPeakCluster
{
	string sKey;

public:
	LessMassSpectralPeakCluster();
	LessMassSpectralPeakCluster(string sKeyInput){ sKey = sKeyInput; }
	void setKey( string sKeyInput ) { sKey = sKeyInput; }
	string getKey() { return sKey; }

	/*
	 * the MassSpectralPeak objects can be sorted by a number of member variables:
	 * the sKey specifies which variables to be used to sort
	 * the accepted values for sKay  are
	 * "Mass"
	 * "Intensity"
	 */

	bool operator() ( MassSpectralPeakCluster cluster1, MassSpectralPeakCluster cluster2 ) const;

};

#endif /*MASSSPECTRALPEAKCLUSTER_H_*/
