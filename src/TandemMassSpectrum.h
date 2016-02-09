#ifndef TANDEMMASSSPECTRUM_H_
#define TANDEMMASSSPECTRUM_H_

#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <iostream>
#include <sstream>
#include "MassSpectralPeak.h"
#include "MassSpectralPeakCluster.h"

#define MAX_Z 10

using namespace std;

class TandemMassSpectrum
{
	vector<MassSpectralPeak> vPeakList;
	vector<MassSpectralPeakCluster> vPeakClusterList;
	int iScanNumber;
	double dParentMZ;
	int iParentChargeState;
	
	// the two functions below are for determination of charge state and not in use
	void extractOneCluster(vector<MassSpectralPeak *> & currentPeakList);	
	double findPeaksInCluster(vector<MassSpectralPeak *> & currentPeakList, 
			unsigned int iChargeState, MassSpectralPeak * maxPeak, bool bAddCluster);
	
public:
	TandemMassSpectrum();
	~TandemMassSpectrum();
	
	int getScanNumber();
	double getParentMZ();
	int getParentChargeState();
	
	// note that parent mass could be zero if the charge state is zero.
	double getParentMass();
	
	void setScanNumber(int iScanNumberInput);
	void setParentMZ(double dParentMZInput);
	void setParentChargeState(int iParentChargeStateInput);
	
	void addPeak(MassSpectralPeak MassSpectralPeakInput);
	void addPeak(double dMZInput, double dIntensityInput, double dResolutionInput, 
			double dBaselineInput, double dNoiseInput, int iChargeStateInput);
	
	void clearPeakList();
	
	// populate vPeakClusterList
	void clusterPeaks(); 
	vector<MassSpectralPeakCluster> getPeakClusterList();
	
	void clearPeakClusterList();
	
	void clearMassSpectrum();	
	
	void print();
	
};

#endif /*TANDEMMASSSPECTRUM_H_*/
