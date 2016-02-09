#ifndef MASSSPECTRALPEAK_H_
#define MASSSPECTRALPEAK_H_
#include <string>
#include "proNovoConfig.h"
using namespace std;

class MassSpectralPeak
{
	
	double dMZ;
	double dIntensity;
	double dResolution;
	double dBaseline;
	double dNoise;
	int iChargeState;
	
public:
	MassSpectralPeak();
	MassSpectralPeak(double dMZInput, double dIntensityInput, double dResolutionInput, 
			double dBaselineInput, double dNoiseInput, int iChargeStateInput);
	~MassSpectralPeak();
	
	double getMZ();
	double getIntensity();
	double getResolution();
	double getBaseline();
	double getNoise();
	int getChargeState();
	
	double getNeutralMass();
	double getSignalIntensity();
	double getSignalToNoiseRatio();
	
	void setMZ(double dMZInput);
	void setIntensity(double dIntensityInput);
	void setResolution(double dResolutionInput);
	void setBaseline(double dBaselineInput);
	void setNoise(double dNoiseInput);
	void setChargeState(int iChargeStateInput);
	
	
};

class LessMassSpectralPeak
{
	string sKey;

public:
	LessMassSpectralPeak();
	LessMassSpectralPeak(string sKeyInput){ sKey = sKeyInput; }
	void setKey( string sKeyInput ) { sKey = sKeyInput; }
	string getKey() { return sKey; }
	
	/*
	 * the MassSpectralPeak objects can be sorted by a number of member variables:
	 * the sKey specifies which variables to be used to sort
	 * the accepted values for sKay  are
	 * "MZ"			
	 * "Intensity"	
	 */
		
	bool operator() ( MassSpectralPeak peak1, MassSpectralPeak peak2 ) const;	
	bool operator() ( MassSpectralPeak * peak1, MassSpectralPeak * peak2 ) const;	
};

#endif /*MASSSPECTRALPEAK_H_*/
