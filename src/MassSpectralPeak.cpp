#include "MassSpectralPeak.h"

MassSpectralPeak::MassSpectralPeak()
{
	dMZ = 0;
	dIntensity = 0;
	dResolution = 0;
	dBaseline = 0;
	dNoise = 0;
	iChargeState = 0;
}

MassSpectralPeak::MassSpectralPeak(double dMZInput, double dIntensityInput, double dResolutionInput, double dBaselineInput, double dNoiseInput, int iChargeStateInput)
{
	dMZ = dMZInput;
	dIntensity = dIntensityInput;
	dResolution = dResolutionInput;
	dBaseline = dBaselineInput;
	dNoise = dNoiseInput;
	iChargeState = iChargeStateInput;
}

MassSpectralPeak::~MassSpectralPeak()
{
}

double MassSpectralPeak::getMZ()
{
	return dMZ;
}

double MassSpectralPeak::getIntensity()
{
	return dIntensity;	
}

double MassSpectralPeak::getResolution()
{
	return dResolution;	
}

double MassSpectralPeak::getBaseline()
{
	return dBaseline;
}

double MassSpectralPeak::getNoise()
{
	return dNoise;
}

int MassSpectralPeak::getChargeState()
{
	return iChargeState;
}

double MassSpectralPeak::getNeutralMass()
{
	return (iChargeState * dMZ - iChargeState * ProNovoConfig::getProtonMass());
}

double MassSpectralPeak::getSignalIntensity()
{
	return (dIntensity - dBaseline);
}

double MassSpectralPeak::getSignalToNoiseRatio()
{
	if(dNoise == 0 )
	{
		return 0;
	}
	return (getSignalIntensity()/dNoise);
	
}
	
void MassSpectralPeak::setMZ(double dMZInput)
{
	dMZ = dMZInput;
}

void MassSpectralPeak::setIntensity(double dIntensityInput)
{
	dIntensity = dIntensityInput;
}

void MassSpectralPeak::setResolution(double dResolutionInput)
{
	dResolution = dResolutionInput;
}

void MassSpectralPeak::setBaseline(double dBaselineInput)
{
	dBaseline = dBaselineInput;	
}

void MassSpectralPeak::setNoise(double dNoiseInput)
{
	dNoise = dNoiseInput;	
}

void MassSpectralPeak::setChargeState(int iChargeStateInput)
{
	iChargeState = iChargeStateInput;	
}

bool LessMassSpectralPeak::operator() ( MassSpectralPeak peak1, MassSpectralPeak peak2 ) const
{
	if( sKey == "MZ" )
	{
		if( peak1.getMZ() < peak2.getMZ() )
			return true;
		else
			return false;
	}
	else if ( sKey == "Intensity" )
	{
		if( peak1.getIntensity() < peak2.getIntensity() )
			return true;
		else
			return false;
	}
	else
	{
		if( peak1.getMZ() < peak2.getMZ() )
			return true;
		else
			return false;
	}
}

bool LessMassSpectralPeak::operator() ( MassSpectralPeak * peak1, MassSpectralPeak * peak2 ) const
{
	if( sKey == "MZ" )
	{
		if( peak1->getMZ() < peak2->getMZ() )
			return true;
		else
			return false;
	}
	else if ( sKey == "Intensity" )
	{
		if( peak1->getIntensity() < peak2->getIntensity() )
			return true;
		else
			return false;
	}
	else
	{
		if( peak1->getMZ() < peak2->getMZ() )
			return true;
		else
			return false;
	}	
}

