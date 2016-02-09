#include "MassSpectralPeakCluster.h"

MassSpectralPeakCluster::MassSpectralPeakCluster()
{
}

MassSpectralPeakCluster::~MassSpectralPeakCluster()
{
}

bool MassSpectralPeakCluster::testPeakValidity( MassSpectralPeak & peak )
{
	if( peak.getChargeState() > 0 )
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool MassSpectralPeakCluster::testPeakCompatibility( MassSpectralPeak & peak1, MassSpectralPeak & peak2 )
{
	double dMassDifference = fabs(peak1.getNeutralMass() - peak2.getNeutralMass());
	bool bSameChargeState = false;
	if(peak1.getChargeState() == peak2.getChargeState())
	{
		bSameChargeState = true;
	}
	else
	{
		bSameChargeState = false;
	}
	if( bSameChargeState && (fabs(dMassDifference - ProNovoConfig::getNeutronMass()) < ProNovoConfig::getMassAccuracyFragmentIon()
	|| dMassDifference < ProNovoConfig::getMassAccuracyFragmentIon() ))
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool MassSpectralPeakCluster::addPeak( MassSpectralPeak peakInput )
{

	if( !testPeakValidity( peakInput ) )
	{
		return false;
	}

	// if the vPeakList is empty, seed this cluster with this peak
	if(vPeakList.size() == 0)
	{
		vPeakList.push_back( peakInput );
		iChargeState = peakInput.getChargeState();
		return true;
	}


	for( unsigned int i = 0; i < vPeakList.size(); i++ )
	{
		if(testPeakCompatibility(peakInput, vPeakList[i]))
		{
			// this peak can be added into this cluster
			vPeakList.push_back( peakInput );
			return true;
		}
	}
	return false;

}


void MassSpectralPeakCluster::insertPeak( const vector<MassSpectralPeak> & vPeakListInput)
{

	for (unsigned int i = 0; i < vPeakListInput.size(); ++i)
	{
		vPeakList.push_back(vPeakListInput[i]);
	}

}

void MassSpectralPeakCluster::getPeakList(vector<MassSpectralPeak> & vPeakListOutput)
{
	vPeakListOutput = vPeakList;
}

void MassSpectralPeakCluster::sortPeaks()
{
	LessMassSpectralPeak lessPeak("MZ");
	sort(vPeakList.begin(), vPeakList.end(), lessPeak);
}

MassSpectralPeak MassSpectralPeakCluster::getMonoisotopicPeak()
{
	sortPeaks();
	return *vPeakList.begin();
}

double MassSpectralPeakCluster::computeAverageNeutralMonoisotopicMass()
{
	sortPeaks();
	// this is a simplistic method that assumes the monoisotopic mass is the lowest in the cluster
	double dMass = 0.0;
	for( unsigned int i = 0; i < vPeakList.size(); i++ )
	{
		// i equal the number of 1 Da shift
		dMass = dMass + ( vPeakList[i].getNeutralMass() - (double)i )*vPeakList[i].getSignalIntensity();
	}
	dMass = dMass / computeTotalSignalIntensity();

	return dMass;
}

double MassSpectralPeakCluster::computeTotalSignalIntensity()
{
	double dIntensity = 0;
	for( unsigned int i = 0; i < vPeakList.size(); i++ )
	{
		dIntensity = dIntensity + (vPeakList[i].getSignalIntensity() );
	}
	return dIntensity;
}

double MassSpectralPeakCluster::computeAggregateSignalToNoiseRatio()
{
	// take the total signal intensity divided by the average noise level
	double dAggregateSNR = 0;
	double dAverageNoise = 0;
	for( unsigned int i = 0; i < vPeakList.size(); i++ )
	{
		dAverageNoise = dAverageNoise + vPeakList[i].getNoise();
	}
	dAverageNoise = dAverageNoise / (double)vPeakList.size();
	dAggregateSNR = computeTotalSignalIntensity() / dAverageNoise ;

	return dAggregateSNR;
}

double MassSpectralPeakCluster::computeAverageResolution()
{
	double dResolution = 0;
	for( unsigned int i = 0; i < vPeakList.size(); i++ )
	{
		dResolution = dResolution + (vPeakList[i].getResolution())*(vPeakList[i].getSignalIntensity());
	}
	dResolution = dResolution / computeTotalSignalIntensity();
	return dResolution;
}

int MassSpectralPeakCluster::getPeakCount()
{
	return vPeakList.size();
}

void MassSpectralPeakCluster::print()
{
	cout << "###" << endl;
	for (unsigned int i = 0; i < vPeakList.size(); ++i)
	{
		cout << "(" << vPeakList[i].getNeutralMass() <<  ",	" << vPeakList[i].getMZ() << ",	"
		<< vPeakList[i].getIntensity()<<  ",	" << vPeakList[i].getChargeState()   << ")" << endl;
	}

}

bool LessMassSpectralPeakCluster::operator() ( MassSpectralPeakCluster cluster1, MassSpectralPeakCluster cluster2 ) const
{
	if( sKey == "Mass" )
	{
		if( cluster1.computeAverageNeutralMonoisotopicMass() < cluster2.computeAverageNeutralMonoisotopicMass() )
			return true;
		else
			return false;
	}
	else if ( sKey == "Intensity" )
	{
		if( cluster1.computeTotalSignalIntensity() < cluster2.computeTotalSignalIntensity() )
			return true;
		else
			return false;
	}
	else
	{
		if( cluster1.computeAverageNeutralMonoisotopicMass() < cluster2.computeAverageNeutralMonoisotopicMass() )
			return true;
		else
			return false;
	}
}
