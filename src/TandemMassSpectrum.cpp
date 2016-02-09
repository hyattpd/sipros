#include "TandemMassSpectrum.h"

TandemMassSpectrum::TandemMassSpectrum()
{
	iScanNumber = 0;
	dParentMZ = 0;
	iParentChargeState = 0;
}

TandemMassSpectrum::~TandemMassSpectrum()
{
}

int TandemMassSpectrum::getScanNumber()
{
	return iScanNumber;
}

double TandemMassSpectrum::getParentMZ()
{
	return dParentMZ;
}

int TandemMassSpectrum::getParentChargeState()
{
	return iParentChargeState;
}

double TandemMassSpectrum::getParentMass()
{
	double dParentMass = dParentMZ*iParentChargeState - ProNovoConfig::getProtonMass()*iParentChargeState;
	return dParentMass;
}

void TandemMassSpectrum::setScanNumber(int iScanNumberInput)
{
	iScanNumber = iScanNumberInput;
	clearPeakList();
}

void TandemMassSpectrum::setParentMZ(double dParentMZInput)
{
	dParentMZ = dParentMZInput;
}

void TandemMassSpectrum::setParentChargeState(int iParentChargeStateInput)
{
	iParentChargeState = iParentChargeStateInput;
}

void TandemMassSpectrum::addPeak(MassSpectralPeak MassSpectralPeakInput)
{
	vPeakList.push_back(MassSpectralPeakInput);
}

void TandemMassSpectrum::addPeak(double dMZInput, double dIntensityInput, double dResolutionInput, double dBaselineInput, double dNoiseInput, int iChargeStateInput)
{
	MassSpectralPeak MassSpectralPeakInput(dMZInput, dIntensityInput, dResolutionInput, dBaselineInput, dNoiseInput, iChargeStateInput);
	vPeakList.push_back(MassSpectralPeakInput);
}

void TandemMassSpectrum::clearPeakList()
{
	vPeakList.clear();
}


void TandemMassSpectrum::clusterPeaks()
{

	// sort the vPeakList by MZ
	LessMassSpectralPeak lessPeak("MZ");
	sort(vPeakList.begin(), vPeakList.end(), lessPeak);


	unsigned int i;
	unsigned int j;

	bool bAdd2ExistingCluster = false;
	// cluster peaks, starting from the lowest MZ peak to the highest
	for (i = 0; i < vPeakList.size(); ++i)
	{
		bAdd2ExistingCluster = false;

		// try if this peak can be added to any of existing clusters
		for( j = 0; j < vPeakClusterList.size(); j++ )
		{
			if( vPeakClusterList[j].addPeak( vPeakList[i] ) )
			{
				bAdd2ExistingCluster = true;
				break;
			}
		}

		if( !bAdd2ExistingCluster )
		{
			// try if this peak can be used to seed a new cluster
			MassSpectralPeakCluster newCluster;
			if( newCluster.addPeak(  vPeakList[i] ) )
			{
				vPeakClusterList.push_back( newCluster );
			}
		}

	}

//	unsigned int i;
//	unsigned int j;
//
//	vector<MassSpectralPeak *> currentPeakList;
//	for (i = 0; i < vPeakList.size(); ++i)
//	{
//		vPeakList[i].setChargeState(0);
//		currentPeakList.push_back( & vPeakList[i] );
//	}
//
//	while( currentPeakList.size() > 0 )
//	{
//		// currentPeakList is shortened
//		// one or no cluster is added to vPeakClusterList
//		extractOneCluster(currentPeakList);
//	}
//


	// consolidate peak cluster list e.g. with different charge state
	/*
	vector<MassSpectralPeakCluster> vPeakClusterListCopy = vPeakClusterList;
	vPeakClusterList.clear();
	bool bMerged;
	double dMassDifference;
	for (i = 0; i < vPeakClusterListCopy.size(); ++i)
	{
		bMerged = false;
		for (j = 0; j < vPeakClusterList.size(); ++j)
		{
			dMassDifference = vPeakClusterListCopy[i].computeAverageNeutralMonoisotopicMass()
							    - vPeakClusterList[j].computeAverageNeutralMonoisotopicMass();
			if( fabs( dMassDifference ) < ProNovoConfig::getMassAccuracyFragmentIon()/2 )
			{
				vector<MassSpectralPeak> vPeakList;
				vPeakClusterListCopy[i].getPeakList( vPeakList );
				vPeakClusterList[j].insertPeak(vPeakList);
				bMerged = true;
				break;
			}
		}
		if(!bMerged)
		{
			vPeakClusterList.push_back(vPeakClusterListCopy[i]);
		}
	}
	*/

	// sort peaks in each cluster
	for( j = 0; j < vPeakClusterList.size(); j++ )
	{
		vPeakClusterList[j].sortPeaks();
	}

	// sort clusters
	LessMassSpectralPeakCluster lessCluster("Mass");
	sort(vPeakClusterList.begin(), vPeakClusterList.end(), lessCluster);

}

void TandemMassSpectrum::extractOneCluster(vector<MassSpectralPeak *> & currentPeakList)
{
	LessMassSpectralPeak lessPeak("Intensity");
	sort(currentPeakList.begin(), currentPeakList.end(), lessPeak);
	reverse(currentPeakList.begin(), currentPeakList.end());
	vector<MassSpectralPeak *>::iterator iterMax = max_element(currentPeakList.begin(),currentPeakList.end(),lessPeak);
	MassSpectralPeak * maxPeak = *iterMax;
	currentPeakList.erase(iterMax);

	unsigned int iOptimumZ = 0;
	double dMaxScore = -1;
	unsigned int iCurrentZ;
	double dCurrentScore = 0;
	for ( iCurrentZ = 1; iCurrentZ < MAX_Z; ++iCurrentZ)
	{
		dCurrentScore = findPeaksInCluster(currentPeakList, iCurrentZ, maxPeak, false);
		if(dCurrentScore > dMaxScore)
		{
			dMaxScore = dCurrentScore;
			iOptimumZ = iCurrentZ;
		}
	}

	if(dMaxScore > 0)
	{
		findPeaksInCluster(currentPeakList, iOptimumZ, maxPeak, true);
	}
}

double TandemMassSpectrum::findPeaksInCluster(vector<MassSpectralPeak *> & currentPeakList,
		unsigned int iChargeState, MassSpectralPeak * maxPeak, bool bAddCluster)
{

	double dMassDelta = ProNovoConfig::getNeutronMass() / iChargeState;
	double dMaxPeakMZ = maxPeak->getMZ();
	double dMaxPeakIntensity = maxPeak->getIntensity();

	unsigned int i;
	unsigned int j;
	// find the peaks that are on the right side of maxPeak in the cluster
	int iPeakIndex = 0;
	double dPreviousPeakIntensity  = dMaxPeakIntensity;
	bool bFoundNextPeak = false;
	vector<unsigned int> viPeakListIndex;
	while(dPreviousPeakIntensity > dMaxPeakIntensity*0.1)
	{
		bFoundNextPeak = false;
		iPeakIndex++;
		for (i = 0; i < currentPeakList.size(); ++i)
		{
			if(currentPeakList[i]->getMZ() < dMaxPeakMZ + iPeakIndex*dMassDelta + ProNovoConfig::getMassAccuracyFragmentIon()
			&& currentPeakList[i]->getMZ() > dMaxPeakMZ + iPeakIndex*dMassDelta - ProNovoConfig::getMassAccuracyFragmentIon())
			{
				if(currentPeakList[i]->getIntensity() < dPreviousPeakIntensity*1.1
				&& currentPeakList[i]->getIntensity() > dMaxPeakIntensity*0.1 )
				{
					viPeakListIndex.push_back(i);
					dPreviousPeakIntensity = currentPeakList[i]->getIntensity();
					bFoundNextPeak = true;
					break;
				}
			}
		}
		if( !bFoundNextPeak )
		{
			dPreviousPeakIntensity = -1;
		}
	}

	// find the peaks that are on the left side of maxPeak in the cluster
	iPeakIndex = 0;
	dPreviousPeakIntensity  = dMaxPeakIntensity;
	bFoundNextPeak = false;
	while (dPreviousPeakIntensity > dMaxPeakIntensity*0.1)
	{
		bFoundNextPeak = false;
		iPeakIndex++;
		for (i = 0; i < currentPeakList.size(); ++i)
		{
			if(currentPeakList[i]->getMZ() < dMaxPeakMZ - iPeakIndex*dMassDelta + ProNovoConfig::getMassAccuracyFragmentIon()
			&& currentPeakList[i]->getMZ() > dMaxPeakMZ - iPeakIndex*dMassDelta - ProNovoConfig::getMassAccuracyFragmentIon())
			{
				if(currentPeakList[i]->getIntensity() < dPreviousPeakIntensity*1.1
				&& currentPeakList[i]->getIntensity() > dMaxPeakIntensity*0.1 )
				{
					viPeakListIndex.push_back(i);
					dPreviousPeakIntensity = currentPeakList[i]->getIntensity();
					bFoundNextPeak = true;
					break;
				}
			}
		}
		if( !bFoundNextPeak )
		{
			dPreviousPeakIntensity = -1;
		}
	}

	// get the results

	double dScore = 0;
	for (i = 0; i < viPeakListIndex.size(); ++i)
	{
		dScore = dScore + currentPeakList[viPeakListIndex[i]]->getIntensity();
	}

	if(bAddCluster && dScore > 0)
	{
		MassSpectralPeakCluster currentCluster;

		maxPeak->setChargeState(iChargeState);
		currentCluster.insertPeak(*maxPeak);
		bool bToBeAdded;
		vector<MassSpectralPeak *> currentPeakListCopy = currentPeakList;
		currentPeakList.clear();
		for (i = 0; i < currentPeakListCopy.size(); ++i)
		{
			bToBeAdded = false;
			for (j = 0; j < viPeakListIndex.size(); ++j)
			{
				if( i == viPeakListIndex[j] )
				{
					bToBeAdded = true;
					break;
				}
			}

			if(bToBeAdded)
			{
				currentPeakListCopy[i]->setChargeState(iChargeState);
				currentCluster.insertPeak(*currentPeakListCopy[i]);
			}
			else
			{
				currentPeakList.push_back(currentPeakListCopy[i]);
			}
		}
		vPeakClusterList.push_back(currentCluster);
	}

	return dScore;
}

vector<MassSpectralPeakCluster> TandemMassSpectrum::getPeakClusterList()
{
	return vPeakClusterList;
}

void TandemMassSpectrum::clearPeakClusterList()
{
	vPeakClusterList.clear();
}

void TandemMassSpectrum::clearMassSpectrum()
{
	clearPeakList();
	clearPeakClusterList();
	setScanNumber(0);
	setParentMZ(0);
	setParentChargeState(0);

}

void TandemMassSpectrum::print()
{
	cout << "iScanNumber = " << iScanNumber << endl;
	cout << "dParentMZ = " << dParentMZ << endl;
	cout << "iParentChargeState = " << iParentChargeState << endl;

	unsigned int i;

	for(i = 0; i < vPeakList.size(); i++)
	{
		cout << vPeakList[i].getMZ() << '\t'
		<< vPeakList[i].getIntensity() << '\t'
		<< vPeakList[i].getResolution() << '\t'
		<< vPeakList[i].getBaseline() << '\t'
		<< vPeakList[i].getNoise() << '\t'
		<< vPeakList[i].getChargeState() << endl;
	}

//		cout << "M/Z" << '\t'
//		<< "Num" << '\t'
//		<< "Mass" << '\t'
//		<< "Signal" << '\t'
//		<< "Resolution" << '\t'
//		<< "SNR" << endl;

	for(i = 0; i < vPeakClusterList.size(); i++)
	{
//		vPeakClusterList[i].print();
//		cout << vPeakClusterList[i].getMonoisotopicPeak().getMZ() << '\t'
//				<< vPeakClusterList[i].getPeakCount() << '\t'
//		<< vPeakClusterList[i].computeAverageNeutralMonoisotopicMass() << '\t'
//		<< vPeakClusterList[i].computeTotalSignalIntensity() << '\t'
//		<< vPeakClusterList[i].computeAverageResolution() << '\t'
//		<< vPeakClusterList[i].computeAggregateSignalToNoiseRatio() << endl;
	}

}
