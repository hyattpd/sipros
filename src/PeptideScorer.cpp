#include "PeptideScorer.h"

ProductIon::ProductIon()
{
	cIonType = 'x';
	iIonNumber = 0;
	iCharge = 1;
	dMostAbundantMass = 0;
	dMostAbundantMZ = 0;
	dMZError = 0;
	dMassError = 0;
	dScoreWeight = 1.0;
	bComplementaryFragmentObserved = false;
}

ProductIon::~ProductIon()
{

}

void ProductIon::setProductIon(char cIonTypeInput, int iIonNumberInput, int iChargeInput)
{
	cIonType = cIonTypeInput;
	iIonNumber = iIonNumberInput;
	iCharge = iChargeInput;
}

void ProductIon::setObservedInfo(double dMZErrorInput, double dWeightInput,
		double dMostAbundantMZInput, int iMostAbundantPeakIndexInput)
{
	dMZError = dMZErrorInput;
	dMassError = dMZError * iCharge;
	dScoreWeight = dWeightInput;
	dMostAbundantMZ = dMostAbundantMZInput;
	double dProtonMass = 1.007825;
	dMostAbundantMass = dMostAbundantMZ * iCharge - dProtonMass * iCharge;
	iMostAbundantPeakIndex = iMostAbundantPeakIndexInput;
}

void ProductIon::setComplementaryFragmentObserved(bool bComplementaryFragmentObservedInput)
{
	bComplementaryFragmentObserved = bComplementaryFragmentObservedInput;
}


bool LessProductIon::operator() (ProductIon ion1, ProductIon ion2 ) const
{
	if( sKey == "mass" )
	{
		if( ion1.getMostAbundantMass() < ion2.getMostAbundantMass() )
			return true;
		else
			return false;
	}
	else if ( sKey == "ion" )
	{
		// sort by ion series
		if( ion1.getIonType() < ion2.getIonType() )
			return true;
		else if( ion1.getIonType() > ion2.getIonType() )
			return false;
		else
		{
			if(ion1.getIonNumber() < ion2.getIonNumber())
				return true;
			else
				return false;
		}
	}
	else
	{
		if( ion1.getMostAbundantMass() < ion2.getMostAbundantMass() )
			return true;
		else
			return false;
	}
}

PeptideScorer::PeptideScorer()
{
	dMassTolerance = ProNovoConfig::getMassAccuracyFragmentIon();

	dNtermMass = ProNovoConfig::getTerminusMassN();
	dCtermMass = ProNovoConfig::getTerminusMassC();
	dProtonMass = 1.007825;

	vector<string> vsSingleResidueNames = ProNovoConfig::vsSingleResidueNames;
	vector<double> vdSingleResidueMasses = ProNovoConfig::vdSingleResidueMasses;

	for (unsigned int n = 0; n < vsSingleResidueNames.size(); ++n)
	{
		mapResidueMass[ vsSingleResidueNames[n][0] ] = vdSingleResidueMasses[n];
	}

}

PeptideScorer::~PeptideScorer()
{

}

bool PeptideScorer::setMS2(
			const vector<double> & vdMZinput,
			const vector<double> & vdIntensityInput,
			const vector<int> & viZinputTemp,
			const double & dParentMZinput,
			const int & iParentChargeStateInput)
{
	int i;
	// assign MS2 data
	vdMZ.clear();
	vdIntensity.clear();
	viZinput.clear();

	vdMZ = vdMZinput;
	vdIntensity = vdIntensityInput;
	viZinput = viZinputTemp;
	dParentMZ = dParentMZinput;
	iParentChargeState = iParentChargeStateInput;

	if(vdMZinput.size() != vdIntensityInput.size() || vdMZinput.size() != viZinput.size())
	{
		cout << "ERROR: Problem with the input mass spectrum" << endl;
		return false;
	}

	// bubble sort the peak list by MZ in ascending order
	int iPeakCount = (int)vdMZ.size();

	int pass;
	double dCurrentMZ;
	double dCurrentInt;
	int iCurrentZ;

	for (pass=1; pass < iPeakCount; pass++) {  // count how many times
		// This next loop becomes shorter and shorter
       for (i=0; i < iPeakCount-pass; i++) {
           if (vdMZ[i] > vdMZ[i+1]) {
               // exchange
               dCurrentMZ = vdMZ[i];
               dCurrentInt = vdIntensity[i];
               iCurrentZ = viZinput[i];

               vdMZ[i] = vdMZ[i+1];
               vdIntensity[i] = vdIntensity[i+1];
               viZinput[i] = viZinput[i+1];

               vdMZ[i+1] = dCurrentMZ;
               vdIntensity[i+1] = dCurrentInt;
               viZinput[i+1] = iCurrentZ;
           }
       }
	}

	// change intensity to relative intensity
	double dMaxInt = vdIntensity[getMaxValueIndex(vdIntensity)];
	for (i = 0; i < (int)vdIntensity.size(); ++i)
	{
		vdIntensity[i] = (vdIntensity[i]/dMaxInt)*100;
	}

	// compute normalized intensity
	// vdMzIntensity[i] is the intensity at M/Z of i (rounded)
	vector<double> vdMzIntensity(5000);
	fill(vdMzIntensity.begin(), vdMzIntensity.end(), 0.0);
	for (i = 0; i < iPeakCount; ++i)
	{
		// if there are multiple peaks in this M/Z bin, use the intensity of the highest peak
		if(vdIntensity[i] > vdMzIntensity[(int)vdMZ[i]] )
		{
			vdMzIntensity[(int)vdMZ[i]] = vdIntensity[i];
		}
	}

	// vdMzIntensity[i] is the maximum intensity at M/Z window of i plus and minus iMzRange
	vector<double> vdMaxMzIntensity(5000);
	fill(vdMaxMzIntensity.begin(), vdMaxMzIntensity.end(), 1.0);
	int iMzRange = 100;
	int iLowerBound = 0;
	int iUpperBound = 0;
	for (i = 0; i < (int)vdMaxMzIntensity.size(); ++i)
	{
		iLowerBound = max(0, i - iMzRange);
		iUpperBound = min(i + iMzRange, (int)vdMzIntensity.size());
		vdMaxMzIntensity[i] = *max_element(vdMzIntensity.begin() + iLowerBound, vdMzIntensity.begin() + iUpperBound );
	}

	vdNormalizedIntensity.clear();
	for (i = 0; i < (int)vdIntensity.size(); ++i)
	{
		vdNormalizedIntensity.push_back(vdIntensity[i]/vdMaxMzIntensity[(int)vdMZ[i]]);
//		cout << vdMZ[i] << "		" << vdIntensity[i] << "		" << vdNormalizedIntensity[i] << endl;
	}

	// populate vbPeakPresenceBins

	if(vdMZ.front() <= 0)
	{
		cout << "ERROR: negative MZ value = " << vdMZ.front() << endl;
		return false;
	}

	vbPeakPresenceBins.clear();
// compute the number of vbPeakPresenceBins needed, add 10 to make sure there is a bin for the last peak
	unsigned long int iBinNumber = (unsigned long int)(( vdMZ.back() + 10.0 ) * BIN_RES);

	// populate all vbPeakPresenceBins with false
//	vbPeakPresenceBins.resize(iBinNumber, false);
	vbPeakPresenceBins.resize(iBinNumber, -1);

	unsigned long int iTarget;
	unsigned long int k;

	unsigned long int iBinRange = (unsigned long int)(dMassTolerance * BIN_RES) + 1;
	for(i = 0; i < (int)vdMZ.size(); i++)
	{
		iTarget = (unsigned long int)(vdMZ[i]*BIN_RES);

		for(k = iTarget-iBinRange; k <= iTarget+iBinRange; k++)
		{
			if(vbPeakPresenceBins[k] == -1)
			{
				vbPeakPresenceBins[k] = i;
			}
			else
			{
				if( vdNormalizedIntensity[vbPeakPresenceBins[k]] < vdNormalizedIntensity[i] )
				{
					vbPeakPresenceBins[k] = i;
				}
			}
		//	vbPeakPresenceBins[k] = true;
		}

//		fill(vbPeakPresenceBins.begin()+iTarget-iBinRange, vbPeakPresenceBins.begin()+iTarget+iBinRange+1, 'T');
	}

	// process MS2 for prelimary scoring
	vdMonoMass.clear();
	vdMonoInt.clear();
	vdMonoZ.clear();


	TandemMassSpectrum currentSpectrum;
	currentSpectrum.setParentMZ(dParentMZ);
	currentSpectrum.setParentChargeState(iParentChargeState);
	for (i = 0; i < (int)vdMZ.size(); ++i)
	{
//		currentSpectrum.addPeak(vdMZ[i], vdIntensity[i], 0, 0, 0, viZinput[i]);
		currentSpectrum.addPeak(vdMZ[i], vdNormalizedIntensity[i], 0, 0, 0, viZinput[i]);
	}
	currentSpectrum.clusterPeaks();

	// vClusterList is sorted by their mass in ascending order
	vector< MassSpectralPeakCluster > vClusterList = currentSpectrum.getPeakClusterList();

	for (i = 0; i < (int)vClusterList.size(); ++i)
	{
		vdMonoMass.push_back(vClusterList[i].computeAverageNeutralMonoisotopicMass());
		vdMonoInt.push_back(vClusterList[i].computeTotalSignalIntensity());
		vdMonoZ.push_back(vClusterList[i].getChargeState());

	//	cout << vdMonoMass[i] << "\t" << vdMonoInt[i] << "\t" << vdMonoZ[i] << endl;
	//	vClusterList[i].print();
	}


	// processing MS2 for primary scoring



	return true;

}

double PeptideScorer::prelimaryScorePeptide3( const string & sPeptide)
{
	vector<double> vdYmass;
	vector<double> vdBmass;
	double dPeptideMass;
	int iNumFragments = 0;
	double dCurrentScore;
	double dYweight;
	double dBweight;
	int n;

	// compute b and y ion masses
	computeYBionMass(sPeptide, vdYmass, vdBmass, dPeptideMass);
	iNumFragments = vdYmass.size();
	dYweight = 0;
	dBweight = 0;
	int iIndex4MostAbundunt = 0;
	double dTargetMZ;
	double dFragmentScore = 0;
	double dBonus;
	dCurrentScore = 0;
	for (n = 0; n < iNumFragments; ++n)
	{
		dYweight = 0;
		dBweight = 0;

		dTargetMZ = vdYmass[n] + ProNovoConfig::getProtonMass();
		if( searchMZ( dTargetMZ, iIndex4MostAbundunt) )
		{
			dFragmentScore = ProNovoConfig::scoreError(fabs(vdMZ[iIndex4MostAbundunt]-dTargetMZ) ) + vdNormalizedIntensity[iIndex4MostAbundunt];
			dBonus = dFragmentScore;
			if(viZinput[iIndex4MostAbundunt] == 1)
			{
				dFragmentScore += dBonus;
			}
			dYweight += dFragmentScore;

			// search for H2O loss and NH3 loss, only one will be considered, because they are off by only 1 Da
			if(searchMZ(dTargetMZ-ProNovoConfig::getMassH2O(), iIndex4MostAbundunt))
			{
				dFragmentScore += dBonus * 0.2;
			}
			else if(searchMZ(dTargetMZ-ProNovoConfig::getMassNH3(), iIndex4MostAbundunt))
			{
				dFragmentScore += dBonus * 0.2;
			}

			dYweight += dFragmentScore;
		}
		dTargetMZ = vdBmass[iNumFragments - n - 1]+ ProNovoConfig::getProtonMass();
		if( searchMZ( dTargetMZ, iIndex4MostAbundunt) )
		{
			dFragmentScore = ProNovoConfig::scoreError(fabs(vdMZ[iIndex4MostAbundunt]-dTargetMZ) ) + vdNormalizedIntensity[iIndex4MostAbundunt];
			dBonus = dFragmentScore;
			if(viZinput[iIndex4MostAbundunt] == 1)
			{
				dFragmentScore += dBonus;
			}

			if(searchMZ(dTargetMZ-ProNovoConfig::getMassH2O(), iIndex4MostAbundunt))
			{
				dFragmentScore += dBonus * 0.2;
			}
			else if(searchMZ(dTargetMZ-ProNovoConfig::getMassNH3(), iIndex4MostAbundunt))
			{
				dFragmentScore += dBonus * 0.2;
			}

			dBweight += dFragmentScore;
		}

		if( iParentChargeState > 1 )
		{

			dTargetMZ = vdYmass[n]/2.0 + ProNovoConfig::getProtonMass();
			if( searchMZ( dTargetMZ, iIndex4MostAbundunt) )
			{
				dFragmentScore = ProNovoConfig::scoreError(fabs(vdMZ[iIndex4MostAbundunt]-dTargetMZ) ) + vdNormalizedIntensity[iIndex4MostAbundunt];
				dBonus = dFragmentScore;

				if(viZinput[iIndex4MostAbundunt] == 2)
				{
					dFragmentScore += dBonus;
				}

				if(searchMZ(dTargetMZ-ProNovoConfig::getMassH2O()/2.0, iIndex4MostAbundunt))
				{
					dFragmentScore += dBonus * 0.2;
				}
				else if(searchMZ(dTargetMZ-ProNovoConfig::getMassNH3()/2.0, iIndex4MostAbundunt))
				{
					dFragmentScore += dBonus * 0.2;
				}

				dYweight += dFragmentScore;
			}
			dTargetMZ = vdBmass[iNumFragments - n - 1]/2.0 + ProNovoConfig::getProtonMass();
			if( searchMZ( dTargetMZ, iIndex4MostAbundunt) )
			{
				dFragmentScore = ProNovoConfig::scoreError(fabs(vdMZ[iIndex4MostAbundunt]-dTargetMZ) ) + vdNormalizedIntensity[iIndex4MostAbundunt];
				dBonus = dFragmentScore;

				if(viZinput[iIndex4MostAbundunt] == 2)
				{
					dFragmentScore += dBonus;
				}

				// search for H2O loss and NH3 loss, only one will be considered, because they are off by only 1 Da
				if(searchMZ(dTargetMZ-ProNovoConfig::getMassH2O()/2.0, iIndex4MostAbundunt))
				{
					dFragmentScore += dBonus * 0.2;
				}
				else if(searchMZ(dTargetMZ-ProNovoConfig::getMassNH3()/2.0, iIndex4MostAbundunt))
				{
					dFragmentScore += dBonus * 0.2;
				}

				dBweight += dFragmentScore;
			}

		}


		if(dYweight > 0.001 && dBweight > 0.001 )
		{
			// bonus point for finding complimentary B and Y ions
			dCurrentScore = dCurrentScore + (dYweight + dBweight)*2;
		}
		else
		{
			dCurrentScore = dCurrentScore + (dYweight + dBweight);
		}



	}

	return dCurrentScore;
}

double PeptideScorer::prelimaryScorePeptide2( const string & sPeptide)
{
	vector<double> vdYmass;
	vector<double> vdBmass;
	double dPeptideMass;
	int iNumFragments = 0;
	double dCurrentScore;
	double dYweight;
	double dBweight;
	int n;

	// compute b and y ion masses
	computeYBionMass(sPeptide, vdYmass, vdBmass, dPeptideMass);
	iNumFragments = vdYmass.size();

	dCurrentScore = 0;
	dYweight = 0;
	dBweight = 0;
	vector<int> viIndex4Found;
	unsigned int i;
	int iIndex4MostAbundunt = 0;
	double dCurrentAbundance = 0;
	for (n = 0; n < iNumFragments; ++n)
	{

		if( binarySearch( vdYmass[n], vdMonoMass, dMassTolerance, viIndex4Found) )
		{
			dCurrentAbundance = 0;
			for (i = 0; i < viIndex4Found.size(); ++i)
			{
				if(vdMonoInt[viIndex4Found[i]] > dCurrentAbundance)
				{
					dCurrentAbundance = vdMonoInt[viIndex4Found[i]];
					iIndex4MostAbundunt = viIndex4Found[i];
				}
			}
			dYweight = 2*(ProNovoConfig::scoreError(fabs(vdMonoMass[iIndex4MostAbundunt]-vdYmass[n]) ) + dCurrentAbundance);
		}
		else if (searchMZ( vdYmass[n] + ProNovoConfig::getProtonMass(), iIndex4MostAbundunt) )
		{
			dYweight = ProNovoConfig::scoreError(fabs(vdMZ[iIndex4MostAbundunt]-vdYmass[n]-ProNovoConfig::getProtonMass()) ) + vdNormalizedIntensity[iIndex4MostAbundunt];
		}
		else
		{
			dYweight = 0;
		}

		if( binarySearch( vdBmass[iNumFragments - n - 1], vdMonoMass, dMassTolerance, viIndex4Found) )
		{
			dCurrentAbundance = 0;
			for (i = 0; i < viIndex4Found.size(); ++i)
			{
				if(vdMonoInt[viIndex4Found[i]] > dCurrentAbundance)
				{
					dCurrentAbundance = vdMonoInt[viIndex4Found[i]];
					iIndex4MostAbundunt = viIndex4Found[i];
				}
			}
			dBweight = 2 * ( ProNovoConfig::scoreError(fabs(vdMonoMass[iIndex4MostAbundunt]-vdBmass[iNumFragments - n - 1])) + dCurrentAbundance );
		}
		else if ( searchMZ( vdBmass[iNumFragments - n - 1]+ ProNovoConfig::getProtonMass(), iIndex4MostAbundunt) )
		{
			dBweight = ProNovoConfig::scoreError(fabs(vdMZ[iIndex4MostAbundunt]-vdBmass[iNumFragments - n - 1]-ProNovoConfig::getProtonMass())) + vdNormalizedIntensity[iIndex4MostAbundunt];
		}
		else
		{
			dBweight = 0;
		}

		if(dYweight > 0 && dBweight > 0)
		{
			// bonus point for finding complimentary B and Y ions
			dCurrentScore = dCurrentScore + (dYweight + dBweight)*2;
		}
		else
		{
			dCurrentScore = dCurrentScore + (dYweight + dBweight);
		}
	}

	return dCurrentScore;
}

double PeptideScorer::prelimaryScorePeptide( const string & sPeptide)
{
	vector<double> vdYmass;
	vector<double> vdBmass;
	double dPeptideMass;
	int iNumFragments = 0;
	double dCurrentScore;
	double dYweight;
	double dBweight;
	int n;

	// compute b and y ion masses
	computeYBionMass(sPeptide, vdYmass, vdBmass, dPeptideMass);
	iNumFragments = vdYmass.size();

	dCurrentScore = 0;
	dYweight = 0;
	dBweight = 0;
	vector<int> viIndex4Found;
	unsigned int i;
	int iIndex4MostAbundunt = 0;
	double dCurrentAbundance = 0;
	for (n = 0; n < iNumFragments; ++n)
	{

		if( binarySearch( vdYmass[n], vdMonoMass, dMassTolerance, viIndex4Found) )
		{
			dCurrentAbundance = 0;
			for (i = 0; i < viIndex4Found.size(); ++i)
			{
				if(vdMonoInt[viIndex4Found[i]] > dCurrentAbundance)
				{
					dCurrentAbundance = vdMonoInt[viIndex4Found[i]];
					iIndex4MostAbundunt = viIndex4Found[i];
				}
			}
			dYweight = ProNovoConfig::scoreError(fabs(vdMonoMass[iIndex4MostAbundunt]-vdYmass[n]));
		}
		else
		{
			dYweight = 0;
		}

		if( binarySearch( vdBmass[iNumFragments - n - 1], vdMonoMass, dMassTolerance, viIndex4Found) )
		{
			dCurrentAbundance = 0;
			for (i = 0; i < viIndex4Found.size(); ++i)
			{
				if(vdMonoInt[viIndex4Found[i]] > dCurrentAbundance)
				{
					dCurrentAbundance = vdMonoInt[viIndex4Found[i]];
					iIndex4MostAbundunt = viIndex4Found[i];
				}
			}
			dBweight = ProNovoConfig::scoreError(fabs(vdMonoMass[iIndex4MostAbundunt]-vdBmass[iNumFragments - n - 1]));		}
		else
		{
			dBweight = 0;
		}

		if(dYweight > 0 && dBweight > 0)
		{
			// bonus point for finding complimentary B and Y ions
			dCurrentScore = dCurrentScore + (dYweight + dBweight)*2;
		}
		else
		{
			dCurrentScore = dCurrentScore + (dYweight + dBweight);
		}
	}

	return dCurrentScore;
}

bool PeptideScorer::searchMZ(const double & dTarget, int & iIndex4Found )
{

	if(dTarget > vdMZ.back())
		return false;
	iIndex4Found = vbPeakPresenceBins[(unsigned long int)(dTarget*BIN_RES)];
	if( iIndex4Found == -1)
		return false;
	else
		return true;


////	if(!vbPeakPresenceBins[(unsigned long int)(dTarget*BIN_RES)])////		return false;


//	if(!(vbPeakPresenceBins[targetBin] || vbPeakPresenceBins[targetBin-1] || vbPeakPresenceBins[targetBin+1]))
//	{
//		return false;
//	}
	/*/*
	// real binary search
	int low = 0;
	int high = vdMZ.size() - 1;
	int mid;
	double dCurrentMaxAbundance = 0;
	while(low <= high)
	{
		mid = (low + high) / 2;
		if( vdMZ[mid] > dTarget + dMassTolerance )
		{
			high = mid - 1;
		}
		else if (vdMZ[mid] < dTarget - dMassTolerance )
		{
			low = mid + 1;
		}
		else
		{
			// found vdMZ[mid] is within the dMassTolerance range of dTarget
			iIndex4Found = mid;
			dCurrentMaxAbundance = vdNormalizedIntensity[mid];

			// now find all other values that is within the range
			int upperBound = mid + 1;
			int lowerBound = mid - 1;

			while(upperBound < (int)vdMZ.size() && vdMZ[upperBound] < dTarget + dMassTolerance)
			{
				if(vdNormalizedIntensity[upperBound] > dCurrentMaxAbundance)
				{
					iIndex4Found = upperBound;
					dCurrentMaxAbundance = vdNormalizedIntensity[upperBound];
				}
				upperBound += 1;
			}

			while(lowerBound >= 0 && vdMZ[lowerBound] > dTarget - dMassTolerance)
			{
				if(vdNormalizedIntensity[lowerBound] > dCurrentMaxAbundance)
				{
					iIndex4Found = lowerBound;
					dCurrentMaxAbundance = vdNormalizedIntensity[lowerBound];
				}
				lowerBound -= 1;
			}

			return true;
		}
	}
	return fals
	*/
}

bool PeptideScorer::binarySearch(const double & dTarget,
		const vector<double> & vdList, const double & dTolerance,
		vector<int> & viIndex4Found )
{
	viIndex4Found.clear();
	int low = 0;
	int high = vdList.size() - 1;
	int mid;
	while(low <= high)
	{
		mid = (low + high) / 2;
		if( vdList[mid] > dTarget + dTolerance )
		{
			high = mid - 1;
		}
		else if (vdList[mid] < dTarget - dTolerance )
		{
			low = mid + 1;
		}
		else
		{
			// found vdList[mid] is within the dTolerance range of dTarget
			viIndex4Found.push_back(mid);

			// now find all other values that is within the range
			int upperBound = mid + 1;
			int lowerBound = mid - 1;

			while(upperBound < (int)vdList.size() && vdList[upperBound] < dTarget + dTolerance)
			{
				viIndex4Found.push_back(upperBound);
				upperBound += 1;
			}

			while(lowerBound >= 0 && vdList[lowerBound] > dTarget - dTolerance)
			{
				viIndex4Found.push_back(lowerBound);
				lowerBound -= 1;
			}

			sort(viIndex4Found.begin(), viIndex4Found.end());
			return true;
		}	}
	return false;
}




bool PeptideScorer::computeYBionMass(const string & sPeptide, vector<double> & vdYmass,
				vector<double> & vdBmass, double & dPeptideMass)
{
	map<char,double>::iterator iter;
	dPeptideMass = 0;
	vdYmass.clear();
	vdYmass.reserve(sPeptide.length());
	vdBmass.clear();
	vdBmass.reserve(sPeptide.length());
	int i;
	if(sPeptide[0] != '[')
	{
		cout << "ERROR: peptide sequence must start with '[' as N terminus; Invalid sequence = " << sPeptide << endl;
	}
	for (i = 1; i < (int)sPeptide.length(); ++i)
	{
		if(sPeptide[i] == ']')
		{
			// hit the C terminus and break out of the loop
			break;
		}
		iter = mapResidueMass.find(sPeptide[i]);
		if(iter == mapResidueMass.end() )
		{
			cout << "WARNING: Residue " << sPeptide[i] << " Peptide " << sPeptide << " is not defined in the config." << endl;
		}
		else
		{
			if(isalpha(sPeptide[i]))
			{
				// this is an amino acid residue
				dPeptideMass = dPeptideMass + iter->second;
				vdBmass.push_back(dPeptideMass);
			}
			else
			{
				// this is a PTM and its mass should be added to the proceding residue
			//	cout << "PTM " << sPeptide[i] << " # " << iter->second << endl;
				dPeptideMass = dPeptideMass + iter->second;
				if(vdBmass.size()>0)
				{
					vdBmass.back() += iter->second;
				}
				else
				{
					// this symbol is represents a PTM to the N terminus, whose mass will be added to the next residue.
				}
			}

		}
	}
	for (i = i+1; i < (int)sPeptide.length(); ++i)
	{
		if(!isalpha(sPeptide[i]))
		{
			// this is PTM to the C terminus
			iter = mapResidueMass.find(sPeptide[i]);
			if(iter == mapResidueMass.end() )
			{
				cout << "WARNING: Residue " << sPeptide[i] << " Peptide " << sPeptide << " is not defined in the config." << endl;
			}
			else
			{
				dPeptideMass = dPeptideMass + iter->second;
			}
		}
	}
	dPeptideMass = dPeptideMass + dCtermMass + dNtermMass;
	vdBmass.pop_back();

	for (i = vdBmass.size() - 1; i >= 0; --i)
	{
		vdYmass.push_back(dPeptideMass - vdBmass[i]);
	}
	/*
	for (i = 0; i < (int)vdBmass.size(); ++i)
	{
		cout << " simple b and y " << i << " =\t" << vdBmass[i] << "\t"<< vdYmass[i]<< endl;
	}
	*/
	return true;

}

double PeptideScorer::primaryScorePetide( const string & sPeptide)
{
	int iPeptideLength = 0;
	for (int i = 0; i < (int)sPeptide.length(); ++i)
	{
		if(isalpha(sPeptide[i]))
		{
			iPeptideLength = iPeptideLength + 1;
		}
	}
	vector< vector<double> > vvdYionMass;
	vector< vector<double> > vvdYionProb;
	vector< vector<double> > vvdBionMass;
	vector< vector<double> > vvdBionProb;
	ProNovoConfig::configIsotopologue.computeProductIon( sPeptide,
			vvdYionMass, vvdYionProb,
			vvdBionMass, vvdBionProb);

	/*
	// print out b and y ion masses
	for (int m = 0; m < (int)vvdBionMass.size(); ++m)
	{
		cout << "b " << m << " = " << vvdBionMass[m][0] << "	y " << m << " = " << vvdYionMass[m][0] << endl;
	}
	*/

	int n; // Ion number starting from one
	int z; // charge state
	vector<ProductIon> vFoundIons;
	double dScoreWeight = 0;
	double dMZError = 1;
	double dMostAbundantObservedMZ = 0;
	int iMostAbundantPeakIndex = 0;
	for (n = 0; n < (int)vvdYionMass.size(); ++n)
	{
		for (z = 1; z <= iParentChargeState; ++z)
		{
			ProductIon currentIon;
			currentIon.setProductIon('y', n+1, z);
			if(findProductIon(vvdYionMass[n], vvdYionProb[n], z,
					dScoreWeight, dMZError, dMostAbundantObservedMZ, iMostAbundantPeakIndex))
			{
				currentIon.setObservedInfo(dMZError, dScoreWeight, dMostAbundantObservedMZ, iMostAbundantPeakIndex);
				vFoundIons.push_back(currentIon);
			}
		}
	}
	for (n = 0; n < (int)vvdBionMass.size(); ++n)
	{
		for (z = 1; z <= iParentChargeState; ++z)
		{
			ProductIon currentIon;
			currentIon.setProductIon('b', n+1, z);
			if(findProductIon(vvdBionMass[n], vvdBionProb[n], z,
					dScoreWeight, dMZError, dMostAbundantObservedMZ, iMostAbundantPeakIndex))
			{
				currentIon.setObservedInfo(dMZError, dScoreWeight, dMostAbundantObservedMZ, iMostAbundantPeakIndex);
				vFoundIons.push_back(currentIon);
			}
		}
	}
	int i;
	int j;
	for (i = 0; i < (int)vFoundIons.size(); ++i)
	{
		vFoundIons[i].setComplementaryFragmentObserved(false);
	}
	for (i = 0; i < (int)vFoundIons.size(); ++i)
	{
		for (j = i+1; j < (int)vFoundIons.size(); ++j)
		{
			if(vFoundIons[i].getIonNumber() + vFoundIons[j].getIonNumber() == iPeptideLength)
			{
				if((vFoundIons[i].getIonType() == 'y' && vFoundIons[j].getIonType() == 'b')
				 ||(vFoundIons[i].getIonType() == 'b' && vFoundIons[j].getIonType() == 'y'))
				{
					vFoundIons[i].setComplementaryFragmentObserved(true);
					vFoundIons[j].setComplementaryFragmentObserved(true);
				}
			}

		}
	}

//	print out found product ions

	/*
//	LessProductIon lessIon("mass");
//	sort(vFoundIons.begin(), vFoundIons.end(), lessIon);
	cout << "Ion	Charge	MZError	Weight	MZ	Int" << endl;
	for (i = 0; i < (int)vFoundIons.size(); ++i)
	{

		cout << vFoundIons[i].getIonType() << vFoundIons[i].getIonNumber() << "\t"
		<< vFoundIons[i].getCharge() << "\t" <<
//		(vFoundIons[i].getMZError()/vFoundIons[i].getMostAbundantMZ())*1000000 <<
		vFoundIons[i].getMZError() <<
		"\t" << vFoundIons[i].getScoreWeight() << "\t"
		<< vFoundIons[i].getMostAbundantMZ() << "\t"
		<< vFoundIons[i].getMostAbundantInt() << endl;
	}
	*/


	double dAverageMZError = 0;
	for (i = 0; i < (int)vFoundIons.size(); ++i)
	{
		dAverageMZError += vFoundIons[i].getMZError();
	}
	dAverageMZError = dAverageMZError / (double)vFoundIons.size();

	double dPrimaryScore = 0;
	double dBonus4ComplementaryFragmentObserved = 1.0;
	for (i = 0; i < (int)vFoundIons.size(); ++i)
	{
		if(vFoundIons[i].getComplementaryFragmentObserved())
		{
			dBonus4ComplementaryFragmentObserved = 2.0;
		}
		else
		{
			dBonus4ComplementaryFragmentObserved = 1.0;
		}

		dPrimaryScore += ProNovoConfig::scoreError(fabs(vFoundIons[i].getMZError()-dAverageMZError))*vFoundIons[i].getScoreWeight()*dBonus4ComplementaryFragmentObserved;
	}

	return dPrimaryScore;
}

double PeptideScorer::primaryScorePetide2( const string & sPeptide)
{
	int iPeptideLength = 0;
	for (int i = 0; i < (int)sPeptide.length(); ++i)
	{
		if(isalpha(sPeptide[i]))
		{
			iPeptideLength = iPeptideLength + 1;
		}
	}
	vector< vector<double> > vvdYionMass;
	vector< vector<double> > vvdYionProb;
	vector< vector<double> > vvdBionMass;
	vector< vector<double> > vvdBionProb;
	ProNovoConfig::configIsotopologue.computeProductIon( sPeptide,
			vvdYionMass, vvdYionProb,
			vvdBionMass, vvdBionProb);
/*
	for (int m = 0; m < (int)vvdBionMass.size(); ++m)
	{
		cout << "b " << m << " = " << vvdBionMass[m][0] << endl;
	}
*/

	int n; // Ion number starting from one
	int z; // charge state
	vector<ProductIon> vFoundIons;
	double dScoreWeight = 0;
	double dMZError = 1;
	double dMostAbundantObservedMZ = 0;
	int iMostAbundantPeakIndex = 0;
	for (n = 0; n < (int)vvdYionMass.size(); ++n)
	{
		for (z = 1; z <= iParentChargeState; ++z)
		{
			ProductIon currentIon;
			currentIon.setProductIon('y', n+1, z);
			if(findProductIon(vvdYionMass[n], vvdYionProb[n], z,
					dScoreWeight, dMZError, dMostAbundantObservedMZ, iMostAbundantPeakIndex))
			{
				currentIon.setObservedInfo(dMZError, dScoreWeight, dMostAbundantObservedMZ, iMostAbundantPeakIndex);
				vFoundIons.push_back(currentIon);
			}
		}
	}
	for (n = 0; n < (int)vvdBionMass.size(); ++n)
	{
		for (z = 1; z <= iParentChargeState; ++z)
		{
			ProductIon currentIon;
			currentIon.setProductIon('b', n+1, z);
			if(findProductIon(vvdBionMass[n], vvdBionProb[n], z,
					dScoreWeight, dMZError, dMostAbundantObservedMZ, iMostAbundantPeakIndex))
			{
				currentIon.setObservedInfo(dMZError, dScoreWeight, dMostAbundantObservedMZ, iMostAbundantPeakIndex);
				vFoundIons.push_back(currentIon);
			}
		}
	}
	int i;
	int j;
	for (i = 0; i < (int)vFoundIons.size(); ++i)
	{
		vFoundIons[i].setComplementaryFragmentObserved(false);
	}
	for (i = 0; i < (int)vFoundIons.size(); ++i)
	{
		for (j = i+1; j < (int)vFoundIons.size(); ++j)
		{
			if(vFoundIons[i].getIonNumber() + vFoundIons[j].getIonNumber() == iPeptideLength)
			{
				if((vFoundIons[i].getIonType() == 'y' && vFoundIons[j].getIonType() == 'b')
				 ||(vFoundIons[i].getIonType() == 'b' && vFoundIons[j].getIonType() == 'y'))
				{
					vFoundIons[i].setComplementaryFragmentObserved(true);
					vFoundIons[j].setComplementaryFragmentObserved(true);
				}
			}

		}
	}

//	print out found product ions

	/*
//	LessProductIon lessIon("mass");
//	sort(vFoundIons.begin(), vFoundIons.end(), lessIon);
	cout << "Ion	Charge	MZError	Weight	MZ	Int" << endl;
	for (i = 0; i < (int)vFoundIons.size(); ++i)
	{

		cout << vFoundIons[i].getIonType() << vFoundIons[i].getIonNumber() << "\t"
		<< vFoundIons[i].getCharge() << "\t" <<
//		(vFoundIons[i].getMZError()/vFoundIons[i].getMostAbundantMZ())*1000000 <<
		vFoundIons[i].getMZError() <<
		"\t" << vFoundIons[i].getScoreWeight() << "\t"
		<< vFoundIons[i].getMostAbundantMZ() << "\t"
		<< vFoundIons[i].getMostAbundantInt() << endl;
	}
	*/


	double dAverageMZError = 0;
	for (i = 0; i < (int)vFoundIons.size(); ++i)
	{
		dAverageMZError += vFoundIons[i].getMZError();
	}
	dAverageMZError = dAverageMZError / (double)vFoundIons.size();

	double dPrimaryScore = 0;
	double dIonScore = 0;
	double dBonus4ComplementaryFragmentObserved = 1.0;
	for (i = 0; i < (int)vFoundIons.size(); ++i)
	{
		if(vFoundIons[i].getComplementaryFragmentObserved())
		{
			dBonus4ComplementaryFragmentObserved = 2.0;
		}
		else
		{
			dBonus4ComplementaryFragmentObserved = 1.0;
		}


		dIonScore = ProNovoConfig::scoreError(fabs(vFoundIons[i].getMZError()-dAverageMZError)) + vdNormalizedIntensity[ vFoundIons[i].getMostAbundantPeakIndex() ];
		dPrimaryScore += dIonScore * vFoundIons[i].getScoreWeight() * dBonus4ComplementaryFragmentObserved;
	}

	return dPrimaryScore;
}

double PeptideScorer::primaryScorePetide3( const string & sPeptide)
{
	int iPeptideLength = 0;
	for (int i = 0; i < (int)sPeptide.length(); ++i)
	{
		if(isalpha(sPeptide[i]))
		{
			iPeptideLength = iPeptideLength + 1;
		}
	}
	vector< vector<double> > vvdYionMass;
	vector< vector<double> > vvdYionProb;
	vector< vector<double> > vvdBionMass;
	vector< vector<double> > vvdBionProb;

	vector<double> vdYmass;
	vector<double> vdBmass;
	// compute b and y ion masses
	double dPeptideMass;
	computeYBionMass(sPeptide, vdYmass, vdBmass, dPeptideMass);
	vector<double> vdTempMass;
	vector<double> vdTempProb;
	for (unsigned int k = 0; k < vdYmass.size(); ++k)
	{
		ProNovoConfig::configIsotopologue.getPolyAveragine( k+1, vdYmass[k],
				vdTempMass, vdTempProb);
		vvdYionMass.push_back(vdTempMass);
		vvdYionProb.push_back(vdTempProb);
	}
	for (unsigned int k = 0; k < vdBmass.size(); ++k)
	{
		ProNovoConfig::configIsotopologue.getPolyAveragine( k+1, vdBmass[k],
				vdTempMass, vdTempProb);
		vvdBionMass.push_back(vdTempMass);
		vvdBionProb.push_back(vdTempProb);
	}
//	ProNovoConfig::configIsotopologue.computeProductIon( sPeptide,
//			vvdYionMass, vvdYionProb,
//			vvdBionMass, vvdBionProb);
/*
	for (int m = 0; m < (int)vvdBionMass.size(); ++m)
	{
		cout << "b " << m << " = " << vvdBionMass[m][0] << endl;
	}
*/

	int n; // Ion number starting from one
	int z; // charge state
	vector<ProductIon> vFoundIons;
	double dScoreWeight = 0;
	double dMZError = 1;
	double dMostAbundantObservedMZ = 0;
	int iMostAbundantPeakIndex = 0;
	for (n = 0; n < (int)vvdYionMass.size(); ++n)
	{
		for (z = 1; z <= iParentChargeState; ++z)
		{
			ProductIon currentIon;
			currentIon.setProductIon('y', n+1, z);
			if(findProductIon(vvdYionMass[n], vvdYionProb[n], z,
					dScoreWeight, dMZError, dMostAbundantObservedMZ, iMostAbundantPeakIndex))
			{
				currentIon.setObservedInfo(dMZError, dScoreWeight, dMostAbundantObservedMZ, iMostAbundantPeakIndex);
				vFoundIons.push_back(currentIon);
			}
		}
	}
	for (n = 0; n < (int)vvdBionMass.size(); ++n)
	{
		for (z = 1; z <= iParentChargeState; ++z)
		{
			ProductIon currentIon;
			currentIon.setProductIon('b', n+1, z);
			if(findProductIon(vvdBionMass[n], vvdBionProb[n], z,
					dScoreWeight, dMZError, dMostAbundantObservedMZ, iMostAbundantPeakIndex))
			{
				currentIon.setObservedInfo(dMZError, dScoreWeight, dMostAbundantObservedMZ, iMostAbundantPeakIndex);
				vFoundIons.push_back(currentIon);
			}
		}
	}
	int i;
	int j;
	for (i = 0; i < (int)vFoundIons.size(); ++i)
	{
		vFoundIons[i].setComplementaryFragmentObserved(false);
	}
	for (i = 0; i < (int)vFoundIons.size(); ++i)
	{
		for (j = i+1; j < (int)vFoundIons.size(); ++j)
		{
			if(vFoundIons[i].getIonNumber() + vFoundIons[j].getIonNumber() == iPeptideLength)
			{
				if((vFoundIons[i].getIonType() == 'y' && vFoundIons[j].getIonType() == 'b')
				 ||(vFoundIons[i].getIonType() == 'b' && vFoundIons[j].getIonType() == 'y'))
				{
					vFoundIons[i].setComplementaryFragmentObserved(true);
					vFoundIons[j].setComplementaryFragmentObserved(true);
				}
			}

		}
	}

//	print out found product ions

	/*
//	LessProductIon lessIon("mass");
//	sort(vFoundIons.begin(), vFoundIons.end(), lessIon);
	cout << "Ion	Charge	MZError	Weight	MZ	Int" << endl;
	for (i = 0; i < (int)vFoundIons.size(); ++i)
	{

		cout << vFoundIons[i].getIonType() << vFoundIons[i].getIonNumber() << "\t"
		<< vFoundIons[i].getCharge() << "\t" <<
//		(vFoundIons[i].getMZError()/vFoundIons[i].getMostAbundantMZ())*1000000 <<
		vFoundIons[i].getMZError() <<
		"\t" << vFoundIons[i].getScoreWeight() << "\t"
		<< vFoundIons[i].getMostAbundantMZ() << "\t"
		<< vFoundIons[i].getMostAbundantInt() << endl;
	}
	*/


	double dAverageMZError = 0;
	for (i = 0; i < (int)vFoundIons.size(); ++i)
	{
		dAverageMZError += vFoundIons[i].getMZError();
	}
	dAverageMZError = dAverageMZError / (double)vFoundIons.size();

	double dPrimaryScore = 0;
	double dIonScore = 0;
	double dBonus4ComplementaryFragmentObserved = 1.0;
	for (i = 0; i < (int)vFoundIons.size(); ++i)
	{
		if(vFoundIons[i].getComplementaryFragmentObserved())
		{
			dBonus4ComplementaryFragmentObserved = 2.0;
		}
		else
		{
			dBonus4ComplementaryFragmentObserved = 1.0;
		}


		dIonScore = ProNovoConfig::scoreError(fabs(vFoundIons[i].getMZError()-dAverageMZError)) + vdNormalizedIntensity[ vFoundIons[i].getMostAbundantPeakIndex() ];
		dPrimaryScore += dIonScore * vFoundIons[i].getScoreWeight() * dBonus4ComplementaryFragmentObserved;
	}

	return dPrimaryScore;
}

double PeptideScorer::primaryScorePetide4( const string & sPeptide)
{
	int iPeptideLength = 0;
	for (int i = 0; i < (int)sPeptide.length(); ++i)
	{
		if(isalpha(sPeptide[i]))
		{
			iPeptideLength = iPeptideLength + 1;
		}
	}
	vector< vector<double> > vvdYionMass;
	vector< vector<double> > vvdYionProb;
	vector< vector<double> > vvdBionMass;
	vector< vector<double> > vvdBionProb;

	vector<double> vdYmass;
	vector<double> vdBmass;
	// compute b and y ion masses
	double dPeptideMass;
	computeYBionMass(sPeptide, vdYmass, vdBmass, dPeptideMass);
	vector<double> vdTempMass;
	vector<double> vdTempProb;
	for (unsigned int k = 0; k < vdYmass.size(); ++k)
	{
		ProNovoConfig::configIsotopologue.getPolyAveragine( k+1, vdYmass[k],
				vdTempMass, vdTempProb);
		vvdYionMass.push_back(vdTempMass);
		vvdYionProb.push_back(vdTempProb);
	}
	for (unsigned int k = 0; k < vdBmass.size(); ++k)
	{
		ProNovoConfig::configIsotopologue.getPolyAveragine( k+1, vdBmass[k],
				vdTempMass, vdTempProb);
		vvdBionMass.push_back(vdTempMass);
		vvdBionProb.push_back(vdTempProb);
	}
//	ProNovoConfig::configIsotopologue.computeProductIon( sPeptide,
//			vvdYionMass, vvdYionProb,
//			vvdBionMass, vvdBionProb);
/*
	for (int m = 0; m < (int)vvdBionMass.size(); ++m)
	{
		cout << "b " << m << " = " << vvdBionMass[m][0] << endl;
	}
*/

	int n; // Ion number starting from one
	int z; // charge state
	vector<ProductIon> vFoundIons;
	double dScoreWeight = 0;
	double dMZError = 1;
	double dMostAbundantObservedMZ = 0;
	int iMostAbundantPeakIndex = 0;
	for (n = 0; n < (int)vvdYionMass.size(); ++n)
	{
		for (z = 1; z <= iParentChargeState; ++z)
		{
			ProductIon currentIon;
			currentIon.setProductIon('y', n+1, z);
			if(findProductIon2(vvdYionMass[n], vvdYionProb[n], z,
					dScoreWeight, dMZError, dMostAbundantObservedMZ, iMostAbundantPeakIndex))
			{
				currentIon.setObservedInfo(dMZError, dScoreWeight, dMostAbundantObservedMZ, iMostAbundantPeakIndex);
				vFoundIons.push_back(currentIon);
			}
		}
	}
	for (n = 0; n < (int)vvdBionMass.size(); ++n)
	{
		for (z = 1; z <= iParentChargeState; ++z)
		{
			ProductIon currentIon;
			currentIon.setProductIon('b', n+1, z);
			if(findProductIon2(vvdBionMass[n], vvdBionProb[n], z,
					dScoreWeight, dMZError, dMostAbundantObservedMZ, iMostAbundantPeakIndex))
			{
				currentIon.setObservedInfo(dMZError, dScoreWeight, dMostAbundantObservedMZ, iMostAbundantPeakIndex);
				vFoundIons.push_back(currentIon);
			}
		}
	}
	int i;
	int j;
	for (i = 0; i < (int)vFoundIons.size(); ++i)
	{
		vFoundIons[i].setComplementaryFragmentObserved(false);
	}
	for (i = 0; i < (int)vFoundIons.size(); ++i)
	{
		for (j = i+1; j < (int)vFoundIons.size(); ++j)
		{
			if(vFoundIons[i].getIonNumber() + vFoundIons[j].getIonNumber() == iPeptideLength)
			{
				if((vFoundIons[i].getIonType() == 'y' && vFoundIons[j].getIonType() == 'b')
				 ||(vFoundIons[i].getIonType() == 'b' && vFoundIons[j].getIonType() == 'y'))
				{
					vFoundIons[i].setComplementaryFragmentObserved(true);
					vFoundIons[j].setComplementaryFragmentObserved(true);
				}
			}

		}
	}

//	print out found product ions

	/*
//	LessProductIon lessIon("mass");
//	sort(vFoundIons.begin(), vFoundIons.end(), lessIon);
	cout << "Ion	Charge	MZError	Weight	MZ	Int" << endl;
	for (i = 0; i < (int)vFoundIons.size(); ++i)
	{

		cout << vFoundIons[i].getIonType() << vFoundIons[i].getIonNumber() << "\t"
		<< vFoundIons[i].getCharge() << "\t" <<
//		(vFoundIons[i].getMZError()/vFoundIons[i].getMostAbundantMZ())*1000000 <<
		vFoundIons[i].getMZError() <<
		"\t" << vFoundIons[i].getScoreWeight() << "\t"
		<< vFoundIons[i].getMostAbundantMZ() << "\t"
		<< vFoundIons[i].getMostAbundantInt() << endl;
	}
	*/


	double dAverageMZError = 0;
	for (i = 0; i < (int)vFoundIons.size(); ++i)
	{
		dAverageMZError += vFoundIons[i].getMZError();
	}
	dAverageMZError = dAverageMZError / (double)vFoundIons.size();

	double dPrimaryScore = 0;
	double dIonScore = 0;
	double dBonus4ComplementaryFragmentObserved = 1.0;
	for (i = 0; i < (int)vFoundIons.size(); ++i)
	{
		if(vFoundIons[i].getComplementaryFragmentObserved())
		{
			dBonus4ComplementaryFragmentObserved = 2.0;
		}
		else
		{
			dBonus4ComplementaryFragmentObserved = 1.0;
		}


		dIonScore = ProNovoConfig::scoreError(fabs(vFoundIons[i].getMZError()-dAverageMZError)) + vdNormalizedIntensity[ vFoundIons[i].getMostAbundantPeakIndex() ];
		dPrimaryScore += dIonScore * vFoundIons[i].getScoreWeight() * dBonus4ComplementaryFragmentObserved;
	}

	return dPrimaryScore;
}
int PeptideScorer::getMaxValueIndex(const vector<double> & vdData)
{
	int iMaxIndex = 0;
	double dMaxValue = 0;
	for( unsigned int i = 0; i < vdData.size(); ++i )
	{
		if( dMaxValue < vdData[i] )
		{
			dMaxValue = vdData[i];
			iMaxIndex = (int)i;
		}
	}
	return iMaxIndex;

}

int PeptideScorer::getMaxValueIndex(const vector<double> & vdData, const vector<int> & viIndexRange)
{
	int iIndex = 0;
	int iMaxIndex = 0;
	double dMaxValue = 0;
	if( viIndexRange.size() == 0)
	{
		cout << "Error: viIndexRange is empty" << endl;
	}
	for( unsigned int i = 0; i < viIndexRange.size(); ++i )
	{
		iIndex = viIndexRange[i];
		if(0 <= iIndex && iIndex < (int)vdData.size())
		{
			if( dMaxValue < vdData[iIndex] )
			{
				dMaxValue = vdData[iIndex];
				iMaxIndex = iIndex;
			}
		}
		else
			cout << "Error: viIndexRange is not valid" << endl;
	}
	return iMaxIndex;

}

bool PeptideScorer::findProductIon(const vector<double> & vdIonMass,
						 const vector<double> & vdIonProb,
						 const int & iCharge,
						 double & dScoreWeight,
						 double & dAverageMZError,
						 double & dMostAbundantObservedMZ,
						 int & iMostAbundantPeakIndex)
{
	int iIndex4MaxInt = getMaxValueIndex(vdIonProb);
	double dMaxIntExpectedMZ = (vdIonMass[iIndex4MaxInt]/(double)iCharge) + dProtonMass;
	int iIndex4SelectedFound = 0;
	dScoreWeight = 0;
	dAverageMZError = dMassTolerance;

	// search for the most abundant peak
	if(!searchMZ(dMaxIntExpectedMZ, iIndex4SelectedFound))
	{
		return false;
	}
	iMostAbundantPeakIndex = iIndex4SelectedFound;
	dMostAbundantObservedMZ = vdMZ[iIndex4SelectedFound];
	double dMostAbundantObservedIntensity = vdIntensity[iIndex4SelectedFound];
	double dMostAbundantMZError = dMostAbundantObservedMZ - dMaxIntExpectedMZ;


	// compute expected MZ and intensity for this product ion
	vector<bool> vbObserved(vdIonProb.size(), false);
	vector<double> vdObservedMZ(vdIonProb.size(), 0);
	vector<double> vdObservedRelativeInt(vdIonProb.size(), 0);
	vector<double> vdMZError(vdIonProb.size(), dMassTolerance);
	vector<double> vdExpectedMZ(vdIonProb.size(), 0);
	vector<double> vdExpectedRelativeInt(vdIonProb.size(), 0);
	// a expected ion have to exceed dMinRelativeExpectedInt to be considered
	int i;
	for (i = 0; i < (int)vdIonProb.size(); ++i)
	{
		vdExpectedMZ[i]= vdIonMass[i]/(double)iCharge + dProtonMass;
		vdExpectedRelativeInt[i] = vdIonProb[i]/vdIonProb[iIndex4MaxInt];
	}

	// set max int value
	vbObserved[iIndex4MaxInt] = true;
	vdObservedMZ[iIndex4MaxInt] = dMostAbundantObservedMZ;
	vdMZError[iIndex4MaxInt] = dMostAbundantMZError;
	vdObservedRelativeInt[iIndex4MaxInt] = 1.0;

	if(vbObserved.size() == 1)
	{
		// there is only one expected peak in the isotopic distribution
		dScoreWeight = 1.0;
		dAverageMZError = dMostAbundantMZError;
		return true;
	}

	// search for  ions on the left of the most abundant peak
	vector<int> viIndex4Found;
	double dShiftedExpectedMZ;
	unsigned int j;
	double dMinRelativeIntRatio = 0.5;
	double dCurrentRelativeInt = 0;
	for (i = iIndex4MaxInt+1; i < (int)vdExpectedMZ.size(); ++i)
	{
		// shift expected MZ by the error of the most abundant ion and search for it
		dShiftedExpectedMZ = vdExpectedMZ[i] + dMostAbundantMZError;
		if(binarySearch(dShiftedExpectedMZ, vdMZ, dMassTolerance/2, viIndex4Found))
		{
			// the observed peaks may be noise peak
			// filter them by requiring their observed relative intensity larger than
			// a minimum ratio of their expected relative intensity
			// and less than a 150% or five times of their expected relative int, whichever smaller
			for (j = 0; j < viIndex4Found.size(); ++j)
			{
				iIndex4SelectedFound = viIndex4Found[j];
				dCurrentRelativeInt = vdIntensity[iIndex4SelectedFound]/dMostAbundantObservedIntensity;
				if(dCurrentRelativeInt > vdExpectedRelativeInt[i]*dMinRelativeIntRatio
				&& dCurrentRelativeInt < min(vdExpectedRelativeInt[i]*5, 1.5))
				{
					vbObserved[i] = true;
					vdObservedMZ[i] = vdMZ[iIndex4SelectedFound];
					vdMZError[i] = vdObservedMZ[i] - vdExpectedMZ[i];
					vdObservedRelativeInt[i] = dCurrentRelativeInt;
					break;
				}
			}
		}
		else
			break;
	}
	// search for ions on the left of the most abundant peak
	// identical to the above function except the index
	for (i = iIndex4MaxInt-1; i >= 0; --i)
	{
		// shift expected MZ by the error of the most abundant ion and search for it
		dShiftedExpectedMZ = vdExpectedMZ[i] + dMostAbundantMZError;
		if(binarySearch(dShiftedExpectedMZ, vdMZ, dMassTolerance/2, viIndex4Found))
		{
			// the observed peaks may be noise peak
			// filter them by requiring their observed relative intensity larger than
			// a minimum ratio of their expected relative intensity
			// and less than a 150% or five times of their expected relative int, whichever smaller
			for (j = 0; j < viIndex4Found.size(); ++j)
			{
				iIndex4SelectedFound = viIndex4Found[j];
				dCurrentRelativeInt = vdIntensity[iIndex4SelectedFound]/dMostAbundantObservedIntensity;
				if(dCurrentRelativeInt > vdExpectedRelativeInt[i]*dMinRelativeIntRatio
				&& dCurrentRelativeInt < min(vdExpectedRelativeInt[i]*5, 1.5))
				{
					vbObserved[i] = true;
					vdObservedMZ[i] = vdMZ[iIndex4SelectedFound];
					vdMZError[i] = vdObservedMZ[i] - vdExpectedMZ[i];
					vdObservedRelativeInt[i] = dCurrentRelativeInt;
					break;
				}
			}
		}
		else
			break;
	}


	// calculate score weight for this product ion
	// this formula is still very ad hoc empirical
	dScoreWeight = 0;
	vector<double> vdTempExpectedRelativeInt = vdExpectedRelativeInt;
	vdTempExpectedRelativeInt[iIndex4MaxInt] = 0;
	int iIndex4SecondHighestInt = getMaxValueIndex(vdTempExpectedRelativeInt);
	double dDetectionLimit4RelativeIntensity = 0.5;
	if(vbObserved[iIndex4SecondHighestInt])
	{
		// if the second highest peak is found
		// to be implemented
		dScoreWeight = 2.0;
	}
	else
	{
		if(vdExpectedRelativeInt[iIndex4SecondHighestInt]>dDetectionLimit4RelativeIntensity)
		{
			// if the second highest peak is not found and is expected to
			// be found because its relative intensity is higher than the detection limit
			dScoreWeight = 0.5;
		}
		else
		{
			// if the second highest peak is not found and is expected to not be found
			dScoreWeight = 1.0;
		}
	}

	// test whether the iCharge is consistant with viZinput
	// if not, lower the dScoreWeight
	if(viZinput[iMostAbundantPeakIndex] != 0 )
	{
		if(viZinput[iMostAbundantPeakIndex] != iCharge)
		{
			dScoreWeight = dScoreWeight / 2;
		}
	}



	// calculate average mass error for this product ion
	double dTotalRelativeIntensity = 0;
	double dTotalMZError = 0;
	for (i = 0; i < (int)vdMZError.size(); ++i)
	{
		if(vbObserved[i])
		{
			dTotalMZError += vdMZError[i] * vdExpectedRelativeInt[i];
			dTotalRelativeIntensity += vdExpectedRelativeInt[i];
		}
	}
	dAverageMZError = dTotalMZError / dTotalRelativeIntensity;
//	dAverageMZError = dMostAbundantMZError;


	return true;
}

bool PeptideScorer::findProductIon2(const vector<double> & vdIonMass,
						 const vector<double> & vdIonProb,
						 const int & iCharge,
						 double & dScoreWeight,
						 double & dAverageMZError,
						 double & dMostAbundantObservedMZ,
						 int & iMostAbundantPeakIndex)
{
	int iIndex4MaxInt = getMaxValueIndex(vdIonProb);
	double dMaxIntExpectedMZ = (vdIonMass[iIndex4MaxInt]/(double)iCharge) + dProtonMass;
	int iIndex4SelectedFound = 0;
	dScoreWeight = 0;
	dAverageMZError = dMassTolerance;

	// search for the most abundant peak
	if(!searchMZ(dMaxIntExpectedMZ, iIndex4SelectedFound))
	{
		return false;
	}
	iMostAbundantPeakIndex = iIndex4SelectedFound;
	dMostAbundantObservedMZ = vdMZ[iIndex4SelectedFound];
	double dMostAbundantObservedIntensity = vdIntensity[iIndex4SelectedFound];
	double dMostAbundantMZError = dMostAbundantObservedMZ - dMaxIntExpectedMZ;


	// compute expected MZ and intensity for this product ion
	vector<bool> vbObserved(vdIonProb.size(), false);
	vector<double> vdObservedMZ(vdIonProb.size(), 0);
	vector<double> vdObservedRelativeInt(vdIonProb.size(), 0);
	vector<double> vdMZError(vdIonProb.size(), dMassTolerance);
	vector<double> vdExpectedMZ(vdIonProb.size(), 0);
	vector<double> vdExpectedRelativeInt(vdIonProb.size(), 0);
	// a expected ion have to exceed dMinRelativeExpectedInt to be considered
	int i;
	for (i = 0; i < (int)vdIonProb.size(); ++i)
	{
		vdExpectedMZ[i]= vdIonMass[i]/(double)iCharge + dProtonMass;
		vdExpectedRelativeInt[i] = vdIonProb[i]/vdIonProb[iIndex4MaxInt];
	}

	// set max int value
	vbObserved[iIndex4MaxInt] = true;
	vdObservedMZ[iIndex4MaxInt] = dMostAbundantObservedMZ;
	vdMZError[iIndex4MaxInt] = dMostAbundantMZError;
	vdObservedRelativeInt[iIndex4MaxInt] = 1.0;

	if(vbObserved.size() == 1)
	{
		// there is only one expected peak in the isotopic distribution
		dScoreWeight = 1.0;
		dAverageMZError = dMostAbundantMZError;
		return true;
	}

	// search for  ions on the left of the most abundant peak
	vector<int> viIndex4Found;
	double dShiftedExpectedMZ;
	unsigned int j;
	double dMinRelativeIntRatio = 0.5;
	double dCurrentRelativeInt = 0;
	for (i = iIndex4MaxInt+1; i < (int)vdExpectedMZ.size(); ++i)
	{
		// shift expected MZ by the error of the most abundant ion and search for it
		dShiftedExpectedMZ = vdExpectedMZ[i] + dMostAbundantMZError;
		if(binarySearch(dShiftedExpectedMZ, vdMZ, dMassTolerance/2, viIndex4Found))
		{
			// the observed peaks may be noise peak
			// filter them by requiring their observed relative intensity larger than
			// a minimum ratio of their expected relative intensity
			// and less than a 150% or five times of their expected relative int, whichever smaller
			for (j = 0; j < viIndex4Found.size(); ++j)
			{
				iIndex4SelectedFound = viIndex4Found[j];
				dCurrentRelativeInt = vdIntensity[iIndex4SelectedFound]/dMostAbundantObservedIntensity;
				if(dCurrentRelativeInt > vdExpectedRelativeInt[i]*dMinRelativeIntRatio
				&& dCurrentRelativeInt < min(vdExpectedRelativeInt[i]*5, 1.5))
				{
					vbObserved[i] = true;
					vdObservedMZ[i] = vdMZ[iIndex4SelectedFound];
					vdMZError[i] = vdObservedMZ[i] - vdExpectedMZ[i];
					vdObservedRelativeInt[i] = dCurrentRelativeInt;
					break;
				}
			}
		}
		else
			break;
	}
	// search for ions on the left of the most abundant peak
	// identical to the above function except the index
	for (i = iIndex4MaxInt-1; i >= 0; --i)
	{
		// shift expected MZ by the error of the most abundant ion and search for it
		dShiftedExpectedMZ = vdExpectedMZ[i] + dMostAbundantMZError;
		if(binarySearch(dShiftedExpectedMZ, vdMZ, dMassTolerance/2, viIndex4Found))
		{
			// the observed peaks may be noise peak
			// filter them by requiring their observed relative intensity larger than
			// a minimum ratio of their expected relative intensity
			// and less than a 150% or five times of their expected relative int, whichever smaller
			for (j = 0; j < viIndex4Found.size(); ++j)
			{
				iIndex4SelectedFound = viIndex4Found[j];
				dCurrentRelativeInt = vdIntensity[iIndex4SelectedFound]/dMostAbundantObservedIntensity;
				if(dCurrentRelativeInt > vdExpectedRelativeInt[i]*dMinRelativeIntRatio
				&& dCurrentRelativeInt < min(vdExpectedRelativeInt[i]*5, 1.5))
				{
					vbObserved[i] = true;
					vdObservedMZ[i] = vdMZ[iIndex4SelectedFound];
					vdMZError[i] = vdObservedMZ[i] - vdExpectedMZ[i];
					vdObservedRelativeInt[i] = dCurrentRelativeInt;
					break;
				}
			}
		}
		else
			break;
	}


	// calculate score weight for this product ion
	// this formula is still very ad hoc empirical
	dScoreWeight = 0;
	vector<double> vdTempExpectedRelativeInt = vdExpectedRelativeInt;
	vdTempExpectedRelativeInt[iIndex4MaxInt] = 0;
	int iIndex4SecondHighestInt = getMaxValueIndex(vdTempExpectedRelativeInt);
	double dDetectionLimit4RelativeIntensity = 0.5;
	if(vbObserved[iIndex4SecondHighestInt])
	{
		// if the second highest peak is found
		// to be implemented
		dScoreWeight = 2.0;
	}
	else
	{
		if(vdExpectedRelativeInt[iIndex4SecondHighestInt]>dDetectionLimit4RelativeIntensity)
		{
			// if the second highest peak is not found and is expected to
			// be found because its relative intensity is higher than the detection limit
			dScoreWeight = 0.5;
		}
		else
		{
			// if the second highest peak is not found and is expected to not be found
			dScoreWeight = 1.0;
		}
	}

	// test whether the iCharge is consistant with viZinput
	// if not, lower the dScoreWeight
	if(viZinput[iMostAbundantPeakIndex] != 0 )
	{
		if(viZinput[iMostAbundantPeakIndex] != iCharge)
		{
			dScoreWeight = dScoreWeight / 2;
		}
	}

	// search for H2O loss and NH3 loss, only one will be considered, because they are off by only 1 Da
	int iIndex4MostAbundunt;
	if(searchMZ(dMostAbundantObservedMZ-ProNovoConfig::getMassH2O()/(double)iCharge, iIndex4MostAbundunt))
	{
		dScoreWeight += 0.2;
	}
	else if(searchMZ(dMostAbundantObservedMZ-ProNovoConfig::getMassNH3()/(double)iCharge, iIndex4MostAbundunt))
	{
		dScoreWeight += 0.2;
	}

	// calculate average mass error for this product ion
	double dTotalRelativeIntensity = 0;
	double dTotalMZError = 0;
	for (i = 0; i < (int)vdMZError.size(); ++i)
	{
		if(vbObserved[i])
		{
			dTotalMZError += vdMZError[i] * vdExpectedRelativeInt[i];
			dTotalRelativeIntensity += vdExpectedRelativeInt[i];
		}
	}
	dAverageMZError = dTotalMZError / dTotalRelativeIntensity;
//	dAverageMZError = dMostAbundantMZError;


	return true;
}
