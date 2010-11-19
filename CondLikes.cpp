#include "CondLikes.h"

CondLikes::CondLikes(Alignment *ap) {

	// initialize some important variables
	alignmentPtr       = ap;
	numTaxa            = alignmentPtr->getNumTaxa();
	numSites           = alignmentPtr->getNumCodonSites();
	numNodes           = 2 * numTaxa - 2;
	numStates          = alignmentPtr->getNumSenseCodons();
	
	// allocate memory
	int oneNodeSize = numSites * numStates;
	int oneClSize   = numNodes * oneNodeSize;
	cls = new double[2 * oneClSize];
	for (int i = 0; i < 2 * oneClSize; i++)
		cls[i] = 0.0;
	clsUpL = new double[2 * oneClSize];
	for (int i = 0; i < 2 * oneClSize; i++)
		clsUpL[i] = 0.0;
	clsUpR = new double[2 * oneClSize];
	for (int i = 0; i < 2 * oneClSize; i++)
		clsUpR[i] = 0.0;
		
	// initialize pointer information
	for (int n = 0; n < 2; n++)
		{
		clsPtr[n] = new double**[numNodes];
		clsPtr[n][0] = new double*[numNodes * numSites];
		for (int i=1; i < numNodes; i++)
			clsPtr[n][i] = clsPtr[n][i - 1] + numSites;
		for (int i=0; i < numNodes; i++)
			for (int j = 0; j < numSites; j++)
				clsPtr[n][i][j] = &cls[n * oneClSize + i * oneNodeSize + j * numStates];

		clsPtrUpL[n] = new double**[numNodes];
		clsPtrUpL[n][0] = new double*[numNodes * numSites];
		for (int i=1; i < numNodes; i++)
			clsPtrUpL[n][i] = clsPtrUpL[n][i-1] + numSites;
		for (int i = 0; i < numNodes; i++)
			for (int j = 0; j < numSites; j++)
				clsPtrUpL[n][i][j] = &clsUpL[n * oneClSize + i * oneNodeSize + j * numStates];

		clsPtrUpR[n] = new double**[numNodes];
		clsPtrUpR[n][0] = new double*[numNodes * numSites];
		for (int i = 1; i < numNodes; i++)
			clsPtrUpR[n][i] = clsPtrUpR[n][i - 1] + numSites;
		for (int i = 0; i < numNodes; i++)
			for (int j = 0; j < numSites; j++)
				clsPtrUpR[n][i][j] = &clsUpR[n * oneClSize + i * oneNodeSize + j * numStates];
		}
	
	// initialize tip conditional likelihoods
	for (int i = 0; i < numTaxa; i++)
	{
		for (int j = 0; j < numSites; j++)
		{
			double* p0 = clsPtr[0][i][j];
			double* p1 = clsPtr[1][i][j];
			double* uL0 = clsPtrUpL[0][i][j];
			double* uL1 = clsPtrUpL[1][i][j];
			double* uR0 = clsPtrUpR[0][i][j];
			double* uR1 = clsPtrUpR[1][i][j];
			
			MbBitfield *bf = alignmentPtr->getCodonMatrixEntry(i, j);
			for (int s = 0; s < alignmentPtr->getNumSenseCodons(); s++)
			{
				if (bf->isBitSet(s) == true)
				{
					p0[s] = 1.0;
					p1[s] = 1.0;
					uL0[s] = 1.0;
					uL1[s] = 1.0;
					uR0[s] = 1.0;
					uR1[s] = 1.0;
				}
			}
		}
	}
		
	// allocate memory for the scalers
	clScalers = new double[2 * numNodes * numSites];
	for (int i = 0; i < 2 * numNodes * numSites; i++)
		clScalers[i] = 0.0;
	for (int n = 0; n < 2; n++)
		{
		clScalerPtr[n] = new double**[numNodes];
		clScalerPtr[n][0] = new double*[numNodes * numSites];
		for (int i=1; i<numNodes; i++)
			clScalerPtr[n][i] = clScalerPtr[n][i - 1] + numSites;
		for (int i = 0; i < numNodes; i++)
			for (int j = 0; j < numSites; j++)
				clScalerPtr[n][i][j] = &clScalers[n * (numNodes * numSites) + i * (numSites) + j];
		}

	clScalersUpL = new double[2 * numNodes * numSites];
	for (int i = 0; i < 2 * numNodes * numSites; i++)
		clScalersUpL[i] = 0.0;
	for (int n = 0; n < 2; n++)
		{
		clScalerPtrUpL[n] = new double**[numNodes];
		clScalerPtrUpL[n][0] = new double*[numNodes * numSites];
		for (int i = 1; i < numNodes; i++)
			clScalerPtrUpL[n][i] = clScalerPtrUpL[n][i - 1] + numSites;
		for (int i = 0; i < numNodes; i++)
			for (int j = 0; j < numSites; j++)
				clScalerPtrUpL[n][i][j] = &clScalersUpL[n * (numNodes * numSites) + i * (numSites) + j];
		}
		
	clScalersUpR = new double[2 * numNodes * numSites];
	for (int i = 0; i < 2 * numNodes * numSites; i++)
		clScalersUpR[i] = 0.0;
	for (int n = 0; n < 2; n++)
		{
		clScalerPtrUpR[n] = new double**[numNodes];
		clScalerPtrUpR[n][0] = new double*[numNodes * numSites];
		for (int i = 1; i < numNodes; i++)
			clScalerPtrUpR[n][i] = clScalerPtrUpR[n][i - 1] + numSites;
		for (int i = 0; i < numNodes; i++)
			for (int j = 0; j < numSites; j++)
				clScalerPtrUpR[n][i][j] = &clScalersUpR[n * (numNodes * numSites) + i * (numSites) + j];
		}
	
	// print();
}

CondLikes::~CondLikes(void) {

	delete [] cls;
	for (int n=0; n<2; n++)
		{
		delete [] clsPtr[n][0];
		delete [] clsPtr[n];
		}
	delete [] clScalers;
	for (int n=0; n<2; n++)
		{
		delete [] clScalerPtr[n][0];
		delete [] clScalerPtr[n];
		}

	delete [] clsUpL;
	for (int n=0; n<2; n++)
		{
		delete [] clsPtrUpL[n][0];
		delete [] clsPtrUpL[n];
		}
	delete [] clScalersUpL;
	for (int n=0; n<2; n++)
		{
		delete [] clScalerPtrUpL[n][0];
		delete [] clScalerPtrUpL[n];
		}

	delete [] clsUpR;
	for (int n=0; n<2; n++)
		{
		delete [] clsPtrUpR[n][0];
		delete [] clsPtrUpR[n];
		}
	delete [] clScalersUpR;
	for (int n=0; n<2; n++)
		{
		delete [] clScalerPtrUpR[n][0];
		delete [] clScalerPtrUpR[n];
		}
}

void CondLikes::print(void) {

	for (int j = 0; j < numSites; j++)
		{
		std::cout << std::setw(5) << j << " -- ";
		for (int i = 0; i < numTaxa; i++)
			{
			double* p0 = getClPtr(0, i, j);
			for (int s = 0; s < numStates; s++)
				std::cout << std::fixed << std::setprecision(0) << p0[s];
			std::cout << " ";
			}
		std::cout << std::endl;
		}
}


