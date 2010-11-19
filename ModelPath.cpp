/*
 * ModelPath.cpp
 *
 *  Created on: Jan 24, 2010
 *      Author: mikee
 */

#include "ModelPath.h"

#define ZERO_ENTRY		0
#define SYN_TI_ENTRY	1
#define SYN_TV_ENTRY	2
#define NONSYN_TI_ENTRY	3
#define NONSYN_TV_ENTRY	4

#define RESCALE_LIKES

ModelPath::ModelPath(Alignment *ap, MbRandom *rp, Settings* sp, Tree* tp) {

	std::cout << "INITIALIZING MODEL: PATH SAMPLER\n";

	// keep track of pointers to important objects
	alignmentPtr	= ap;
	ranPtr			= rp;
	treePtr			= tp;
	settingsPtr		= sp;
	numCodonSites	= alignmentPtr->getNumCodonSites();

	// allocate a vector holding the stationary probabilities
	stationaryProbs = new double[numCodonSites];

	if (settingsPtr->getUseECM() == false)
	{
		// add the parameters to the phylogenetic model
		parameters.push_back(treePtr);
		parameters.push_back(new CodonFrequencies(ranPtr, "Codon Frequencies", treePtr, alignmentPtr));
		parameters.push_back(new Kappa(ranPtr, treePtr, "Kappa"));
		parameters.push_back(new Omega(ranPtr, treePtr, "Omega"));
		std::cout << "\n";

		// initialize the proposal probabilities
		proposalProbs.push_back( 0.0 ); // probability of proposing a change to the tree
		proposalProbs.push_back( 1.0 ); // probability of proposing a change to the codon frequencies
		proposalProbs.push_back( 1.0 ); // probability of proposing a change to the transition/transversion rate ratio
		proposalProbs.push_back( 1.0 ); // probability of proposing a change to the dN/dS rate ratio

		double sum = 0.0;
		for (int i = 0; i < (int)proposalProbs.size(); i++)
			sum += proposalProbs[i];
		for (int i = 0; i < (int)proposalProbs.size(); i++)
			proposalProbs[i] /= sum;

		// set up the qTemplate rate matrix (which runs fillInRateMatrix() in the process
		initializeRateMatrix(alignmentPtr->getNumSenseCodons());
	}

	else if (settingsPtr->getUseECM() == true)
	{
		// set up pi_j
		std::cout << "INITIALIZING EMPIRICAL CODON MATRIX: CODON FREQUENCIES\n";
		std::vector<double> ECMfreqs = initializeEmpiricalCodonFreqs(settingsPtr->getCodonFreqFileName(), alignmentPtr->getNumSenseCodons());

		// add the parameters to the phylogenetic model
		parameters.push_back(treePtr);
		parameters.push_back(new CodonFrequencies(ECMfreqs, ranPtr, "Codon Frequencies", treePtr, alignmentPtr));
		std::cout << "\n";

		// initialize the proposal probabilities
		proposalProbs.push_back( 0.0 ); // probability of proposing a change to the tree
		proposalProbs.push_back( 1.0 ); // probability of proposing a change to the codon frequencies

		double sum = 0.0;
		for (int i = 0; i < (int)proposalProbs.size(); i++)
			sum += proposalProbs[i];
		for (int i = 0; i < (int)proposalProbs.size(); i++)
			proposalProbs[i] /= sum;

		// set up q_ij = s_ij * pi_j
		std::cout << "INITIALIZING EMPIRICAL CODON MATRIX: EXCHANGIBILITY VALUES\n";
		//initializeEmpiricalCodonMatrix(settingsPtr->getCodonMatrixFileName(), alignmentPtr->getNumSenseCodons());

	}

	// instantiate the conditional likelihoods
	condLikes = new CondLikes(alignmentPtr);

	// instantiate the transition probabilities
	tiMatrix = new MbTransitionMatrix(q, false);
	tiProbs = new TiProbs( 2 * alignmentPtr->getNumTaxa() - 2, tiMatrix, treePtr);

#	if defined(RESCALE_LIKES)
	// allocate a vector that holds (temporarily) the log scalers for log likelihood calculations
	lnScalers = new double[numCodonSites];
#	endif

	// print the parameters
	//for (std::vector<Parm*>::iterator p=parameters.begin(); p != parameters.end(); p++)
	//	(*p)->print();

	std::cout << "\n";
}

ModelPath::~ModelPath(void) {

	delete condLikes;
	delete tiProbs;
	for (std::vector<Parm*>::iterator p=parameters.begin(); p != parameters.end(); p++)
		delete (*p);
	delete [] stationaryProbs;
#	if defined(RESCALE_LIKES)
	delete [] lnScalers;
#	endif
}

void ModelPath::initializeRateMatrix(int nc) {

	numSenseCodons		= nc;
	numStates 			= nc;
	uniformizationRate	= 0.0;

	q = MbMatrix<double>(numSenseCodons, numSenseCodons);
	r = MbMatrix<double>(numSenseCodons, numSenseCodons);
	for (int i = 0; i < numSenseCodons; i++)
	{
		for (int j = 0; j < numSenseCodons; j++)
		{
			q[i][j] = 0.0;
		}
	}

	// fill in the template
	qTemplate.resize(numSenseCodons);
	for (int i = 0; i < numSenseCodons; i++)
	{
		qTemplate[i].resize(numSenseCodons, ZERO_ENTRY);
	}

	for (int i=0; i<numSenseCodons; i++)
	{
		int a1 = alignmentPtr->getGeneticCode()->getCodonNuc(i, 0);
		int a2 = alignmentPtr->getGeneticCode()->getCodonNuc(i, 1);
		int a3 = alignmentPtr->getGeneticCode()->getCodonNuc(i, 2);
		int aAa = alignmentPtr->getGeneticCode()->getTransCodon(i);

		for (int j=0; j<numSenseCodons; j++)
		{
			int b1 = alignmentPtr->getGeneticCode()->getCodonNuc(j, 0);
			int b2 = alignmentPtr->getGeneticCode()->getCodonNuc(j, 1);
			int b3 = alignmentPtr->getGeneticCode()->getCodonNuc(j, 2);
			int bAa = alignmentPtr->getGeneticCode()->getTransCodon(j);

			int numChanges = 0;
			if (a1 != b1)
				numChanges++;
			if (a2 != b2)
				numChanges++;
			if (a3 != b3)
				numChanges++;
			if (numChanges == 1)
			{
				bool isSynonymous = true;
				if (aAa != bAa)
					isSynonymous = false;

				bool isTransition = true;
				if ( isTransversion(a1, b1) == true )
					isTransition = false;

				if ( isSynonymous == true && isTransition == true )
					qTemplate[i][j] = SYN_TI_ENTRY;
				else if ( isSynonymous == true && isTransition == false )
					qTemplate[i][j] = SYN_TV_ENTRY;
				else if ( isSynonymous == false && isTransition == true )
					qTemplate[i][j] = NONSYN_TI_ENTRY;
				else if ( isSynonymous == false && isTransition == false )
					qTemplate[i][j] = NONSYN_TV_ENTRY;
			}
		}
	}

	// fill in the rate matrix
	fillInRateMatrix();
	printQ();
	printR();

}

void ModelPath::fillInRateMatrix(void) {

	// get the parameters of the model
	std::vector<double> f	= getActiveFreqs()->getFreqs();
	double w				= getActiveDnds()->getDnds();
	double kappa			= getActiveTiTv()->getRate();

	// set diagonal elements of rate matrix to zero
	for (int i = 0; i < numSenseCodons; i++)
		q[i][i] = 0.0;

	// fill in the non-diagonal elements of the matrix, corresponding to substitutions
	double averageRate = 0.0;
	for (int i=0; i<numSenseCodons; i++)
	{
		for (int j=0; j<numSenseCodons; j++)
		{
			if ( qTemplate[i][j] != ZERO_ENTRY )
			{
				if ( qTemplate[i][j] == SYN_TV_ENTRY )
				{
					double val = f[j];
					q[i][j] = val;
					q[i][i] -= val;
					averageRate += f[i] * val;
				}
				else if ( qTemplate[i][j] == SYN_TI_ENTRY )
				{
					double val = f[j] * kappa;
					q[i][j] = val;
					q[i][i] -= val;
					averageRate += f[i] * val;
				}
				else if ( qTemplate[i][j] == NONSYN_TV_ENTRY )
				{
					double val = f[j];
					double x = val * w;
					q[i][j] = x;
					q[i][i] -= x;
					averageRate += f[i] * x;
				}
				else if ( qTemplate[i][j] == NONSYN_TI_ENTRY )
				{
					double val = f[j] * kappa;
					double x = val * w;
					q[i][j] = x;
					q[i][i] -= x;
					averageRate += f[i] * x;
				}
			}
		}
	}

	// calculate uniformization rate for R matrix
	for (int i = 0; i < numSenseCodons; i++) {
		if (q[i][i] < uniformizationRate) {
			uniformizationRate = q[i][i];
		}
	}
	uniformizationRate = -uniformizationRate;

	// fills in stochastic matrix as follows:
	// R = (1/mu) * Q * I
	for (int i = 0; i < numSenseCodons; i++) {
		for (int j = 0; j < numSenseCodons; j++) {
			r[i][j] = ((1 / uniformizationRate) * q[i][j]);
			if (i == j) r[i][i] += 1.0;
		}
	}
}


void ModelPath::fillInStationaryProbs(void) {

	std::vector<double> f = getActiveFreqs()->getFreqs();
	double p = getActiveDnds()->getDnds();
	double *x = &stationaryProbs[0];
	for (int i=0; i<numSenseCodons; i++)
	{
		(*x) = f[i] * p;
		x++;
	}
}


void ModelPath::initializeEmpiricalCodonMatrix(std::string ecmFilePath, int nc) {


	numSenseCodons		= nc;
	numStates 			= nc;
	uniformizationRate	= 0.0;

	q = MbMatrix<double>(numSenseCodons, numSenseCodons);
	r = MbMatrix<double>(numSenseCodons, numSenseCodons);
	for (int i = 0; i < numSenseCodons; i++)
	{
		for (int j = 0; j < numSenseCodons; j++)
		{
			q[i][j] = 0.0;
		}
	}

	FileMgr codonMatrixFileMgr(ecmFilePath);
	std::ifstream seqStream;

	if (codonMatrixFileMgr.openFile(seqStream) == false) {
		std::cerr << "Cannot open file \"" + codonMatrixFileMgr.getFileName() + "\"" << std::endl;
		exit(1);
	}

	// temp variables
	std::vector<double> codonFreq = getActiveFreqs()->getFreqs();
	std::string lineStr = "";
	char* tempStr = NULL;
	char* rateStr = NULL;
	double rateVal = 0.0;
	int i = 0;
	int j = 0;

	// INITIALIZE EMPIRICAL CODON MATRIX EXCHANGABILITY PARAMETERS
	while (getline(seqStream, lineStr).good() && i < numSenseCodons)
	{
		// populate a matrix row with ECM.dat values
		tempStr = (char*)lineStr.c_str();
		rateStr = strtok(tempStr, " ");
		while (rateStr != NULL && j < numSenseCodons)
		{
			//std::cout << "j: " << j << "\n";
			rateVal = atof(rateStr);
			//std::cout << "rateVal: " << rateVal << "\n";
			if (i != j)
			{
				/*
				q[i][j] += rateVal * codonFreq[j];
				q[j][i] += rateVal * codonFreq[i];
				q[i][i] -= rateVal * codonFreq[i];
				std::cout << "q[" << i << "][" << j << "]: " << q[i][j] << ", ";
				std::cout << "q[" << i << "][" << i << "]: " << q[i][i] << "\n";
				*/
				q[i][j] = rateVal;
				q[j][i] = rateVal;
			}

			else if (i == j)
				q[i][i] = rateVal;

			//std::cout << "rateStr: " << rateStr << "\n";
			rateStr = strtok(NULL, " ");
			j++;
		}
		j = 0;
		i++;
	}

	// TURN S_IJ MATRIX INTO Q_IJ MATRIX -- Q_ij = S_ij * P_j
	for (int i = 0; i < numSenseCodons; i++)
	{
		for (int j = 0; j < numSenseCodons; j++)
		{
			//if ()
			//q[i][j] =
		}
	}


	// GET UNIFORMIZATION RATE
	for (int i = 0; i < numSenseCodons; i++)
	{
		if (q[i][i] < uniformizationRate) {
			uniformizationRate = q[i][i];
		}
	}
	uniformizationRate = -uniformizationRate;

	// fills in stochastic matrix as follows:
	// R = (1/mu) * Q * I
	for (int i = 0; i < numSenseCodons; i++) {
		for (int j = 0; j < numSenseCodons; j++) {
			r[i][j] = ((1 / uniformizationRate) * q[i][j]);
			if (i == j) r[i][i] += 1.0;
		}
	}

	printQ();
	printR();
}

std::vector<double> ModelPath::initializeEmpiricalCodonFreqs(std::string ecmFilePath, int nc) {

	numSenseCodons = nc;

	FileMgr codonFreqFileMgr(ecmFilePath);
	std::ifstream seqStream;

	if (codonFreqFileMgr.openFile(seqStream) == false) {
		std::cerr << "Cannot open file \"" + codonFreqFileMgr.getFileName() + "\"" << std::endl;
		exit(1);
	}

	// temp variables
	std::vector<double> codonFreq(numSenseCodons);
	std::string lineStr = "";
	char* tempStr = NULL;
	char* rateStr = NULL;
	double rateVal = 0.0;
	double sumVal = 0.0;
	int i = 0;

	// INITIALIZE CODON FREQUENCIES
	while (getline(seqStream, lineStr).good())
	{
		tempStr = (char*)lineStr.c_str();
		rateStr = strtok(tempStr, " ");
		while (rateStr != NULL)
		{
			rateVal = atof(rateStr);
			rateStr = strtok(NULL, " ");
			codonFreq[i] = rateVal;
			sumVal += rateVal;
			if (rateVal > 1.0)
			{
				std::cout << "ERROR: sum of empirical codon frequencies exceeds 1.0.\n";
				std::exit(1);
			}
			i++;
		}
	}

	return codonFreq;
}

void ModelPath::initializeCondLikes(int codonSite) {

	// get a pointer to the current tree topology
	Topology* t = getActiveTopology();

	// update the transition probabilities, if necessary
	tiProbs->setTransitionProbs();

	// allocate some matrices to temporarily hold the transition matrices
	MbMatrix<double> tiL;
	MbMatrix<double> tiR;
	MbMatrix<double> tiA;

	// calculate the conditional likelihoods down the tree
	for (int n = 0; n < t->getNumNodes(); n++)
	{
		Node* p = t->getDownPassNode(n);
		if (p->getLft() != NULL && p->getRht() != NULL && p->getAnc() != NULL)
		{
#			if defined(RESCALE_LIKES)
			double* clS = condLikes->getClScalerPtr( p->getActiveCl(), p->getIndex(), codonSite );
#			endif
			int lftIdx = p->getLft()->getIndex();
			int rhtIdx = p->getRht()->getIndex();
			tiL = tiProbs->getTiMatrix( p->getLft()->getActiveTi(), lftIdx );
			tiR = tiProbs->getTiMatrix( p->getRht()->getActiveTi(), rhtIdx );
			double* clL = condLikes->getClPtr( p->getLft()->getActiveCl(), lftIdx, codonSite );
			double* clR = condLikes->getClPtr( p->getRht()->getActiveCl(), rhtIdx, codonSite );
			double* clP = condLikes->getClPtr( p->getActiveCl(), p->getIndex(), codonSite );

#			if defined(RESCALE_LIKES)
			double maxCl = 0.0;
#			endif
			for (int i = 0; i < numStates; i++)
			{
				double sumL = 0.0;
				double sumR = 0.0;
				for (int j = 0; j < numStates; j++)
				{
					sumL += tiL[i][j] * clL[j];
					sumR += tiR[i][j] * clR[j];
				}
				clP[i] = sumL * sumR;
#				if defined(RESCALE_LIKES)
				if (clP[i] > maxCl)
					maxCl = clP[i];
#				endif
			}
#			if defined(RESCALE_LIKES)
			double scaler = 1.0 / maxCl;
			for (int i = 0; i < numStates; i++)
				clP[i] *= scaler;
			clS[codonSite] = log(maxCl);
#			endif
		}
	}

	// calculate the conditional likelihoods up the tree
	for (int n = t->getNumNodes() - 1; n >= 0; n--)
	{
		Node* p = t->getDownPassNode(n);
		if (p->getLft() != NULL && p->getRht() != NULL && p->getAnc() != NULL)
		{
			bool pToLeft = true;
			if (p->getAnc()->getLft() != p)
				pToLeft = false;
#			if defined(RESCALE_LIKES)
			double* clSUpL = condLikes->getClScalerPtrUpL( p->getActiveCl(), p->getIndex(), codonSite );
			double* clSUpR = condLikes->getClScalerPtrUpR( p->getActiveCl(), p->getIndex(), codonSite );
#			endif
			int lftIdx = p->getLft()->getIndex();
			int rhtIdx = p->getRht()->getIndex();
			int ancIdx = p->getAnc()->getIndex();
			tiL = tiProbs->getTiMatrix( p->getLft()->getActiveTi(), lftIdx );
			tiR = tiProbs->getTiMatrix( p->getRht()->getActiveTi(), rhtIdx );
			tiA = tiProbs->getTiMatrix( p->getActiveTi(), p->getIndex() );
			double* clL = condLikes->getClPtr( p->getLft()->getActiveCl(), lftIdx, codonSite );
			double* clR = condLikes->getClPtr( p->getRht()->getActiveCl(), rhtIdx, codonSite );
			double* clA;
			if (pToLeft == true)
				clA = condLikes->getClPtrUpR( p->getAnc()->getActiveCl(), ancIdx, codonSite );
			else
				clA = condLikes->getClPtrUpL( p->getAnc()->getActiveCl(), ancIdx, codonSite );
			double* clPUpL = condLikes->getClPtrUpL( p->getActiveCl(), p->getIndex(), codonSite );
			double* clPUpR = condLikes->getClPtrUpR( p->getActiveCl(), p->getIndex(), codonSite );

			// LEFT SIDE
#			if defined(RESCALE_LIKES)
			double maxCl = 0.0;
#			endif

			for (int i = 0; i < numStates; i++)
			{
				double sumL = 0.0;
				double sumA = 0.0;
				for (int j = 0; j < numStates; j++)
				{
					sumL += tiL[i][j] * clL[j];
					sumA += tiA[i][j] * clA[j];
					/*
					std::cout << "tiL[" << i << "][" << j << "] = " << tiL[i][j] << "; " << std::endl;
					std::cout << "clL[" << j << "] = " << clL[j] << "; " << std::endl;
					std::cout << "tiA[" << i << "][" << j << "] = " << tiA[i][j] << "; " << std::endl;
					std::cout << "clA[" << j << "] = " << clA[j] << "; " << std::endl;
					std::cout << std::endl << std::endl;
					*/
				}
				clPUpL[i] = sumL * sumA;

#				if defined(RESCALE_LIKES)
				if (clPUpL[i] > maxCl)
					maxCl = clPUpL[i];
#				endif
			}

#			if defined(RESCALE_LIKES)
			double scaler = 1.0 / maxCl;
			for (int i=0; i<numStates; i++)
				clPUpL[i] *= scaler;
			clSUpL[codonSite] = log(maxCl);
#			endif


			// RIGHT SIDE
#			if defined(RESCALE_LIKES)
			maxCl = 0.0;
#			endif

			for (int i = 0; i < numStates; i++)
			{
				double sumR = 0.0;
				double sumA = 0.0;
				for (int j = 0; j < numStates; j++)
				{
					sumR += tiR[i][j] * clR[j];
					sumA += tiA[i][j] * clA[j];
				}
				clPUpR[i] = sumR * sumA;

#				if defined(RESCALE_LIKES)
				if (clPUpR[i] > maxCl)
					maxCl = clPUpR[i];
#				endif
			}

#			if defined(RESCALE_LIKES)
			scaler = 1.0 / maxCl;
			for (int i=0; i<numStates; i++)
				clPUpR[i] *= scaler;
			clSUpR[codonSite] = log(maxCl);
#			endif

		}
	}
}


void ModelPath::createPathMappings() {

	Topology* t = getActiveTopology();

	for (int i = 0; i < numCodonSites; i++)
		initializeCondLikes(i);

	fillInStationaryProbs();
	std::cout << "\n";

	for (int n = t->getNumNodes() - 1; n >= 0; n--)
	{
		Node* p = t->getDownPassNode(n);
		std::cout << "Sampling NODE STATE for " << p->getIndex() << "\n";
		//std::cout << "n = " << n << "p->getDnPassIndex() = " << p->getDnPassIndex() << "\n";
		for (int codonSite = 0; codonSite < numCodonSites; codonSite++)
		{
			sampleNodeSiteState(p, codonSite);
		}
	}

	for (int n = t->getNumNodes() - 1; n >= 0; n--)
	{
		Node* p = t->getDownPassNode(n);
		if (p->getAnc() != NULL) {
			std::cout << "Sampling BRANCH PATH for " << p->getAnc()->getIndex() << "->" << p->getIndex() << "\n";
			for (int i = 0; i < numCodonSites; i++)
			{
				sampleNodeSitePath(p, i);
			}
		}
	}

	for (int n = t->getNumNodes() - 1; n >= 0; n--)
	{
		Node* p = t->getDownPassNode(n);
		if (p->getAnc() != NULL) {
			std::cout << "PATH MAPPING, " << p->getAnc()->getIndex() << "->" << p->getIndex() << std::endl;
			p->getActivePath()->printPathHistory();
			p->getAnc()->getActivePath()->printSeqState();
			std::cout << "t" << std::setw(24) << std::setprecision(16) << p->getV() << ": PATH TRANSFORMATION COMPLETE\n";
			std::cout << std::endl;
		}
	}

	// store initial path mappings as both accepted and proposed states
	keepPathMappings();
	std::cout << "createPathMappings() complete!" << std::endl;
}


void ModelPath::keepPathMappings() {

	Topology* t = getActiveTopology();

	for (int n = t->getNumNodes() - 1; n >= 0; n--)
	{
		Node* p = t->getDownPassNode(n);
		if (p->getAnc() != NULL)
			p->getHistory()->keep();
	}
}

void ModelPath::restorePathMappings() {

	Topology* t = getActiveTopology();

	for (int n = t->getNumNodes() - 1; n >= 0; n--)
	{
		Node* p = t->getDownPassNode(n);
		if (p->getAnc() != NULL)
			p->getHistory()->restore();
	}
}

void ModelPath::changePathMappings() {

	Topology* t = getActiveTopology();

	double u = (double)rand() / (double)RAND_MAX;
	int randSite = (int)((double)rand() / (double)RAND_MAX) * numCodonSites;
	int randNode = (int)((double)rand() / (double)RAND_MAX) * t->getNumNodes();

	// branch site resampling occurs at pr(resampleRate)
	// node site & adjacent branch site resampling occurs at pr (1.0 - resampleRate)
	double resampleRate = 0.75;

	Node* p = t->getDownPassNode(randNode);

	// resample path along branch for a site
	if (u < resampleRate) {
		sampleNodeSitePath(p, randSite);
	}

	// resample node and any adjacent branches for a site
	else {
		sampleNodeSiteState(p, randSite);
		if (p->getAnc() != NULL)
			sampleNodeSitePath(p, randSite);
		if (p->getLft() != NULL)
			sampleNodeSitePath(p->getLft(), randSite);
		if (p->getRht() != NULL)
			sampleNodeSitePath(p->getRht(), randSite);
	}
}


void ModelPath::sampleNodeSiteState(Node* p, int codonSite) {

	// STEP ONE of three : sample states at internal nodes
	//
	//        accomplished by first sampling the state at the root node according to the product of the
	//        stationary probability and the conditional likelihood of each state, followed by a
	//        recursive sampling along the tree, with each case conditioned on the state at the
	//        immediate ancestral node (Nielsen, 2002)

	// sample the root state
	double u = (double)rand() / (double)RAND_MAX;
	int ancState = -1;
	double pDenSum = 0.0;
	double pNumSum = 0.0;
	double* clP;
	MbMatrix<double> tiP;

	if (p->getAnc() != NULL)
	{
		p->getActivePath()->clearSitePathHistory(codonSite);
		std::vector<std::vector<Step> > ancPath = p->getAnc()->getActivePath()->getPathHistory();

		u = (double)rand() / (double)RAND_MAX;

		ancState = ancPath[codonSite][0].getState();
		clP = condLikes->getClPtr( p->getActiveCl(), p->getIndex(), codonSite);
		tiP = tiProbs->getTiMatrix( p->getActiveTi(), p->getIndex() );

		// denominator sum
		// sum of cl(root-1, k) * P(root-1, i->k)
		for (int k = 0; k < numStates; k++)
		{
			pDenSum += clP[k] * tiP[ancState][k];
		}

		// numerator sum
		// cl(root-1, j) * P(root-1, i->j)
		for (int j = 0; j < numStates; j++)
		{
			pNumSum += clP[j] * tiP[ancState][j];
			if (u < (pNumSum / pDenSum))
			{
				p->getActivePath()->pushStep(Step(p->getV(), j, codonSite), codonSite);
				break;
			}
			else if (u > pNumSum / pDenSum && j == numStates)
			{
				std::cout << "ERROR: no new state set!\n";
			}
		}
	}

	else if (p->getAnc() == NULL)
	{
		clP = condLikes->getClPtr( p->getActiveCl(), p->getIndex(), codonSite);
		u = (double)rand() / (double)RAND_MAX;

		// denominator sum
		for (int i = 0; i < numStates; i++)
		{
			pDenSum += clP[i] * stationaryProbs[i];
		}

		// numerator sum
		for (int i = 0; i < numStates; i++) {
			pNumSum += clP[i] * stationaryProbs[i];
			if (u < (pNumSum / pDenSum)) {
				p->getActivePath()->pushStep(Step(0.0, i, codonSite), codonSite);
				break;
			}
			else if (u > pNumSum / pDenSum && i == numStates)
			{
				std::cout << "no new state set!\n";
			}
		}
	}
}

void ModelPath::sampleNodeSitePath(Node* p, int codonSite) {


	// initialize variables
	int srcState = (p->getAnc()->getActivePath()->getSitePathHistory(codonSite))[0].getState();
	int destState = (p->getActivePath()->getSitePathHistory(codonSite))[0].getState();
	double brlen = p->getV();
	double muV = uniformizationRate * brlen;

	std::vector< MbMatrix<double> > tempRExponents;
	MbMatrix<double> tempR = MbMatrix<double>(numStates, numStates, 0.0);
	MbMatrix<double> tiP = tiProbs->getTiMatrix( p->getActiveTi(), p->getIndex() );

	std::vector<double> eventTimes;
	int eventCount = 0;
	int eventCountMax = 0;

	double eventSum = 0.0;
	double u = 0.0;
	double g = tiP[srcState][destState] * (double)rand() / (double)RAND_MAX;
	//std::cout << "g: " << g << "\n";

	double stepSum = 0.0;
	double stepSample = 0.0;

	// initalize the temporary R stochastic matrix
	for (int i = 0; i < numStates; i++) {
		for (int j = 0; j < numStates; j++) {
			if (i != j) tempR[i][j] = 0.0;
		}
		tempR[i][i] = 1.0;
	}

	// sample number of events to describe path history
	//std::cout << codonSite << ": " << "path n" << p->getAnc()->getIndex() << ":s" << srcState;
	//std::cout << " -> n" << p->getIndex() << ":s" << destState << " p->getV()=" << p->getV() << std::endl;


	for (int i = 0; eventSum < g; i++) {
		//	eventSum = eventSum + ( tempR[srcState][destState] * ranPtr->poissonRv(muV) );
		double poissonVal = ranPtr->poissonProb(muV, i);
		//std::cout << "R[" << srcState << "][" << destState << "] = " << std::setw(10) << tempR[srcState][destState];
		//std::cout << "   *  Poisson(mu*v) = " << poissonVal;
		eventSum = eventSum + ( tempR[srcState][destState] * poissonVal );
		tempR = r * tempR;	// R^n is augmented after the calculation so as to begin with R^0 = I
		eventCount = i;

		//std::cout << "   = eventSum = " << std::setw(10) << eventSum;
		//std::cout << "   < g = " << g << "\n";


	}
	//std::cout << "event count: " << eventCount << "\n";

	// sample uniformly distributed event times, then sort
	for (int i = 0; i < eventCount; i++) {
		eventTimes.push_back(brlen * (double)rand() / (double)RAND_MAX);
	}
	sort(eventTimes.begin(), eventTimes.end());

	// helps make matrix exponentiation more efficient
	if (eventCount > eventCountMax) eventCountMax = eventCount;

	// sample the nature of the events in order, marginalized over their exact timing
	//      s_1 ~     p(s_1 = l | s_0 = a, s_n = b)     := [R]_al * [R^(n-1)]_lb
	//      s_2 ~     p(s_2 = m | s_1 = l, s_n = b)     := [R]_lm * [R^(n-2)]_mb
	//      ...
	//      s_k ~    p(s_k = q | s_(k-1) = r, s_n = b)  := [R]_rq * [R^(n-k)]_qb

	// initialize variables
	for (int i = 0; i < numStates; i++) {
		for (int j = 0; j < numStates; j++) {
			if (i != j) tempR[i][j] = 0.0;
		}
		tempR[i][i] = 1.0;
	}

	// perform additional matrix exponentiations as needed
	while(tempRExponents.size() < (unsigned int)eventCountMax) {
		tempRExponents.push_back(tempR);
		tempR = r * tempR;
	}

	p->getActivePath()->clearSitePathHistory(codonSite);
	p->getActivePath()->pushStep(Step(p->getV(), destState, codonSite), codonSite);

	//      add s_k = q
	for (int i = 0; i < eventCount; i++)
	{
		// initialize variables for event
		stepSum = 0.0;
		stepSample = 0.0;

		// calculate state space for sampling
		for (int j = 0; j < numStates; j++)
			stepSum += r[srcState][j] * tempRExponents[eventCount - i - 1][j][destState];

		// sample state from state space
		u = ((double)rand() / (double)RAND_MAX) * stepSum;
		for (int j = 0; j < numStates; j++)
		{
			stepSample += r[srcState][j] * tempRExponents[eventCount - i - 1][j][destState];
			if (stepSample > u)
			{
				if (r[srcState][j] == 0.0) std::cout << "ERROR: illegal state step committed!\n";
				srcState = j;
				j = numStates;       // break loop
			}
			else if (stepSample < u && srcState == numStates) std::cout << "stepSample < u*stepSum. state never found!" << std::endl;
		}

		// commit step to path (event & time)
		p->getActivePath()->addNextStep(Step(eventTimes[i], srcState, codonSite), codonSite);

		//std::cout << "event#: " << i << ", t: " << eventTimes[i] << ", state: " << srcState << "\n";
	}
	//std::cout << "\n";
}

bool ModelPath::isTransversion(int a, int b) {

	if (a == b)
		return false;
	int pA = (int)pow(2.0, (double)a);
	int pB = (int)pow(2.0, (double)b);
	if ( pA + pB == 5 || pA + pB == 10)
		return false;
	return true;
}

Freqs* ModelPath::getActiveFreqs(void) {

	for (std::vector<Parm*>::iterator p=parameters.begin(); p != parameters.end(); p++)
		{
		CodonFrequencies* derivedPtr = dynamic_cast<CodonFrequencies*> (*p);
		if (derivedPtr != 0)
			{
			return derivedPtr->getActiveFreqs();
			break;
			}
		}
	return NULL;
}

Dnds* ModelPath::getActiveDnds(void) {

	for (std::vector<Parm*>::iterator p=parameters.begin(); p != parameters.end(); p++)
	{
		Omega* derivedPtr = dynamic_cast<Omega*> (*p);
		if (derivedPtr != 0)
		{
			return derivedPtr->getActiveDnds();
			break;
		}
	}
	return NULL;
}

TiTv* ModelPath::getActiveTiTv(void) {

	for (std::vector<Parm*>::iterator p=parameters.begin(); p != parameters.end(); p++)
	{
		Kappa* derivedPtr = dynamic_cast<Kappa*> (*p);
		if (derivedPtr != 0)
		{
			return derivedPtr->getActiveTiTv();
			break;
		}
	}
	return NULL;
}

Topology* ModelPath::getActiveTopology(void) {

	for (std::vector<Parm*>::iterator p=parameters.begin(); p != parameters.end(); p++)
	{
		Tree* derivedPtr = dynamic_cast<Tree*> (*p);
		if (derivedPtr != 0)
		{
			return derivedPtr->getActiveTopology();
			break;
		}
	}
	return NULL;
}

Tree* ModelPath::getActiveTree(void) {

	for (std::vector<Parm*>::iterator p=parameters.begin(); p != parameters.end(); p++)
	{
		Tree* derivedPtr = dynamic_cast<Tree*> (*p);
		if (derivedPtr != 0)
		{
			return derivedPtr;
			break;
		}
	}
	return NULL;
}

Parm* ModelPath::pickParmAtRandom(void) {

	double u = ranPtr->uniformRv();
	double sum = 0.0;
	for (int i=0; i < (int)proposalProbs.size(); i++) {
		sum += proposalProbs[i];
		if (u < sum)
			return parameters[i];
	}
	return NULL;
}

void ModelPath::updateAllFlags(void) {

	Topology* t = getActiveTopology();
	t->updateAllCls(true);
	t->updateAllTis(true);
	t->flipAllActiveCls();
	t->flipAllActiveTis();
	tiProbs->setTransitionProbs();
	t->updateAllTis(true);
}

void ModelPath::updateQ(Parm* p) {

	bool shouldUpdate = false;
	{
		CodonFrequencies* derivedPtr = dynamic_cast<CodonFrequencies*> (p);
		if (derivedPtr != 0)
		{
			shouldUpdate = true;
			goto funcEnd;
		}
	}

	{
		Kappa* derivedPtr = dynamic_cast<Kappa*> (p);
		if (derivedPtr != 0)
		{
			shouldUpdate = true;
			goto funcEnd;
		}
	}

	{
		Omega* derivedPtr = dynamic_cast<Omega*> (p);
		if (derivedPtr != 0)
		{
			shouldUpdate = true;
			goto funcEnd;
		}
	}

	funcEnd:
		if (shouldUpdate == true)
		{
			fillInRateMatrix();
			tiMatrix->updateQ(q);
		}
}

void ModelPath::restoreQ(Parm* p) {

	bool shouldUpdate = false;
	{
		CodonFrequencies* derivedPtr = dynamic_cast<CodonFrequencies*> (p);
		if (derivedPtr != 0)
			{
			shouldUpdate = true;
			goto funcEnd;
			}
	}

	{
		Kappa* derivedPtr = dynamic_cast<Kappa*> (p);
		if (derivedPtr != 0)
			{
			shouldUpdate = true;
			goto funcEnd;
		}
	}

	{
		Omega* derivedPtr = dynamic_cast<Omega*> (p);
		if (derivedPtr != 0)
			{
			shouldUpdate = true;
			goto funcEnd;
		}
	}

	funcEnd:
		if (shouldUpdate == true)
		{
			tiMatrix->restoreQ();
		}
}

void ModelPath::keep(Parm* p) {

	p->keep();

	Tree* derivedPtr = dynamic_cast<Tree*> (p);
	if (derivedPtr == 0)
		{
		getActiveTree()->keep();
		}
}

void ModelPath::restore(Parm* p) {

	p->restore();

	Tree* derivedPtr = dynamic_cast<Tree*> (p);
	if (derivedPtr == 0)
	{
		getActiveTree()->restore();
	}

}

void ModelPath::printQ(void) {

	std::cout << "ModelPath: Q MATRIX\n";
	for (int i=0; i<numStates; i++)
	{
		for (int j=0; j<numStates; j++)
		{
			if (i == j)
				std::cout << std::fixed << std::setprecision(8) << q[i][j] << " ";
			else
				std::cout << std::fixed << std::setprecision(8) << " " << q[i][j] << " ";

		}
		std::cout << "\n";
	}
	std::cout << "\n";
}

void ModelPath::printR(void) {

	std::cout << "ModelPath: R MATRIX\n";
	for (int i=0; i<numStates; i++)
	{
		for (int j=0; j<numStates; j++)
		{
			std::cout << std::fixed << std::setprecision(8) << r[i][j] << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
}


/*
void ModelPath::initializeCondLikes(void) {

	// get a pointer to the current tree topology
	Topology* t = getActiveTopology();

	// update the transition probabilities, if necessary
	tiProbs->setTransitionProbs();
	//stochasticTiProbs->setTransitionProbs();

	// allocate some matrices to temporarily hold the transition matrices
	MbMatrix<double> tiL;
	MbMatrix<double> tiR;
	MbMatrix<double> tiA;

	// calculate the conditional likelihoods down the tree
	for (int n = 0; n < t->getNumNodes(); n++)
	{
		Node* p = t->getDownPassNode(n);
		if (p->getLft() != NULL && p->getRht() != NULL && p->getAnc() != NULL)
		{
#			if defined(RESCALE_LIKES)
			double* clS = condLikes->getClScalerPtr( p->getActiveCl(), p->getIndex(), 0 );
#			endif
			int lftIdx = p->getLft()->getIndex();
			int rhtIdx = p->getRht()->getIndex();
			tiL = tiProbs->getTiMatrix( p->getLft()->getActiveTi(), lftIdx );
			tiR = tiProbs->getTiMatrix( p->getRht()->getActiveTi(), rhtIdx );
			double* clL = condLikes->getClPtr( p->getLft()->getActiveCl(), lftIdx, 0 );
			double* clR = condLikes->getClPtr( p->getRht()->getActiveCl(), rhtIdx, 0 );
			double* clP = condLikes->getClPtr( p->getActiveCl(), p->getIndex(), 0 );
			for (int c = 0; c < numSites; c++)
			{
#				if defined(RESCALE_LIKES)
				double maxCl = 0.0;
#				endif
				for (int i = 0; i < numStates; i++)
				{
					double sumL = 0.0;
					double sumR = 0.0;
					for (int j = 0; j < numStates; j++)
					{
						sumL += tiL[i][j] * clL[j];
						sumR += tiR[i][j] * clR[j];
					}
					clP[i] = sumL * sumR;
#					if defined(RESCALE_LIKES)
					if (clP[i] > maxCl)
						maxCl = clP[i];
#					endif
				}
#				if defined(RESCALE_LIKES)
				double scaler = 1.0 / maxCl;
				for (int i = 0; i < numStates; i++)
					clP[i] *= scaler;
				clS[c] = log(maxCl);
#				endif
				clL += numStates;
				clR += numStates;
				clP += numStates;
			}
		}
	}

	// calculate the conditional likelihoods up the tree
	for (int n = t->getNumNodes() - 1; n >= 0; n--)
	{
		Node* p = t->getDownPassNode(n);
		if (p->getLft() != NULL && p->getRht() != NULL && p->getAnc() != NULL)
		{
			bool pToLeft = true;
			if (p->getAnc()->getLft() != p)
				pToLeft = false;
#			if defined(RESCALE_LIKES)
			double* clSUpL = condLikes->getClScalerPtrUpL( p->getActiveCl(), p->getIndex(), 0 );
			double* clSUpR = condLikes->getClScalerPtrUpR( p->getActiveCl(), p->getIndex(), 0 );
#			endif
			int lftIdx = p->getLft()->getIndex();
			int rhtIdx = p->getRht()->getIndex();
			int ancIdx = p->getAnc()->getIndex();
			tiL = tiProbs->getTiMatrix( p->getLft()->getActiveTi(), lftIdx );
			tiR = tiProbs->getTiMatrix( p->getRht()->getActiveTi(), rhtIdx );
			tiA = tiProbs->getTiMatrix( p->getActiveTi(),           p->getIndex() );
			double* clL = condLikes->getClPtr( p->getLft()->getActiveCl(), lftIdx, 0 );
			double* clR = condLikes->getClPtr( p->getRht()->getActiveCl(), rhtIdx, 0 );
			double* clA;
			if (pToLeft == true)
				clA = condLikes->getClPtrUpR( p->getAnc()->getActiveCl(), ancIdx, 0 );
			else
				clA = condLikes->getClPtrUpL( p->getAnc()->getActiveCl(), ancIdx, 0 );
			double* clPUpL = condLikes->getClPtrUpL( p->getActiveCl(), p->getIndex(), 0 );
			double* clPUpR = condLikes->getClPtrUpR( p->getActiveCl(), p->getIndex(), 0 );
			for (int c=0; c<numSites; c++)
			{
#				if defined(RESCALE_LIKES)
				double maxCl = 0.0;
#				endif
				for (int i=0; i<numStates; i++)
				{
					double sumL = 0.0, sumA = 0.0;
					for (int j=0; j<numStates; j++)
					{
						sumL += tiL[i][j] * clL[j];
						sumA += tiA[i][j] * clA[j];
					}
					clPUpL[i] = sumL * sumA;
#					if defined(RESCALE_LIKES)
					if (clPUpL[i] > maxCl)
						maxCl = clPUpL[i];
#					endif
				}
#				if defined(RESCALE_LIKES)
				double scaler = 1.0 / maxCl;
				for (int i=0; i<numStates; i++)
					clPUpL[i] *= scaler;
				clSUpL[c] = log(maxCl);
#				endif

#				if defined(RESCALE_LIKES)
				maxCl = 0.0;
#				endif
				for (int i=0; i<numStates; i++)
				{
					double sumR = 0.0, sumA = 0.0;
					for (int j=0; j<numStates; j++)
					{
						sumR += tiR[i][j] * clR[j];
						sumA += tiA[i][j] * clA[j];
					}
					clPUpR[i] = sumR * sumA;
#					if defined(RESCALE_LIKES)
					if (clPUpR[i] > maxCl)
						maxCl = clPUpR[i];
#					endif
				}
#				if defined(RESCALE_LIKES)
				scaler = 1.0 / maxCl;
				for (int i=0; i<numStates; i++)
					clPUpR[i] *= scaler;
				clSUpR[c] = log(maxCl);
#				endif

				clL += numStates;
				clR += numStates;
				clA += numStates;
				clPUpL += numStates;
				clPUpR += numStates;
			}
		}
	}
}
*/

/*
void ModelPath::samplePathToNode(Node* p) {

	// TODO: THIS METHOD IS NO LONGER IN USE

	// STEP TWO of three : sample series of events between nodes for each site

	// get a pointer to the current tree topology
	//Topology* t = getActiveTopology();

	// initialize variables
	MbMatrix<double> tempR = MbMatrix<double>(numStates, numStates, 0.0);
	MbMatrix<double> tiP;
	std::vector<double> eventTimes;
	double brlen = 0.0;
	double muV = 0.0;
	double g = 0.0;
	double u = 0.0;
	double eventSum = 0.0;
	int eventCount = 0;
	int eventCountMax = 0;
	int codonSite = 0;
	int srcState;
	int destState;

	// list iterator to get states
//	std::list<Step> srcTempList = p->getAnc()->getPath()->getEventHistory();
//	std::list<Step> destTempList = p->getPath()->getEventHistory();
	std::list<Step>::iterator srcIt = srcTempList.begin();
	std::list<Step>::iterator destIt = destTempList.begin();

	// 2.1: sample the number of events (always including virtual events) marginalized over their nature and timing)
	for (int n = 0; n < numCodonSites; n++)
	{
		// initialize variables
		for (int i = 0; i < numStates; i++) {
			for (int j = 0; j < numStates; j++) {
				if (i != j) tempR[i][j] = 0.0;
			}
			tempR[i][i] = 1.0;
		}
		eventCount = 0;
		eventSum = 0.0;
		eventTimes.clear();

		// update variables
		brlen = p->getV();
		srcState = (*srcIt).getState();
		destState = (*destIt).getState();
		muV = uniformizationRate * brlen;
		tiP = tiProbs->getTiMatrix( p->getActiveTi(), p->getIndex() );
		g = tiP[srcState][destState] * (double)rand() / (double)RAND_MAX;

		//std::cout << n << ": " << "path n" << p->getAnc()->getIndex() << ":s" << srcState;
		//std::cout << " -> n" << p->getIndex() << ":s" << destState << " p->getV()=" << p->getV() << std::endl;

		for (int i = 0; eventSum < g; i++) {
			eventSum = eventSum + ( tempR[srcState][destState] * ranPtr->poissonRv(muV) );
			tempR = r * tempR;	// R^n is incremented after the calculation so as to begin with R^0 = I
			eventCount = i;
		}

		// update event count for path
		// std::cout << "eventCount = " << eventCount << std::endl;
		// p->setEventCount(eventCount);

		// 2.3 sample event times for all events (necessarily precedes step 2.2)
		for (int i = 0; i < eventCount; i++) {
			////eventTimes.push_back(brlen * ranPtr->uniformRv());
			eventTimes.push_back(brlen * (double)rand() / (double)RAND_MAX);
		}

		//std::cout << "Sorting eventTimes" << std::endl;
		sort(eventTimes.begin(), eventTimes.end());
		//reverse(eventTimes.begin(), eventTimes.end());

		//std::cout << "site: " << n << "\n";
		//for (std::vector<double>::iterator it = eventTimes.begin(); it != eventTimes.end(); it++) std::cout << (*it) << std::endl;

		// helps make matrix exponentiation more efficient
		if (eventCount > eventCountMax) eventCountMax = eventCount;


		// 2.2 sample the nature of the events in order, marginalized over their exact timing
		//      s_1 ~     p(s_1 = l | s_0 = a, s_n = b)     := [R]_al * [R^(n-1)]_lb
		//      s_2 ~     p(s_2 = m | s_1 = l, s_n = b)     := [R]_lm * [R^(n-2)]_mb
		//      ...
		//      s_k ~    p(s_k = q | s_(k-1) = r, s_n = b) := [R]_rq * [R^(n-k)]_qb

		// initialize variables
		for (int i = 0; i < numStates; i++) {
			for (int j = 0; j < numStates; j++) {
				if (i != j) tempR[i][j] = 0.0;
			}
			tempR[i][i] = 1.0;
		}
		double stepSum = 0.0;
		double stepSample = 0.0;

		// update variables
		//Step *currStepPtr;
		Step currStep;

		std::list<Step> tempList;
		std::vector< MbMatrix<double> > tempRExponents;

		// perform additional matrix exponentiations as needed
		while(tempRExponents.size() < (unsigned int)eventCountMax) {
			tempRExponents.push_back(tempR);
			tempR = r * tempR;
		}

		// TODO: probably need to be more careful about object instantiation (new & delete) in order to create linked steps...
		// create a currStepPtr.  create NEW nextStepPtr as necessary, so the # of stepPtrs increases as the loop iterates

		//      add s_0 = a
		//nextStepPtr = new Step();
		//currStepPtr = new Step(0.0, srcState, codonSite, nextStepPtr);
		//std::cout << "s_0 " << currStepPtr << std::endl;
		//std::cout << "s_k " << nextStepPtr << std::endl;
		//std::cout << "s_0 -> next " << currStepPtr->getNextStep() << "\n\n";

		//tempList.push_back(*currStepPtr);
		//// tempList.push_back(Step(0.0, srcState, codonSite));



		//      add s_k = q
		for (int i = 0; i < eventCount; i++)
		{
			// initialize variables for event
			stepSum = 0.0;
			stepSample = 0.0;
			//currStepPtr = nextStepPtr;
			//nextStepPtr = new Step();

			// calculate state space for sampling
			for (int j = 0; j < numStates; j++)
				stepSum += r[srcState][j] * tempRExponents[eventCount - i - 1][j][destState];

			// sample state from state space
			////u = ranPtr->uniformRv() * stepSum;
			u = ((double)rand() / (double)RAND_MAX) * stepSum;
			for (int j = 0; j < numStates; j++)
			{
				stepSample += r[srcState][j] * tempRExponents[eventCount - i - 1][j][destState];
				if (stepSample > u)
				{
					if (r[srcState][j] == 0.0) std::cout << "ERROR: illegal state step committed!\n";
					srcState = j;
					j = numStates;       // break loop
				}
				else if (stepSample < u && srcState == numStates) std::cout << "stepSample < u*stepSum. state never found!" << std::endl;
			}

			// commit step to path (event & time)
			// std::cout << "  step# " << i+1 << ", state# " << currStep << ", t: " << eventTimes[i] << std::endl;

			//currStepPtr = Step(eventTimes[i], srcState, codonSite);
			currStep = Step(eventTimes[i], srcState, codonSite);
			//currStepPtr = &(Step(eventTimes[i], srcState, codonSite, nextStepPtr));
			tempList.push_back(currStep);
			//tempList.push_back(*currStepPtr);
		}

		// update the path
		p->getPath()->pushEventHistory(tempList);

		srcIt++;
		destIt++;
		codonSite++;
	}
}
*/


#if 0
// GRAVEYARD
void ModelPath::sampleNodeStates() {

	// get a pointer to the current tree topology
	Topology* t = getActiveTopology();

	// sample the root state
	Node* r = t->getRoot();
	double u = 0.0;

	Step tempStep(0.0, 0, 0);
	std::vector<Step> tempVector;

	for (int codonSite = 0; codonSite < numCodonSites; codonSite++)
	{
		double rootDemSum = 0.0;
		double rootNumSum = 0.0;
		double* clRoot = condLikes->getClPtr( r->getActiveCl(), r->getIndex(), codonSite);
		u = (double)rand() / (double)RAND_MAX;

		for (int i = 0; i < numStates; i++) {
			// denominator sum
			rootDemSum += clRoot[i] * stationaryProbs[i];
		}

		for (int i = 0; i < numStates; i++) {
			// numerator sum
			rootNumSum += clRoot[i] * stationaryProbs[i];
			if (u < (rootNumSum / rootDemSum)) {
				tempStep = Step(0.0, i, codonSite);
				////r->getPath()->pushStep(tempStep, codonSite);
				r->getActivePath()->pushStep(tempStep, codonSite);
				break;
			}
			else if (u > rootNumSum / rootDemSum && i == numStates)
			{
				std::cout << "no new state set!\n";
			}
		}
	}

	// validate root states are sampled correctly
//	r->getActivePath()->printPathHistory();

	// sample all node states, traveling from root to tips
//	int ancState = -1;
//	double pDenSum = 0.0;
//	double pNumSum = 0.0;
//	double* clP = NULL;
	MbMatrix<double> tiP;

	for (int n = t->getNumNodes() - 1; n >= 0; n--)
	{
		Node* p = t->getDownPassNode(n);
	//	if (p->getAnc() != NULL)
		//{
			////std::vector<std::vector<Step> > ancPath = p->getAnc()->getPath()->getPathHistory();
//			std::vector<std::vector<Step> > ancPath = p->getAnc()->getActivePath()->getPathHistory();
			for (int codonSite = 0; codonSite < numCodonSites; codonSite++)
			{
			sampleNodeSiteState(p, codonSite);
				/*
				u = (double)rand() / (double)RAND_MAX;
				pDenSum = 0.0;
				pNumSum = 0.0;

				ancState = ancPath[codonSite][0].getState();
				clP = condLikes->getClPtr( p->getActiveCl(), p->getIndex(), codonSite);
				tiP = tiProbs->getTiMatrix( p->getActiveTi(), p->getIndex() );

				// denominator
				// sum of cl(root-1, k) * P(root-1, i->k)
				for (int k = 0; k < numStates; k++)
				{
					pDenSum += clP[k] * tiP[ancState][k];
				}

				// numerator
				// cl(root-1, j) * P(root-1, i->j)
				for (int j = 0; j < numStates; j++)
				{
					pNumSum += clP[j] * tiP[ancState][j];
					if (u < (pNumSum / pDenSum))
					{
						tempStep = Step(p->getV(), j, codonSite);
						p->getActivePath()->pushStep(tempStep, codonSite);
						break;
					}
				}
				*/
			}

		//}
	}
}
#endif



/*
void ModelPath::sampleRootState() {

	// STEP ONE of three : sample states at internal nodes
	//
	//        accomplished by first sampling the state at the root node according to the product of the
	//        stationary probability and the conditional likelihood of each state, followed by a
	//        recursive sampling along the tree, with each case conditioned on the state at the
	//        immediate ancestral node (Nielsen, 2002)

	// random number generator


	// get a pointer to the current tree topology
	Topology* t = getActiveTopology();

	// sample the root state
	Node* r = t->getRoot();
	double u = 0.0;

	Step tempStep(0.0, 0, 0);
	std::vector<Step> tempVector;

	for (int codonSite = 0; codonSite < numCodonSites; codonSite++)
	{
		double rootDemSum = 0.0;
		double rootNumSum = 0.0;
		double* clRoot = condLikes->getClPtr( r->getActiveCl(), r->getIndex(), codonSite);
		u = (double)rand() / (double)RAND_MAX;

		for (int i = 0; i < numStates; i++) {
			// denominator sum
			rootDemSum += clRoot[i] * stationaryProbs[i];
		}

		for (int i = 0; i < numStates; i++) {
			// numerator sum
			rootNumSum += clRoot[i] * stationaryProbs[i];
			if (u < (rootNumSum / rootDemSum)) {
				tempStep = Step(0.0, i, codonSite);
				////r->getPath()->pushStep(tempStep, codonSite);
				r->getActivePath()->pushStep(tempStep, codonSite);
				break;
			}
		}
	}
}
*/
