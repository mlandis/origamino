/*
 * ModelSeq.cpp
 *
 *  Created on: Jan 24, 2010
 *      Author: mikee
 */

#include "ModelSeq.h"

#define ZERO_ENTRY		0
#define SYN_TI_ENTRY	1
#define SYN_TV_ENTRY	2
#define NONSYN_TI_ENTRY	3
#define NONSYN_TV_ENTRY	4

#define RESCALE_LIKES

ModelSeq::ModelSeq(Alignment *ap, MbRandom *rp, Settings* sp, Tree* tp) {

	std::cout << "INITIALIZING MODEL: SEQUENCE-STRUCTURE FIT\n\n";

	// keep track of pointers to important objects
	alignmentPtr	= ap;
	ranPtr			= rp;
	treePtr			= tp;
	numCodonSites	= alignmentPtr->getNumCodonSites();
	numSampleSeqs	= 100;
	contactDistance = sp->getContactDistance();
	SVESamplesPerCycle = sp->getSVESamplesPerCycle();
	brlenLambda		= sp->getBrlenLambda();

	// allocate a vector holding the stationary probabilities
	stationaryProbs = new double[numCodonSites];

	// add the parameters to the phylogenetic model
	parameters.push_back(treePtr);
	parameters.push_back(new CodonFrequencies(ranPtr, "Codon Frequencies", treePtr, alignmentPtr));
	parameters.push_back(new Kappa(ranPtr, treePtr, "Kappa"));
	parameters.push_back(new Omega(ranPtr, treePtr, "Omega"));
	parameters.push_back(new ScalerU(ranPtr, treePtr, brlenLambda, "SubstRate"));
	parameters.push_back(new SeqStructS( ranPtr, treePtr, "SeqStructS"));
	parameters.push_back(new SeqStructP( ranPtr, treePtr, "SeqStructP"));

	// initialize parametersGibbs
	// parametersGibbs = parameters;

	std::cout << "\n";

	// initialize the proposal probabilities
	proposalProbs.push_back( 0.0 ); // probability of proposing a change to the tree
	proposalProbs.push_back( 1.0 ); // probability of proposing a change to the codon frequencies
	proposalProbs.push_back( 1.0 ); // probability of proposing a change to the transition/transversion rate ratio
	proposalProbs.push_back( 1.0 ); // probability of proposing a change to the dN/dS rate ratio
	proposalProbs.push_back( 5.0 ); // probability of proposing a change to the scaling parameters
	proposalProbs.push_back( 1.0 ); // probability of proposing a change to the solvent accessibility
	proposalProbs.push_back( 1.0 ); // probability of proposing a change to the pairwise structure (contact map?)

	double sum = 0.0;
	for (int i = 0; i < (int)proposalProbs.size(); i++)
		sum += proposalProbs[i];
	for (int i = 0; i < (int)proposalProbs.size(); i++)
		proposalProbs[i] /= sum;

	// set up the qTemplate rate matrix (which runs fillInRateMatrix() in the process
	initializeRateMatrix(alignmentPtr->getNumSenseCodons());
	initializeContactMap(sp->getCoordinatesFileName(), numCodonSites, contactDistance);
	initSVE();

	// instantiate the conditional likelihoods
	condLikes = new CondLikes(alignmentPtr);

	// instantiate the transition probabilities

	////tiMatrix = new MbTransitionMatrix(q, false);
	////tiProbs = new TiProbs( 2 * alignmentPtr->getNumTaxa() - 2, tiMatrix, treePtr);

#	if defined(RESCALE_LIKES)
	// allocate a vector that holds (temporarily) the log scalers for log likelihood calculations
	lnScalers = new double[numCodonSites];
#	endif

	// print the parameters
	//for (std::vector<Parm*>::iterator p=parameters.begin(); p != parameters.end(); p++)
	//	(*p)->print();
	std::cout << "\n";
}

ModelSeq::~ModelSeq(void) {

	delete condLikes;
	delete tiProbs;
	for (std::vector<Parm*>::iterator p=parameters.begin(); p != parameters.end(); p++)
		delete (*p);
	delete [] stationaryProbs;
#	if defined(RESCALE_LIKES)
	delete [] lnScalers;
#	endif
}

void ModelSeq::initializeRateMatrix(int nc) {

	numSenseCodons		= nc;
	numStates 			= nc;

	q = MbMatrix<double>(numSenseCodons, numSenseCodons);
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

}

void ModelSeq::fillInRateMatrix(void) {

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

	// calculate R_v,# values (rate away from v)
	double rateSum = 0.0;
	for (int i = 0; i < numSenseCodons; i++) {
		for (int j = 0; j < numSenseCodons; j++) {
			if (i != j) rateSum += q[i][j];
		}
		qRateAway.push_back(rateSum);
		rateSum = 0.0;
	}



}

void ModelSeq::fillInStationaryProbs(void) {

	std::vector<double> f = getActiveFreqs()->getFreqs();
	double p = getActiveDnds()->getDnds();
	double *x = &stationaryProbs[0];
	for (int i=0; i<numSenseCodons; i++)
	{
		(*x) = f[i] * p;
		x++;
	}
}


void ModelSeq::initializeContactMap(std::string pdbFilePath, int ncs, double contactDistance) {

	numCodonSites = ncs;

	FileMgr coordMatrixFileMgr(pdbFilePath);
	std::ifstream seqStream;

	if (coordMatrixFileMgr.openFile(seqStream) == false) {
		std::cerr << "Cannot open file \"" + coordMatrixFileMgr.getFileName() + "\"" << std::endl;
		exit(1);
	}

	// temp variables
	std::string lineStr = "";
	std::string name = "";
	std::string atom = "";
	std::string residue = "";

	std::vector<std::vector<double> > xCoordSeq(numCodonSites);
	std::vector<std::vector<double> > yCoordSeq(numCodonSites);
	std::vector<std::vector<double> > zCoordSeq(numCodonSites);

	// initialize contactMap
	contactMap.resize(numCodonSites);
	for (int i = 0; i < numCodonSites; i++) {
		contactMap[i].resize(numCodonSites, 0);
	}

	// read values
	std::cout << "PARSING PDB FILE: " << pdbFilePath << "\n";
	while (getline(seqStream, lineStr).good())
	{
		name = lineStr.substr(0,6);
		atom = lineStr.substr(12,4);
		residue = lineStr.substr(22,4);
		int residuePos = atoi(residue.c_str());

		//std::cout << "name: " << name << ", atom: " << atom << ", residue:" << residue << "\n";

		if (name == "ATOM  " && atom != " H  " && residuePos < numCodonSites)
		{
			std::string xstr = lineStr.substr(30,8);
			std::string ystr = lineStr.substr(38,8);
			std::string zstr = lineStr.substr(46,8);
		//	std::cout << "xstr: " << xstr << ", ystr: " << ystr  << ", zstr: " << zstr << "\n";
		//	std::cout << "xstr: " << atof(xstr.c_str()) << ", ystr: " << atof(ystr.c_str()) << ", zstr: " << atof(zstr.c_str()) << "\n";
			xCoordSeq[residuePos].push_back(atof(xstr.c_str()));
			yCoordSeq[residuePos].push_back(atof(ystr.c_str()));
			zCoordSeq[residuePos].push_back(atof(zstr.c_str()));
		}
	}

	// compute contact map using x^2 + y^2 + z^2 = d^2
	std::cout << "CONSTRUCTING CONTACT MAP\n";
	for (int ai = 0; ai < numCodonSites; ai++)
	{
		for (unsigned int aj = 0; aj < xCoordSeq[ai].size(); aj++) {
			double xaVal = xCoordSeq[ai][aj];
			double yaVal = yCoordSeq[ai][aj];
			double zaVal = zCoordSeq[ai][aj];


			for (int bi = ai; bi < numCodonSites; bi++)
			{
				if (ai < bi - 1 || ai > bi + 1)
				{
					for (unsigned int bj = 0; bj < xCoordSeq[bi].size(); bj++)
					{
						double xbVal = xCoordSeq[bi][bj];
						double ybVal = yCoordSeq[bi][bj];
						double zbVal = zCoordSeq[bi][bj];

						// Euclidean distance
						if (pow(pow(xaVal - xbVal, 2) + pow(yaVal - ybVal, 2) + pow(zaVal - zbVal, 2), 0.5) < contactDistance)
						{
							contactMap[ai][bi]++;// = 1;
							contactMap[bi][ai]++;// = 1;
						}
					}
				}
			}
		}
	}

	std::cout << "ModelSeq: CONTACT MAP\n";
	for (std::vector<std::vector<int> >::iterator iti = contactMap.begin(); iti != contactMap.end(); iti++)
	{
		for (std::vector<int>::iterator itj = (*iti).begin(); itj != (*iti).end(); itj++)
		{
			if (*itj == 0)
				std::cout << ".  ";
			else
				std::cout << std::setw(2) << *itj << " ";
		}
		std::cout << "\n";
	}
}

double ModelSeq::calcPseudoEnergy(std::vector<Step> seq)
{
	// TODO: populate this method with eqn from Kleinman et al 2006
	return 0.0;
}

bool ModelSeq::isTransversion(int a, int b) {

	if (a == b)
		return false;
	int pA = (int)pow(2.0, (double)a);
	int pB = (int)pow(2.0, (double)b);
	if ( pA + pB == 5 || pA + pB == 10)
		return false;
	return true;
}

double ModelSeq::lnPathLikelihood(void) {

	//
	// this term is equals pr(j,rho|i,theta) * pr(i|theta)
	// 		where
	// 			pr(j,rho|i,theta) 	= [*=:events, R_i(z-1),i(z) * e^(-R_i(z-1),# * (t(z) - t(z-1))] * e^(-R_i(q),# * (t(q+1) - t(q)))
	//		and
	// 			pr(i|theta) 		= e^(-2*s*E_s(i)-2*p*E_p(i) * [*=1:N(pi_i__m)] /
	//		      			          [+=1:k(e^(-2*s*E_s(k)-2*p*E_p(k) * [*=1:N(pi_k__n)]]
	//
	// 		note: the pr(i,j) cancels out in the Metropolis-Hastings ratio
	//
	// pr(D|rho,theta,beta) = pr(s^(root)|theta)*[*=(2N-2, b=1), pr(s^(b->down),rho_b|s^(b->up),beta_b,theta)]
	// 		where: rho = paths, theta = params, beta = topology
	//
	//
	// get a pointer to the current tree topology
	Topology* t = getActiveTopology();

	// update the transition probabilities, if necessary
	//// tiProbs->setTransitionProbs();

	double unitBrlen = 0.0;
	double unitTimeMarginOfError = 0.0000000000000005;

	int currState = 0, prevState = 0;
	double unitCurrTime = -1.0, unitPrevTime = 0.0, unitFinalTime = 0.0;
	double lnStepProbSum = 0.0; //, stepProbSum = 1.0;
	double rateAwaySum = 0.0;
	double lnStepProb = 0.0, stepProb = 0.0;
	double x = 0.0;
	std::vector<std::vector<Step> > srcPath;
	std::vector<std::vector<Step> > destPath;
	std::vector<int> prevSeq(numCodonSites);
	std::vector<double> u = getActiveSubstRate()->getSubstRate();

	for (int n = 0; n < t->getNumNodes(); n++)
	{
		Node* p = t->getDownPassNode(n);
		if (p->getAnc() != NULL)
		{
			// FIND		pr(j,rho|i,theta) * pr(i|theta)
			// 		where
			// 			pr(j,rho|i,theta) 	= [*=:events, R_i(z-1),i(z) * e^(-R_i(z-1),# * (t(z) - t(z-1))] * e^(-R_i(q),# * (t(q+1) - t(q)))
			//		and
			// 			pr(i|theta) 		= e^(-2*s*E_s(i)-2*p*E_p(i) * [*=1:N(pi_i__m)] /
			//		      			          [+=1:k(e^(-2*s*E_s(k)-2*p*E_p(k) * [*=1:N(pi_k__n)]]


			// pr(j,rho|i,theta) 	= [*=:events, R_i(z-1),i(z) * e^(-R_i(z-1),# * (t(z) - t(z-1))] * e^(-R_i(q),# * (t(q+1) - t(q)))

			unitBrlen = u[n];
			unitCurrTime = -1.0;
			unitPrevTime = 0.0;
			unitFinalTime = 0.0;
			rateAwaySum = 0.0;

			srcPath = p->getAnc()->getActivePath()->getPathHistory();
			destPath = p->getActivePath()->getPathHistory();


			//std::cout << "\n";
			//std::cout << "path describing n" << p->getAnc()->getIndex() << "-> n" << p->getIndex() << "   (u[" << n << "]: " << unitBrlen << ")\n";


			// intialize eventVector
			std::vector<Step> eventVector;
			for (std::vector<std::vector<Step> >::iterator it1 = destPath.begin(); it1 != destPath.end(); it1++)
			{
				for (std::vector<Step>::iterator it2 = (*it1).begin(); it2 != (*it1).end(); it2++)
				{
					eventVector.push_back(*it2);
					//std::cout << std::setw(2) << (*it2).getState() << " ";
				}
			}
			sort(eventVector.begin(), eventVector.end());
			reverse(eventVector.begin(), eventVector.end());

			// intialize prevSeq
			prevSeq.clear();
			for (std::vector<std::vector<Step> >::iterator it1 = srcPath.begin(); it1 != srcPath.end(); it1++)
			{
				Step s = *(*it1).begin();
				prevSeq[s.getSite()] = s.getState();
			}

			// calculate likelihood for events
			for (std::vector<Step>::iterator it = eventVector.begin(); it != eventVector.end(); it++)
			{
				unitCurrTime = (*it).getTime() * unitBrlen;
				if (1.0 - unitTimeMarginOfError < unitCurrTime && 1.0 + unitTimeMarginOfError > unitCurrTime) unitCurrTime = 1.0;
				if (unitCurrTime < 1.0 - unitTimeMarginOfError) unitFinalTime = unitCurrTime;

				// if the scaled event time == 1.0, it is a destination event
				if (unitCurrTime == 1.0) {
					//prevSeq[(*it).getSite()] = (*it).getState();
					unitPrevTime = unitCurrTime;
				}

				// otherwise, calculate the probability of the event occuring
				else if (unitPrevTime != unitCurrTime && unitCurrTime != 1.0)
				{
					rateAwaySum = 0.0;
					//stepProb = 1.0;
					lnStepProb = 0.0;

					prevState = prevSeq[(*it).getSite()];
					currState = (*it).getState();

					if (q[prevState][currState] == 0.0)
					{
						std::cout << "WARNING: q[" << std::setw(2) << prevState << "][" << std::setw(2) << currState << "] == 0.0\n";
						std::exit(1);
					}

					// do not consider virtual events
					//    or are these negative q[i][i] values that will improve likelihood
					if (prevState != currState) {
						// get the sum of rates for all sequences away from prevSeq
						// (numCodonSites) * (numSenseCodons - 1) sequences differ by one nucleotide
						for (int i = 0; i < numCodonSites; i++) {
							for (int j = 0; j < numSenseCodons; j++) {
								if (prevSeq[i] != j) rateAwaySum += q[prevSeq[i]][j];
							}
						}
						rateAwaySum *= unitBrlen;

						stepProb = q[prevState][currState] * unitBrlen;
						lnStepProb = log(stepProb);

						x = -rateAwaySum * (unitCurrTime - unitPrevTime);
						lnStepProb += x;

						/*
						std::cout << prevState << "->" << currState << ":s" << (*it).getSite();
						std::cout << ", lnstepProb: "			<< lnStepProb;
						std::cout << ", R_ij: "				<< q[prevState][currState] * unitBrlen;
						std::cout << ", rateAwaySum: "		<< rateAwaySum;
						std::cout << ", dt: "				<< unitCurrTime - unitPrevTime;
						std::cout << ", -R_i* * dt: "		<< -rateAwaySum * unitBrlen * (unitCurrTime - unitPrevTime);
						std::cout << ", e^(-R_i* * dt): "	<<  exp(-rateAwaySum * unitBrlen * (unitCurrTime - unitPrevTime));
						std::cout << "\n";
						*/

						lnStepProbSum += lnStepProb;
						//std::cout << "lnStepProbSum: " << lnStepProbSum << "\n\n";

						prevSeq[(*it).getSite()] = currState;
						unitCurrTime = unitPrevTime;
					}
				}

				else if (unitPrevTime == unitCurrTime && unitCurrTime != 1.0)
					std::cout << "ERROR: duplicate event times, prevTime: " << unitPrevTime << ", currTime: " << unitCurrTime << "\n";

			}

			// calculate final term
			rateAwaySum = 0.0;
			for (int i = 0; i < numCodonSites; i++) {
				for (int j = 0; j < numSenseCodons; j++) {
					if (prevSeq[i] != j) rateAwaySum += q[prevSeq[i]][j];
				}
			}
			x = -rateAwaySum * unitBrlen * (1.0 - unitFinalTime);
			lnStepProbSum += x;
			//std::cout << "lnStepProbSum: " << lnStepProbSum << "\n";
		}
	}

	return lnStepProbSum;
}


double ModelSeq::lnRootLikelihood()

{

	// update sequence from SVE (single variable exchange)
	sampleSVE();

	// TODO: eventually, turn this into a parameter (or constant?)
	double beta = 0.0;

	Topology* t = getActiveTopology();
	Node* r	= t->getRoot();

	double lnRootStationaryProb = 0.0;
	double lnSampleStationaryProb = 0.0;

	std::vector<double> f = getActiveFreqs()->getFreqs();
	std::vector<Step> rootSeq = r->getActivePath()->getSeqState();

	// calculate root stationary probability from sampled sequences
	//	pr(s^(root)|theta') / pr(s^root)|theta) ~=
	//			e^(-2(p'-p)*E(s^(root)))
	//			* (*=[m=1,N]:(pi'_i__m / pi_i__m))
	//			* (+=[h=1,N]:(e^(-2(p -p*)*E(nu^(h))) * (*=[n=1,N]:(pi _(nu^(h)_n) / pi*_(nu^(h)_n))
	//			/ (+=[h=1,N]:(e^(-2(p'-p*)*E(nu^(h))) * (*=[n=1,N]:(pi'_(nu^(h)_n) / pi*_(nu^(h)_n))

	for (int i = 0; i < numCodonSites; i++) {
		lnRootStationaryProb += log(f[rootSeq[i].getState()]);
		lnSampleStationaryProb += log(f[sequenceSVE[i].getState()]);
	}

	lnRootStationaryProb += exp(-2 * beta * calcPseudoEnergy(rootSeq));
	lnSampleStationaryProb += exp(-2 * beta * calcPseudoEnergy(sequenceSVE));

	return lnRootStationaryProb - lnSampleStationaryProb;
}

double ModelSeq::lnBranchSiteLikelihood(Node* p, int site)
{

	double unitBrlen = 0.0;
	double unitTimeMarginOfError = 0.0000000000000005;

	int currState = 0, prevState = 0;
	double unitCurrTime = -1.0, unitPrevTime = 0.0, unitFinalTime = 0.0;
	double lnStepProbSum = 0.0;
	double rateAwaySum = 0.0;
	double lnStepProb = 0.0, stepProb = 0.0;
	double x = 0.0;
	std::vector<Step> srcPath;
	std::vector<Step> destPath;
	std::vector<double> u = getActiveSubstRate()->getSubstRate();

	if (p->getAnc() != NULL)
	{
		// FIND		pr(j,rho|i,theta) * pr(i|theta)
		// 		where
		// 			pr(j,rho|i,theta) 	= [*=:events, R_i(z-1),i(z) * e^(-R_i(z-1),# * (t(z) - t(z-1))] * e^(-R_i(q),# * (t(q+1) - t(q)))

		unitBrlen = u[p->getDnPassIndex()];

		/*
		std::cout << "\n";
		std::cout << "path describing n" << p->getAnc()->getIndex() << "-> n" << p->getIndex() << "   (u[" << p->getDnPassIndex() << "]: " << unitBrlen << ")\n";
		*/

		// intialize src event and dest events
		std::vector<Step> eventVector = p->getActivePath()->getSitePathHistory(site);
		reverse(eventVector.begin(), eventVector.end());
		int prevSiteState = p->getAnc()->getActivePath()->getSitePathHistory(site)[0].getState();

		/*
		for(std::vector<Step>::iterator it = eventVector.begin(); it != eventVector.end(); it++)
		{
			std::cout << "p           state: " << (*it).getState() << "\n";
		}
		std::cout <<     "p->getAnc() state: " << prevSiteState << "\n";
		*/

		// calculate likelihood for events
		for (std::vector<Step>::iterator it = eventVector.begin(); it != eventVector.end(); it++)
		{
			//std::cout << "site: " << site << ", state: " << (*it).getState() << "\n";

			unitCurrTime = (*it).getTime() * unitBrlen;
			if (1.0 - unitTimeMarginOfError < unitCurrTime && 1.0 + unitTimeMarginOfError > unitCurrTime)
				unitCurrTime = 1.0;

			if (unitCurrTime < 1.0 - unitTimeMarginOfError)
				unitFinalTime = unitCurrTime;

			// if the scaled event time == 1.0, it is a destination event
			if (unitCurrTime == 1.0) {
				prevSiteState = (*it).getState();
				unitPrevTime = unitCurrTime;
			}

			// otherwise, calculate the probability of the event occurring
			else if (unitPrevTime != unitCurrTime && unitCurrTime != 1.0)
			{
				rateAwaySum = 0.0;
				lnStepProb = 0.0;

				prevState = prevSiteState;
				currState = (*it).getState();

				if (q[prevState][currState] == 0.0) {
					std::cout << "WARNING: q[" << prevState << "][" << currState << "] == 0.0\n";
				}

				// do not consider virtual events
				if (prevState != currState) {
					// get the sum of rates for all sequences away from prevSeq
					// (numCodonSites) * (numSenseCodons - 1) sequences differ by one nucleotide
					for (int i = 0; i < numCodonSites; i++) {
						for (int j = 0; j < numSenseCodons; j++) {
							if (prevSiteState != j) rateAwaySum += q[prevSiteState][j];
						}
					}
					rateAwaySum *= unitBrlen;

					stepProb = q[prevSiteState][currState] * unitBrlen;
					lnStepProb = log(stepProb);

					x = -rateAwaySum * (unitCurrTime - unitPrevTime);
					lnStepProb += x;

					/*
					std::cout << prevState << "->" << currState << ":s" << (*it).getSite();
					std::cout << ", lnstepProb: "			<< lnStepProb;
					std::cout << ", R_ij: "				<< q[prevState][currState] * unitBrlen;
					std::cout << ", rateAwaySum: "		<< rateAwaySum;
					std::cout << ", dt: "				<< unitCurrTime - unitPrevTime;
					std::cout << ", -R_i* * dt: "		<< -rateAwaySum * unitBrlen * (unitCurrTime - unitPrevTime);
					std::cout << ", e^(-R_i* * dt): "	<<  exp(-rateAwaySum * unitBrlen * (unitCurrTime - unitPrevTime));
					std::cout << "\n";
					*/

					lnStepProbSum += lnStepProb;

					prevSiteState = currState;
					unitCurrTime = unitPrevTime;
				}
			}

			else if (unitPrevTime == unitCurrTime && unitCurrTime != 1.0)
				std::cout << "ERROR: duplicate event times, prevTime: " << unitPrevTime << ", currTime: " << unitCurrTime << "\n";

		}

		// calculate final term
		rateAwaySum = 0.0;
		for (int j = 0; j < numSenseCodons; j++) {
			if (prevSiteState != j) rateAwaySum += q[prevSiteState][j];
		}

		x = -rateAwaySum * unitBrlen * (1.0 - unitFinalTime);
		lnStepProbSum += x;
		// std::cout << "lnStepProbSum: " << lnStepProbSum << "\n";
	}


	return lnStepProbSum;
}


void ModelSeq::initSVE() {

	//	sample possible alternative root sequences

	double stationaryFreqSum = 0.0;
	double stationaryFreqSample = 0.0;
	double sampleRV = 0.0;
	std::vector<double> f = getActiveFreqs()->getFreqs();
	Step tempStep;

	// initialize sequenceSVE size to hold data
	sequenceSVE.resize(numCodonSites);

	// assign a new state to sampleSeq
	for (int j = 0; j < numSenseCodons; j++)
		stationaryFreqSum += f[j];

	for (int i = 0; i < numCodonSites; i++)
	{
		stationaryFreqSample = 0.0;
		sampleRV = (double)rand() / (double)RAND_MAX * stationaryFreqSum;

		for (int j = 0; j < numSenseCodons; j++)
		{
			stationaryFreqSample += f[j];
			if (stationaryFreqSample > sampleRV)
			{
				sequenceSVE[i] = Step(0.0, j, i);
				break;
			}
		}
	}
}


void ModelSeq::sampleSVE() {

	// update sequenceSVE at a number of random sites equal to SVESamplesPerCycle (default = 5)

	double stationaryFreqSum = 0.0;
	double stationaryFreqSample = 0.0;
	double sampleRV = 0.0;
	int randSite = 0;
	std::vector<double> f = getActiveFreqs()->getFreqs();

	// assign a new state to sampleSeq
	for (int j = 0; j < numSenseCodons; j++)
		stationaryFreqSum += f[j];

	for (int i = 0; i < SVESamplesPerCycle; i++)
	{
		sampleRV = (double)rand() / (double)RAND_MAX * stationaryFreqSum;
		randSite = numCodonSites * ((double)rand() / (double)RAND_MAX);
		for (int j = 0; j < numSenseCodons; j++)
		{
			stationaryFreqSample += f[j];
			if (stationaryFreqSample > sampleRV)
			{
				sequenceSVE[randSite] = Step(0.0, j, randSite);
				break;
			}
		}
	}
}

void ModelSeq::printSVE() {
	for (std::vector<Step>::iterator it = sequenceSVE.begin(); it != sequenceSVE.end(); it++)
	{
		std::cout << std::setw(2) << (*it).getState() << " ";
	}
	std::cout << "\n";
}

Freqs* ModelSeq::getActiveFreqs(void) {

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

/*
Contexts* ModelSeq::getActiveContexts(void) {

	for (std::vector<Parm*>::iterator p=parameters.begin(); p != parameters.end(); p++)
		{
		ContextDependence* derivedPtr = dynamic_cast<ContextDependence*> (*p);
		if (derivedPtr != 0)
			{
			return derivedPtr->getActiveContexts();
			break;
			}
		}
	return NULL;
}
*/

Dnds* ModelSeq::getActiveDnds(void) {

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

TiTv* ModelSeq::getActiveTiTv(void) {

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

Topology* ModelSeq::getActiveTopology(void) {

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

Tree* ModelSeq::getActiveTree(void) {

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

Pairwise* ModelSeq::getActivePairwise(void) {
	for (std::vector<Parm*>::iterator p=parameters.begin(); p != parameters.end(); p++)
	{
		SeqStructP* derivedPtr = dynamic_cast<SeqStructP*> (*p);
		if (derivedPtr != 0)
		{
			return derivedPtr->getActivePairwise();
			break;
		}
	}
	return NULL;
}

SolvAcc* ModelSeq::getActiveSolvAcc(void) {


	for (std::vector<Parm*>::iterator p=parameters.begin(); p != parameters.end(); p++)
	{
		SeqStructS* derivedPtr = dynamic_cast<SeqStructS*> (*p);
		if (derivedPtr != 0)
		{
			return derivedPtr->getActiveSolvAcc();
			break;
		}
	}
	return NULL;
}

SubstRate* ModelSeq::getActiveSubstRate(void) {

	for (std::vector<Parm*>::iterator p=parameters.begin(); p != parameters.end(); p++)
	{
		ScalerU* derivedPtr = dynamic_cast<ScalerU*> (*p);
		if (derivedPtr != 0)
		{
			return derivedPtr->getActiveSubstRate();
			break;
		}
	}
	return NULL;
}


Parm* ModelSeq::pickParmAtRandom(void) {

	double u = ranPtr->uniformRv();
	double sum = 0.0;
	for (int i=0; i < (int)proposalProbs.size(); i++) {
		sum += proposalProbs[i];
		if (u < sum)
			return parameters[i];
	}
	return NULL;
}

void ModelSeq::updateAllFlags(void) {

	Topology* t = getActiveTopology();
	t->updateAllCls(true);
	t->updateAllTis(true);
	t->flipAllActiveCls();
	t->flipAllActiveTis();
	////tiProbs->setTransitionProbs();
	t->updateAllTis(true);
}

void ModelSeq::updateQ(Parm* p) {

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

	{
		SeqStructS* derivedPtr = dynamic_cast<SeqStructS*> (p);
		if (derivedPtr != 0)
		{
			shouldUpdate = true;
			goto funcEnd;
		}
	}

	{
		SeqStructP* derivedPtr = dynamic_cast<SeqStructP*> (p);
		if (derivedPtr != 0)
		{
			shouldUpdate = true;
			goto funcEnd;
		}
	}

	{
		ScalerU* derivedPtr = dynamic_cast<ScalerU*> (p);
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
			////tiMatrix->updateQ(q);
		}
}

void ModelSeq::restoreQ(Parm* p) {

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

	{
		SeqStructS* derivedPtr = dynamic_cast<SeqStructS*> (p);

		if (derivedPtr != 0)
		{
			shouldUpdate = true;
			goto funcEnd;
		}
	}

	{
		SeqStructP* derivedPtr = dynamic_cast<SeqStructP*> (p);
		if (derivedPtr != 0)
		{
			shouldUpdate = true;
			goto funcEnd;
		}
	}

	{
		ScalerU* derivedPtr = dynamic_cast<ScalerU*> (p);
		if (derivedPtr != 0)
		{
			shouldUpdate = true;
			goto funcEnd;
		}
	}

	funcEnd:
		if (shouldUpdate == true)
		{
			////tiMatrix->restoreQ();
		}
}

void ModelSeq::keep(Parm* p) {

	p->keep();

	Tree* derivedPtr = dynamic_cast<Tree*> (p);
	if (derivedPtr == 0)
	{
		getActiveTree()->keep();
	}
}

void ModelSeq::restore(Parm* p) {

	p->restore();

	Tree* derivedPtr = dynamic_cast<Tree*> (p);
	if (derivedPtr == 0)
	{
		getActiveTree()->restore();
	}

}

void ModelSeq::printQ(void) {

	std::cout << "ModelSeq: Q MATRIX\n";
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
