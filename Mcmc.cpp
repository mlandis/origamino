#include "Mcmc.h"

Mcmc::Mcmc(Alignment* ap, Settings* sp, MbRandom* rp, ModelPath* mpp, ModelSeq* msp, FileMgr* om) {

	// remember important objects
	settingsPtr  = sp;
	ranPtr       = rp;
	modelPathPtr = mpp;
	modelSeqPtr  = msp;
	alignmentPtr = ap;
	
	// initialize parameters of the chain
	numCycles				= sp->getChainLength();
	printFrequency			= sp->getPrintFrequency();
	sampleFrequency			= sp->getSampleFrequency();
	summarizeFrequency		= sp->getSummarizeFrequency();
	SVESamplesPerCycle		= sp->getSVESamplesPerCycle();
	siteChangesPerCycle		= sp->getSiteChangesPerCycle();
	
	// open the output files
	openFiles( settingsPtr->getOutputFileName() );
	
	// create object that keeps track of probabilities on the tree
	//treeSummary = new TreeSum(alignmentPtr, modelPtr, settingsPtr->getOutputFileName(), sampleFrequency);
	
	// run the Markov chain
	runChain();
}

Mcmc::~Mcmc(void) {

	// close the output files
	treeFileStrm.close();
	parmFileStrm.close();
}

void Mcmc::runChain(void) {

	std::cout << "Mcmc::runChain() called.\n";

	// get initial likelihood
	modelPathPtr->updateAllFlags();
	modelSeqPtr->updateAllFlags();

	// create mappings
	modelPathPtr->createPathMappings();

	double oldLnL = modelSeqPtr->lnPathLikelihood() + modelSeqPtr->lnRootLikelihood();

	std::cout << "oldLnl = " << oldLnL << "\n";


	// cycle events:
	//
	//	BranchHistory - Resample branch history for one randomly selected branch.
	//	NodeState - Resample node state for randomly selected node. Then, resample branch histories for all adjacent branches. (i.e. 1-3x BranchHistory)
	//	GibbsResample - Resample the M sequences for the Gibbs Sampler. Update theta* as needed.
	//	ParametersPath - Change one parameter for the ModelPath -- these should get tuned, right? shouldn't they relate to ParametersSeq?
	//	ParametersSeq - Change one parameter for the ModelSeq

	for (int n=1; n<=numCycles; n++)
	{
		// PARAMETER UPDATES theta'|theta

		// pick a parameter to change
		Parm* parm = modelSeqPtr->pickParmAtRandom();
		
		// TODO: pick site residues to change -- how often? every 500 cycles??
		// TODO: update relevant paths

		// propose a new state for the parameter
		double lnProposalRatio = parm->change();
		//	std::cout << "lnProposalRatio: " << lnProposalRatio << "\n";

		// calculate the prior ratio 
		double lnPriorRatio = parm->lnPriorRatio();
		//	std::cout << "lnPriorRatio: " << lnPriorRatio << "\n";

		// calculate the likelihood
		modelSeqPtr->updateQ(parm);
		double lnL = modelSeqPtr->lnPathLikelihood() + modelSeqPtr->lnRootLikelihood();
		//	std::cout << "lnPathLikelihood + lnRootLikelihood: " << lnL << "\n";

		double lnLikelihoodRatio = lnL - oldLnL;
		//	std::cout << "lnLikelihoodRatio: " << lnLikelihoodRatio << "\n";

		double r = safeExp( lnLikelihoodRatio + lnPriorRatio + lnProposalRatio );
		
		// accept or reject the proposed state as the next state of the chain
		bool acceptState = false;
		if (ranPtr->uniformRv() < r)
			acceptState = true;
			
		// print information to the screen
		if ( n % printFrequency == 0 )
		{
			std::cout << std::setw(5) << n << " -- ";
			std::cout << std::fixed << std::setprecision(8) << std::setw(16) << oldLnL << " -> " << std::setprecision(8) << std::setw(16) << lnL;
			if (acceptState == true)
				std::cout << " -- Accepted          change to " << parm->getName() << std::endl;
			else
				std::cout << " --          Rejected change to " << parm->getName() << std::endl;
		}

		// update the state of the chain
		if (acceptState == true)
		{
			parm->incrementNumAccepted();
			oldLnL = lnL;
			modelSeqPtr->keep(parm);
		}
		else 
		{
			modelSeqPtr->restore(parm);
			modelSeqPtr->restoreQ(parm);
		}



		// PATH UPDATES rho'|rho

		Topology* t = modelSeqPtr->getActiveTopology();
		double lnOldPathLikelihood = 0.0;
		double lnNewPathLikelihood = 0.0;

		for (int i = 0; i < siteChangesPerCycle; i++)
		{
			double u = (double)rand() / (double)RAND_MAX;
			int randSite = (int)((double)rand() / (double)RAND_MAX * alignmentPtr->getNumCodonSites());
			int randNode = (int)((double)rand() / (double)RAND_MAX * t->getNumNodes());
			//randNode = 8;

			// branch site resampling occurs at pr(resampleRate)
			// node site & adjacent branch site resampling occurs at pr (1.0 - resampleRate)
			double resampleRate = 0.75;
			bool resampleNode = false;

			Node* p = t->getDownPassNode(randNode);

			// resample path along branch for a site


			if (u > resampleRate) resampleNode = true;

			/*
			std::cout << "randSite: " << randSite << ", randNode: " << randNode << ", p->getIndex(): " << p->getIndex();
			if (resampleNode)	std::cout << " --- resample NODE\n";
			else				std::cout << " --- resample BRANCH\n";
			*/

			if (resampleNode == false)
			{

				while (p->getAnc() == NULL)
				{
					randNode = (int)((double)rand() / (double)RAND_MAX * t->getNumNodes());
					p = t->getDownPassNode(randNode);
				}

				if (p->getAnc() != NULL) {
					lnOldPathLikelihood += modelSeqPtr->lnBranchSiteLikelihood(p, randSite);
					modelPathPtr->sampleNodeSitePath(p, randSite);
					lnNewPathLikelihood += modelSeqPtr->lnBranchSiteLikelihood(p, randSite);
				}
			}

			// resample internal node and any adjacent branches for a site
			else if (resampleNode == true)
			{
				// ensure node is an internal node


				while (p->getAnc() == NULL || (p->getLft() == NULL && p->getRht() == NULL))
				{
					randNode = (int)((double)rand() / (double)RAND_MAX * t->getNumNodes());
					p = t->getDownPassNode(randNode);
				}

				if (p->getAnc() != NULL)
				{
					lnOldPathLikelihood += modelSeqPtr->lnBranchSiteLikelihood(p, randSite);
				}

				if (p->getLft() != NULL)
				{
					lnOldPathLikelihood += modelSeqPtr->lnBranchSiteLikelihood(p->getLft(), randSite);
				}

				if (p->getRht() != NULL)
				{
					lnOldPathLikelihood += modelSeqPtr->lnBranchSiteLikelihood(p->getRht(), randSite);
				}

				modelPathPtr->sampleNodeSiteState(p, randSite);

				if (p->getAnc() != NULL)
				{
					modelPathPtr->sampleNodeSitePath(p, randSite);
					lnNewPathLikelihood += modelSeqPtr->lnBranchSiteLikelihood(p, randSite);
				}

				if (p->getLft() != NULL)
				{
					modelPathPtr->sampleNodeSitePath(p->getLft(), randSite);
					lnNewPathLikelihood += modelSeqPtr->lnBranchSiteLikelihood(p->getLft(), randSite);
				}

				if (p->getRht() != NULL)
				{
					modelPathPtr->sampleNodeSitePath(p->getRht(), randSite);
					lnNewPathLikelihood += modelSeqPtr->lnBranchSiteLikelihood(p->getRht(), randSite);
				}
			}
		}

		double lnPathLikelihoodRatio = lnNewPathLikelihood - lnOldPathLikelihood;

		// TODO configure proper Hastings ratio
		r = safeExp( lnPathLikelihoodRatio + lnPriorRatio + lnProposalRatio );

		// accept or reject the proposed state as the next state of the chain
		acceptState = false;
		if (ranPtr->uniformRv() < r)
			acceptState = true;

		// update the state of the chain
		if (acceptState == true)
		{
			for (int n = 0; n < t->getNumNodes(); n++)
			{
				Node* p = t->getDownPassNode(n);
				p->getHistory()->incrementNumAccepted();
				p->getHistory()->keep();
				if (p->getRht() != NULL) p->getRht()->getHistory()->keep();
				if (p->getLft() != NULL) p->getLft()->getHistory()->keep();
			}
		}
		else
		{
			for (int n = 0; n < t->getNumNodes(); n++)
			{
				Node* p = t->getDownPassNode(n);
				p->getHistory()->restore();
				if (p->getRht() != NULL) p->getRht()->getHistory()->restore();
				if (p->getLft() != NULL) p->getLft()->getHistory()->restore();
			}
		}

		if ( n % printFrequency == 0 )
		{
			std::cout << std::setw(5) << n << " -- ";
			std::cout << std::fixed << std::setprecision(8) << std::setw(16) << lnOldPathLikelihood << " -> " << std::setprecision(8) << std::setw(16) << lnNewPathLikelihood;
			if (acceptState == true)
			{
				std::cout << " -- Accepted          change to path mapping\n";
			}
			else
			{
				std::cout << " --          Rejected change to path mapping\n";
			}
		}



		// print out current chain state

		//if ( n == 1 || n % sampleFrequency == 0 || n == numCycles )
		//	printChainState(n, oldLnL);
			
		// print summary of chain
		// if ( n % summarizeFrequency == 0 ) treeSummary->printSummary();

		//std::cout << "\n";
	}

	//treeSummary->print();
	//modelPtr->printAcceptanceInfo();
}

void Mcmc::openFiles(std::string fn) {

	std::string tf = fn + ".t";
	std::string pf = fn + ".p";
	
	treeFileStrm.open( tf.c_str(), std::ios::out );
	if ( !treeFileStrm )
		{
		std::cerr << "ERROR: Problem opening tree output file" << std::endl;
		exit(1);
		}

	parmFileStrm.open( pf.c_str(), std::ios::out );
	if ( !parmFileStrm )
		{
		std::cerr << "ERROR: Problem opening parameter output file" << std::endl;
		exit(1);
		}
}

//TODO: warning; modelPathPtr used simply for compilation. re-evaluate later.
void Mcmc::printChainState(int n, double lnL) {

	if (n == 1)
		{
		std::string pHeaderStr = "";
		std::string tHeaderStr = "";
		for (int i=0; i<modelPathPtr->getNumParameters(); i++)
			{
			Parm* p = modelPathPtr->getParameter(i);
			Tree* derivedPtr = dynamic_cast<Tree*> (p);
			if (derivedPtr != 0)
				tHeaderStr += p->getParameterHeader();
			else 
				pHeaderStr += p->getParameterHeader();
			}

		treeFileStrm << tHeaderStr;
		parmFileStrm << "Cycle\tlnL\t" << pHeaderStr << std::endl;
		}
		
	std::string pStr = "";
	std::string tStr = "";
	for (int i=0; i<modelSeqPtr->getNumParameters(); i++)
		{
		Parm* p = modelPathPtr->getParameter(i);
		Tree* derivedPtr = dynamic_cast<Tree*> (p);
		if (derivedPtr != 0)
			tStr += p->getParameterStr();
		else 
			pStr += p->getParameterStr();
		}
	treeFileStrm << "   tree_" << n << " = " << tStr << ";" << std::endl;
	parmFileStrm << n << '\t' << std::fixed << std::setprecision(2) << lnL << '\t' << pStr << std::endl;
	
	if (n == numCycles)
		{
		std::cout << "Mcmc::runChain() complete!\n";
		treeFileStrm << "end;" << std::endl;
		}
		
	//if (n != 1) treeSummary->addState(n);
}

double Mcmc::safeExp(double lnX) {

	if (lnX < -300.0)
		return 0.0;
	else if (lnX > 0.0)
		return 1.0;
	else
		return exp(lnX);
}

std::string Mcmc::trueFalse(bool tf) {

	if (tf == true)
		return "True";
	return "False";
}
