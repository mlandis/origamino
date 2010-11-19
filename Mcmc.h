#ifndef MCMC_H
#define MCMC_H

#include "MbRandom.h"
#include "Mcmc.h"
#include "ModelPath.h"
#include "ModelSeq.h"
#include "Parm.h"
#include "Parm_tree.h"
#include "Settings.h"
//#include "TreeSum.h"

#include <math.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <string>

class Alignment;
class FileMgr;
class MbRandom;
class ModelPath;
class ModelSeq;
class Parm;
class Settings;
//class TreeSum;

class Mcmc {

	public:
                            Mcmc(Alignment* ap, Settings* sp, MbRandom* rp, ModelPath* mpp, ModelSeq* msp, FileMgr* om);
							~Mcmc(void);

	private:
	                 void   printChainState(int n, double lnL);
	                 void   runChain(void);
				   double   safeExp(double lnX);
			  std::string   trueFalse(bool tf);
			         void   openFiles(std::string fn);
				   ModelPath*   modelPathPtr;
				   ModelSeq*	modelSeqPtr;
				MbRandom*   ranPtr;
				Settings*   settingsPtr;
			   Alignment*   alignmentPtr;
					  int   numCycles;
					  int   printFrequency;
					  int   sampleFrequency;
					  int   summarizeFrequency;
					  int	SVESamplesPerCycle;
					  int	pathSampleFrequency;
					  int	siteChangesPerCycle;
			std::ofstream   treeFileStrm;
			std::ofstream   parmFileStrm;
	//		     TreeSum*   treeSummary;
};


#endif
