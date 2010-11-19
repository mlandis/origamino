/*
 * ModelSeq.h
 *
 *  Created on: Jan 24, 2010
 *      Author: mikee
 */

#ifndef MODEL_SEQ_H_
#define MODEL_SEQ_H_

#include "Alignment.h"
#include "Code.h"
#include "CondLikes.h"
#include "MbRandom.h"
#include "MbTransitionMatrix.h"
#include "Parm_context.h"
#include "Parm_freqs.h"
#include "Parm_kappa.h"
#include "Parm_omega.h"
#include "Parm_tree.h"
#include "Parm_scaling.h"
#include "Parm_solv.h"
#include "Parm_pairwise.h"
#include "Path.h"
#include "Settings.h"
#include "TiProbs.h"
#include "Util.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <vector>

class Alignment;
class CondLikes;
class Contexts;
class Dnds;
class Freqs;
class MbRandom;
class MbTransitionMatrix;
class Node;
class Pairwise;
class Parm;
class Settings;
class SolvAcc;
class Step;
class SubstRate;
class TiProbs;
class TiTv;
class Topology;
class Tree;

class ModelSeq {

public:
						ModelSeq(Alignment *ap, MbRandom *rp, Settings* sp, Tree* tp);
						~ModelSeq(void);
	Contexts*			getActiveContexts(void);
	Dnds*				getActiveDnds(void);
	Freqs*				getActiveFreqs(void);
	Pairwise*			getActivePairwise(void);
	SolvAcc*			getActiveSolvAcc(void);
	SubstRate*			getActiveSubstRate(void);
	TiTv*				getActiveTiTv(void);
	Topology*			getActiveTopology(void);
	Tree*				getActiveTree(void);
	Parm*				pickParmAtRandom(void);
	Parm*				getParameter(int i)				{ return parameters[i]; }
	int					getNumParameters(void)			{ return parameters.size(); }
	double				lnPathLikelihood(void);
	double				lnRootLikelihood(void);
	double				lnBranchSiteLikelihood(Node* p, int site);
	void				initSVE(void);
	void				sampleSVE(void);
	void				printSVE(void);
	void				fillInRateMatrix(void);
	void				printAcceptanceInfo(void);
	void				updateAllFlags(void);
	void				updateQ(Parm* p);
	void				restoreQ(Parm* p);
	void				printQ(void);
	void				keep(Parm* p);
	void				restore(Parm* p);
	void				calcProbs(double **dnds, double **posPr);
	//double			calcDwellTime(int i, int j, double v, int k);

private:
	void							initializeCondLikes(int codonSite);
	void							initializeRateMatrix(int nc);
	void							initializeContactMap(std::string pdbFilePath, int nc, double cd);
	double							calcPseudoEnergy(std::vector<Step> seq);
	bool							isTransversion(int a, int b);
	void							fillInStationaryProbs(void);
	Alignment*						alignmentPtr;
	Settings*						settingsPtr;
	Tree*							treePtr;
	MbRandom*						ranPtr;
	MbTransitionMatrix*				tiMatrix;
	MbMatrix<double>				q;
	std::vector<std::vector<int> >	qTemplate;
	std::vector<double>				qRateAway;
	std::vector<Parm*>				parameters;
	std::vector<double>				proposalProbs;
	std::vector<Step>				sequenceSVE;
	std::vector<std::vector <int> >	contactMap;
	double*							lnScalers;
	double*							stationaryProbs;
	int								numStates;
	int								numSenseCodons;
	int								numCodonSites;
	int								numSampleSeqs;
	double							brlenLambda;
	double							contactDistance;
	int								SVESamplesPerCycle;
	CondLikes*						condLikes;
	TiProbs*						tiProbs;
};

#endif /* MODEL_SEQ_H_ */
