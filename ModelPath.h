/*
 * Model.h
 *
 *  Created on: Jan 24, 2010
 *      Author: mikee
 */

#ifndef MODEL_PATH_H_
#define MODEL_PATH_H_

#include "Alignment.h"
#include "Code.h"
#include "CondLikes.h"
#include "FileMgr.h"
#include "MbRandom.h"
#include "MbTransitionMatrix.h"
#include "Parm_context.h"
#include "Parm_freqs.h"
#include "Parm_kappa.h"
#include "Parm_omega.h"
#include "Parm_tree.h"
#include "Parm_scaling.h"
#include "Path.h"
#include "Settings.h"
#include "TiProbs.h"
#include "Util.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

class Alignment;
class CondLikes;
class Dnds;
class FileMgr;
class Freqs;
class MbRandom;
class MbTransitionMatrix;
class Node;
class Parm;
class Settings;
class Step;
class SubstRate;
class TiProbs;
class TiTv;
class Topology;
class Tree;

class ModelPath {

public:
						ModelPath(Alignment *ap, MbRandom *rp, Settings* sp, Tree* tp);
						~ModelPath(void);
	Dnds*				getActiveDnds(void);
	Freqs*				getActiveFreqs(void);
	SubstRate*			getActiveSubstRate(void);
	TiTv*				getActiveTiTv(void);
	Topology*			getActiveTopology(void);
	Tree*				getActiveTree(void);
	Parm*				pickParmAtRandom(void);
	Parm*				getParameter(int i)				{ return parameters[i]; }
	int					getNumParameters(void)			{ return parameters.size(); }
	void				fillInRateMatrix(void);			// fills in both Q & R
	void				printQ(void);
	void				printR(void);
	void				printAcceptanceInfo(void);
	void				updateAllFlags(void);
	void				updateQ(Parm* p);
	void				restoreQ(Parm* p);
	void				keep(Parm* p);
	void				restore(Parm* p);
	void				calcProbs(double **dnds, double **posPr);
	double				calcDwellTime(int i, int j, double v, int k);
	double				calcExpectedOmega(int i, int j, double v, int k);
	void				createPathMappings(void);
	void				changePathMappings(void);
	void				keepPathMappings(void);
	void				restorePathMappings(void);
	void				sampleNodeSiteState(Node* p, int codonSite);
	void				sampleNodeSitePath(Node* p, int codonSite);

private:
	void				initializeCondLikes(int codonSite);
	void				initializeRateMatrix(int nc);
	void				initializeEmpiricalCodonMatrix(std::string ecmFilePath, int nc);
	std::vector<double> initializeEmpiricalCodonFreqs(std::string ecmFilePath, int nc);
	void				copyECMtoR();
	bool				isTransversion(int a, int b);
	void				fillInStationaryProbs(void);
	Alignment*			alignmentPtr;
	MbRandom*			ranPtr;
	Tree*				treePtr;
	Settings*			settingsPtr;
	MbTransitionMatrix*	tiMatrix;
	std::vector<std::vector<int> >	qTemplate;
	MbMatrix<double>	q;							// used to initialize r
	MbMatrix<double>	r;							// used for path sampling
	std::vector<Parm*>	parameters;
	std::vector<double>	proposalProbs;
	double*				lnScalers;
	double*				stationaryProbs;
	double				uniformizationRate;
	int					numStates;
	int					numSenseCodons;
	int					numCodonSites;
	CondLikes*			condLikes;
	TiProbs*			tiProbs;
};

#endif /* MODEL_PATH_H_ */
