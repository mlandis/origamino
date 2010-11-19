/*
 * Parm_context.h
 *
 *  Created on: Jan 26, 2010
 *      Author: mikee
 */

#ifndef PARM_CONTEXT_H_
#define PARM_CONTEXT_H_

#include "Parm.h"
#include "Parm_tree.h"
#include "MbRandom.h"
#include "ModelPath.h"

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

class Alignment;
class MbRandom;
class Tree;

class Contexts {

public:
								Contexts(MbRandom* rp, Tree* tp);
								Contexts(Contexts &c);
	Contexts					&operator=(const Contexts &ctxts);
	double						change(void);
	std::vector<double>&		getContexts(void)			{ return contextProbs; }
	double						getContextVal(int x)		{ return contextProbs[x]; }
	double						lnProbability(void);
	void						setContextNames();
	void						print(void);

private:
	MbRandom*					ranPtr;
	Tree*						treePtr;
	void						clone(const Contexts &c);
	void						normalizeContextProbs(std::vector<double>& v, double threshold);
	int							numContexts;
	std::vector<double>			a;
	std::vector<double>			contextProbs;
	std::vector<std::string>	contextNames;
};

class ContextDependence : public Parm {

public:
								ContextDependence(MbRandom* rp, Tree* tp, std::string pn);
								~ContextDependence(void);
	Contexts*					getActiveContexts(void)		{ return contexts[activeState]; }
	double						lnPriorRatio(void);
	double						change(void);
	void						print(void);
	void						keep(void);
	void						restore(void);
	std::string					getParameterStr(void);
	std::string					getParameterHeader(void);

private:
	Contexts					*contexts[2];

};

#endif /* PARM_CONTEXT_H_ */
