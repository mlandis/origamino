/*
 * Parm_pairwise.h
 *
 *  Created on: Apr 22, 2010
 *      Author: mikee
 *
 *      pairwise interactions
 */

#ifndef PARM_PAIRWISE_H_
#define PARM_PAIRWISE_H_

#include "MbRandom.h"
#include "Parm.h"
#include "Parm_tree.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

class Alignment;
class MbRandom;
class Tree;

class Pairwise {

public:
				Pairwise(MbRandom* rp, Tree* tp);
				Pairwise(Pairwise &d);
	Pairwise	&operator=(const Pairwise &d);
	double		change(void);
	double		getPairwise(void)				{ return pairwiseVal; }
	double		lnProbability(void);
	void		print(void);

private:
	MbRandom*	ranPtr;
	Tree*		treePtr;
	void		clone(const Pairwise &d);
	double		pairwiseVal;

};

class SeqStructP: public Parm
{

public:
							SeqStructP(MbRandom* rp, Tree* tp, std::string pn);
							~SeqStructP(void);
	Pairwise*				getActivePairwise(void)									{return pairwise[activeState]; }
	double					lnPriorRatio(void);
	double					change(void);
	void					print(void);
	void					keep(void);
	void					restore(void);
	std::string				getParameterStr(void);
	std::string				getParameterHeader(void);

private:
	Pairwise				*pairwise[2];

};

#endif /* PARM_PAIRWISE_H_ */
