/*
 * Parm_solv.h
 *
 *  Created on: Apr 22, 2010
 *      Author: mikee
 *
 *      solvent accessibility
 */


#ifndef PARM_SOLV_H_
#define PARM_SOLV_H_

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

class SolvAcc {

public:
				SolvAcc(MbRandom* rp, Tree *tp);
				SolvAcc(SolvAcc &d);
	SolvAcc		&operator=(const SolvAcc &d);
	double		change(void);
	double		getSolvAcc(void)				{ return solvaccVal; }
	double		lnProbability(void);
	void		print(void);

private:
	MbRandom*	ranPtr;
	Tree*		treePtr;
	void		clone(const SolvAcc &d);
	double		solvaccVal;

};

class SeqStructS: public Parm
{

public:
							SeqStructS(MbRandom* rp, Tree* tp, std::string pn);
							~SeqStructS(void);
	SolvAcc*				getActiveSolvAcc(void)									{return solvacc[activeState]; }
	double					lnPriorRatio(void);
	double					change(void);
	void					print(void);
	void					keep(void);
	void					restore(void);
	std::string				getParameterStr(void);
	std::string				getParameterHeader(void);

private:
	SolvAcc					*solvacc[2];

};

#endif /* PARM_SOLV_H_ */
