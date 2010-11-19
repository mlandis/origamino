/*
 * Parm_scaling.h
 *
 *  Created on: May 4, 2010
 *      Author: mikee
 */

#ifndef PARM_SCALING_H_
#define PARM_SCALING_H_

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

class SubstRate {

public:
							SubstRate(MbRandom* rp, Tree* tp, double ep);
							SubstRate(SubstRate &d);
	SubstRate				&operator=(const SubstRate &d);
	double					change(void);
	std::vector<double>		getSubstRate(void)					{ return u; }
	double					getSubstRateVal(int x)				{ return u[x]; }
	double					lnProbability(void);
	void					print(void);

private:
	MbRandom*				ranPtr;
	Tree*					treePtr;
	double					brlenLambda;
	void					clone(const SubstRate &d);
	std::vector<double> 	u;									// scaling vector stored [0,n] in order of getDownPassNode()

};


class ScalerU : public Parm {

public:
							ScalerU(MbRandom *rp, Tree* tp, double ep, std::string pn);
							~ScalerU(void);
	SubstRate*				getActiveSubstRate(void)			{return substRate[activeState]; }
	double					lnPriorRatio(void);
	double					change(void);
	void					print(void);
	void					keep(void);
	void					restore(void);
	std::string				getParameterStr(void);
	std::string				getParameterHeader(void);

private:
	SubstRate		*substRate[2];

};

#endif /* PARM_SCALING_H_ */
