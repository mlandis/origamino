/*
 * Parm_kappa.h
 *
 *  Created on: Jan 26, 2010
 *      Author: mikee
 */

#ifndef PARM_KAPPA_H_
#define PARM_KAPPA_H_

#include "MbRandom.h"
#include "Parm.h"
#include "Parm_tree.h"

#include <iomanip>
#include <iostream>
#include <string>

class MbRandom;
class Tree;

class TiTv {

public:
				TiTv(MbRandom* rp, Tree* tp);
				TiTv(TiTv &t);
	TiTv		&operator=(const TiTv &k);
	double		change(void);
	double		getRate(void)			{ return k; }
	double		lnProbability(void);
	void		print(void);

private:
	void		clone(const TiTv &t);
	MbRandom*	ranPtr;
	Tree*		treePtr;
	double		k;

};


class Kappa : public Parm {

public:
				Kappa(MbRandom* rp, Tree* tp, std::string pn);
				~Kappa(void);
	TiTv*		getActiveTiTv(void)		{return titvs[activeState]; }
	TiTv*		getInactiveTiTv(void)	{return titvs[getInactiveState()]; }
	double		lnPriorRatio(void);
	double		change(void);
	void		print(void);
	void		keep(void);
	void		restore(void);
	std::string	getParameterStr(void);
	std::string	getParameterHeader(void);

private:
	TiTv		*titvs[2];

};

#endif /* PARM_KAPPA_H_ */
