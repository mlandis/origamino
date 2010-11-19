/*
 * Parm_omega.h
 *
 *  Created on: Jan 26, 2010
 *      Author: mikee
 */

#ifndef PARM_OMEGA_H_
#define PARM_OMEGA_H_

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
class ModelPath;

class Dnds {

public:
				Dnds(MbRandom* rp, Tree* tp);
				Dnds(Dnds &d);
	Dnds		&operator=(const Dnds &d);
	double		change(void);
	double		getDnds(void)				{ return dndsVal; }
	double		lnProbability(void);
	void		print(void);

private:
	MbRandom*	ranPtr;
	Tree*		treePtr;
	void		clone(const Dnds &d);
	double		dndsVal;

};


class Omega : public Parm {

public:
				Omega(MbRandom* rp, Tree* tp, std::string pn);
				~Omega(void);
	Dnds*		getActiveDnds(void)			{return dnds[activeState]; }
	double		lnPriorRatio(void);
	double		change(void);
	void		print(void);
	void		keep(void);
	void		restore(void);
	std::string	getParameterStr(void);
	std::string	getParameterHeader(void);

private:
	Dnds	*dnds[2];

};

#endif /* PARM_OMEGA_H_ */
