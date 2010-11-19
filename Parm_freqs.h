/*
 * Parm_freqs.h
 *
 *  Created on: Jan 26, 2010
 *      Author: mikee
 */

#ifndef PARM_FREQS_H_
#define PARM_FREQS_H_

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "Alignment.h"
#include "MbRandom.h"
#include "Parm.h"
#include "Parm_tree.h"

class Alignment;
class MbRandom;
class ModelPath;
class Tree;

class Freqs {

public:
							Freqs(MbRandom *rp, Tree *tp, Alignment *ap);
							Freqs(std::vector<double> values, MbRandom *rp, Tree *tp, Alignment *ap);
							Freqs(Freqs &c);
	Freqs					&operator=(const Freqs &A);
	double					change(void);
	void					setFreqs(std::vector<double> x);
	std::vector<double>&	getFreqs(void)				{ return f; }
	double					getFreqVal(int x)			{ return f[x]; }
	double					lnProbability(void);
	void					print(void);

private:
	MbRandom*				ranPtr;
	Tree*					treePtr;
	void					clone(const Freqs &c);
	void					normalizeFreqs(std::vector<double>& v, double threshold);
	int						numCodons;
	std::vector<double>		a;
	std::vector<double>		f;

};

class CodonFrequencies : public Parm {

public:
							CodonFrequencies(MbRandom *rp, std::string pn, Tree *tp, Alignment *ap);
							CodonFrequencies(std::vector<double> values, MbRandom *rp, std::string pn, Tree *tp, Alignment *ap);
							~CodonFrequencies(void);
	Freqs*					getActiveFreqs(void)	{ return freqs[activeState]; }
	double					lnPriorRatio(void);
	double					change(void);
	void					print(void);
	void					keep(void);
	void					restore(void);
	std::string				getParameterStr(void);
	std::string				getParameterHeader(void);

private:
	Freqs					*freqs[2];

};

#endif /* PARM_FREQS_H_ */
