/*
 * Parm_pairwise.cpp
 *
 *  Created on: Apr 22, 2010
 *      Author: mikee
 */

#include "Parm_pairwise.h"

Pairwise::Pairwise(MbRandom* rp, Tree* tp) {

	ranPtr = rp;
	treePtr = tp;
	// TODO: pairwise values currently set to 0.0
	// pairwiseVal = ranPtr->exponentialRv(1.0);
	pairwiseVal = 0.0;
	print();
}

Pairwise::Pairwise(Pairwise &d) {

	clone(d);
}

Pairwise& Pairwise::operator=(const Pairwise &d) {

	if (this != &d)
		clone(d);
	return *this;
}

double Pairwise::change(void) {

	// update conditional likelihood and transition probability flags

	Topology* t = treePtr->getActiveTopology();
	t->updateAllCls(true);
	t->updateAllTis(true);
	t->flipAllActiveCls();
	t->flipAllActiveTis();

	double tuning = log(4.0);
	double oldP = pairwiseVal;
	double newP = oldP * exp(tuning * (ranPtr->uniformRv() - 0.5));
	//pairwiseVal = newP;

	// TODO: pairwise values currently set to 0.0
	// commit the change
	pairwiseVal = 0.0;
	return 1.0;
	return log(newP) - log(oldP);
}

double Pairwise::lnProbability(void) {

	return -2.0 * log(1.0 + pairwiseVal);
}

void Pairwise::print(void) {

	std::cout << "Pairwise: ";
	std::cout << std::fixed << std::setprecision(10) << pairwiseVal;
	std::cout << std::endl;
}

void Pairwise::clone(const Pairwise &d) {

	ranPtr = d.ranPtr;
	treePtr = d.treePtr;
	pairwiseVal = d.pairwiseVal;
}

SeqStructP::SeqStructP(MbRandom* rp, Tree* tp, std::string pn) : Parm(rp, pn) {

	pairwise[0] = new Pairwise(rp, tp);
	pairwise[1] = new Pairwise(*pairwise[0]);
}

SeqStructP::~SeqStructP(void) {

	delete pairwise[0];
	delete pairwise[1];
}

double SeqStructP::lnPriorRatio(void) {

	return pairwise[activeState]->lnProbability() - pairwise[getInactiveState()]->lnProbability();
}

double SeqStructP::change(void) {

	numAttemptedChanges++;
	return pairwise[activeState]->change();
}

void SeqStructP::print(void) {

	pairwise[activeState]->print();
}

void SeqStructP::keep(void) {

	*pairwise[getInactiveState()] = *pairwise[activeState];
}

void SeqStructP::restore(void) {

	*pairwise[activeState] = *pairwise[getInactiveState()];
}

std::string SeqStructP::getParameterStr(void) {

	double w = pairwise[activeState]->getPairwise();
	char tempCh[50];
	sprintf(tempCh, "%1.3lf\t", w);
	std::string pStr = tempCh;
	return pStr;
}

std::string SeqStructP::getParameterHeader(void) {

	std::string pStr = "SeqStructP\t";
	return pStr;
}
