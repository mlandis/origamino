/*
 * Parm_solv.cpp
 *
 *  Created on: Apr 22, 2010
 *      Author: mikee
 */

#include "Parm_solv.h"

SolvAcc::SolvAcc(MbRandom* rp, Tree* tp) {

	ranPtr = rp;
	treePtr = tp;
	// TODO: solvacc currently disabled
	//	solvaccVal = ranPtr->exponentialRv(1.0);
	solvaccVal = 0.0;
	print();
}

SolvAcc::SolvAcc(SolvAcc &d) {

	clone(d);
}

SolvAcc& SolvAcc::operator=(const SolvAcc &d) {

	if (this != &d)
		clone(d);
	return *this;
}

double SolvAcc::change(void) {

	// update conditional likelihood and transition probability flags
	Topology *t = treePtr->getActiveTopology();
	t->updateAllCls(true);
	t->updateAllTis(true);
	t->flipAllActiveCls();
	t->flipAllActiveTis();

	double tuning = log(4.0);
	double oldS = solvaccVal;
	double newS = oldS * exp(tuning * (ranPtr->uniformRv() - 0.5));
	//solvaccVal = newS;

	// TODO: solvacc values currently set to 0.0
	// commit the change
	solvaccVal = 0.0;
	return 1.0;
	return log(newS) - log(oldS);
}

double SolvAcc::lnProbability(void) {

	return -2.0 * log(1.0 + solvaccVal);
}

void SolvAcc::print(void) {

	std::cout << "SolvAcc: ";
	std::cout << std::fixed << std::setprecision(10) << solvaccVal;
	std::cout << std::endl;
}

void SolvAcc::clone(const SolvAcc &d) {

	ranPtr = d.ranPtr;
	treePtr = d.treePtr;
	solvaccVal = d.solvaccVal;
}

SeqStructS::SeqStructS(MbRandom* rp, Tree* tp, std::string pn) : Parm(rp, pn) {

	solvacc[0] = new SolvAcc(rp, tp);
	solvacc[1] = new SolvAcc(*solvacc[0]);
}

SeqStructS::~SeqStructS(void) {

	delete solvacc[0];
	delete solvacc[1];
}

double SeqStructS::lnPriorRatio(void) {

	return solvacc[activeState]->lnProbability() - solvacc[getInactiveState()]->lnProbability();
}

double SeqStructS::change(void) {

	numAttemptedChanges++;
	return solvacc[activeState]->change();
}

void SeqStructS::print(void) {

	solvacc[activeState]->print();
}

void SeqStructS::keep(void) {

	*solvacc[getInactiveState()] = *solvacc[activeState];
}

void SeqStructS::restore(void) {

	*solvacc[activeState] = *solvacc[getInactiveState()];
}

std::string SeqStructS::getParameterStr(void) {

	double w = solvacc[activeState]->getSolvAcc();
	char tempCh[50];
	sprintf(tempCh, "%1.3lf\t", w);
	std::string pStr = tempCh;
	return pStr;
}

std::string SeqStructS::getParameterHeader(void) {

	std::string pStr = "SeqStructS\t";
	return pStr;
}
