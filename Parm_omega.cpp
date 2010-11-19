/*
 * Parm_omega.cpp
 *
 *  Created on: Jan 26, 2010
 *      Author: mikee
 */

#include "Parm_omega.h"

Dnds::Dnds(MbRandom* rp, Tree* tp) {

	ranPtr = rp;
	treePtr = tp;
	double dS = ranPtr->exponentialRv(1.0);
	dndsVal = ranPtr->exponentialRv(1.0) / dS;
	print();
}

Dnds::Dnds(Dnds &d) {

	clone(d);
}

Dnds& Dnds::operator=(const Dnds &d) {

	if (this != &d)
		clone(d);
	return *this;
}

double Dnds::change(void) {

	// update conditional likelihood and transition probability flags
	Topology *t = treePtr->getActiveTopology();
	t->updateAllCls(true);
	t->updateAllTis(true);
	t->flipAllActiveCls();
	t->flipAllActiveTis();

	// set the maximum distance of the move
	double x = 0.1;

	// propose new values for omega
	double newVal = 0.0;
	double oldVal = 0.0;
	double offset = 0.0;

	do {
		// pick the value
		offset = x * (ranPtr->uniformRv() - 0.5);

		// calculate the new point
		newVal = dndsVal + offset;

	} while (newVal < 0.0); // decide if the new point is valid

	// commit the change
	oldVal = dndsVal;
	dndsVal = newVal;

	return log(newVal) - log(oldVal);
}

double Dnds::lnProbability(void) {

	return -2.0 * log(1.0 + dndsVal);
}

void Dnds::print(void) {

	std::cout << "Omega: ";
	std::cout << std::fixed << std::setprecision(10) << dndsVal;
	std::cout << std::endl;
}

void Dnds::clone(const Dnds &d) {

	ranPtr = d.ranPtr;
	treePtr = d.treePtr;
	dndsVal = d.dndsVal;
}

Omega::Omega(MbRandom* rp, Tree* tp, std::string pn) : Parm(rp, pn) {

	dnds[0] = new Dnds(rp, tp);
	dnds[1] = new Dnds(*dnds[0]);
}

Omega::~Omega(void) {

	delete dnds[0];
	delete dnds[1];
}

double Omega::lnPriorRatio(void) {

	return dnds[activeState]->lnProbability() - dnds[getInactiveState()]->lnProbability();
}

double Omega::change(void) {

	numAttemptedChanges++;
	return dnds[activeState]->change();
}

void Omega::print(void) {

	dnds[activeState]->print();
}

void Omega::keep(void) {

	*dnds[getInactiveState()] = *dnds[activeState];
}

void Omega::restore(void) {

	*dnds[activeState] = *dnds[getInactiveState()];
}

std::string Omega::getParameterStr(void) {

	double w = dnds[activeState]->getDnds();
	char tempCh[50];
	sprintf(tempCh, "%1.3lf\t", w);
	std::string pStr = tempCh;
	return pStr;
}

std::string Omega::getParameterHeader(void) {

	std::string pStr = "Omega\t";
	return pStr;
}
