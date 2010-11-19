/*
 * Parm_scaling.cpp
 *
 *  Created on: May 4, 2010
 *      Author: mikee
 */

#include "Parm_scaling.h"

SubstRate::SubstRate(MbRandom* rp, Tree* tp, double ep) {

	ranPtr = rp;
	treePtr = tp;
	brlenLambda = ep;
	Topology* t = treePtr->getActiveTopology();
	u.resize(t->getNumNodes());

	// aligned with downpassnode sequence only when the root is skipped (e.g. for p->getAnc() == NULL)
	for (int n = 0; n < t->getNumNodes(); n++)
	{
		Node* p = t->getDownPassNode(n);
		if (p->getAnc() != NULL)
			u[n] = 1.0 / p->getV();
	}

	print();
}

SubstRate::SubstRate(SubstRate &d) {

	clone(d);
}

SubstRate& SubstRate::operator=(const SubstRate &d) {

	if (this != &d)
		clone(d);
	return *this;
}

double SubstRate::change(void) {

	// update conditional likelihood and transition probability flags
	// TODO: probably don't need these
	Topology *t = treePtr->getActiveTopology();
	t->updateAllCls(true);
	t->updateAllTis(true);
	t->flipAllActiveCls();
	t->flipAllActiveTis();

	int n = 0;
	do
	{
		n = (int)(ranPtr->uniformRv() * u.size());
	}
	while(t->getDownPassNode(n)->getAnc() == NULL);

	double newU = 0.0;
	double oldU = u[n];
	double tuning = log(4.0);
	do {
		// pick the value
		newU = oldU * exp(tuning * (ranPtr->uniformRv() - 0.5));

	} while (newU < 0.0); // decide if the new point is valid

	// commit the change
	u[n] = newU;

//	std::cout << "n: " << n << ", oldU: " << oldU << ", newU: " << newU << "\n";
	return log(newU) - log(oldU);
}

double SubstRate::lnProbability(void) {

	double lnP = 0.0;
	Topology* t = treePtr->getActiveTopology();

	for (int n = 0; n < t->getNumNodes(); n++) {
		lnP += ranPtr->lnExponentialPdf(brlenLambda, u[n]);
	}

	//print();
	return lnP;
}

void SubstRate::print(void) {

	std::cout << "Substitution Rates (indexed by DownPassNode): ";
	for (unsigned int i = 0; i < u.size(); i++)
		std::cout << std::setw(2) << i << ":" << std::fixed << std::setprecision(10) << u[i] << " ";
	std::cout << std::endl;
}

void SubstRate::clone(const SubstRate &d) {

	ranPtr = d.ranPtr;
	treePtr = d.treePtr;
	brlenLambda = d.brlenLambda;
	u = d.u;
}

ScalerU::ScalerU(MbRandom* rp, Tree* tp, double ep, std::string pn) : Parm(rp, pn) {

	substRate[0] = new SubstRate(rp, tp, ep);
	substRate[1] = new SubstRate(*substRate[0]);
}

ScalerU::~ScalerU(void) {

	delete substRate[0];
	delete substRate[1];
}

double ScalerU::lnPriorRatio(void) {

	return substRate[activeState]->lnProbability() - substRate[getInactiveState()]->lnProbability();
}

double ScalerU::change(void) {

	numAttemptedChanges++;
	return substRate[activeState]->change();
}

void ScalerU::print(void) {

	substRate[activeState]->print();
}

void ScalerU::keep(void) {

	*substRate[getInactiveState()] = *substRate[activeState];
}

void ScalerU::restore(void) {

	*substRate[activeState] = *substRate[getInactiveState()];
}

std::string ScalerU::getParameterStr(void) {

	std::vector<double> f = substRate[activeState]->getSubstRate();
	std::string pStr = "";
	for(int i = 0; i < (int)f.size(); i++)
	{
		char tempCh[50];
		sprintf(tempCh, "%1.3lf\t", f[i]);
		std::string tempStr = tempCh;
		pStr += tempStr;
	}
	return pStr;
}

std::string ScalerU::getParameterHeader(void) {

	std::string pStr = "SubstRate\t";
	return pStr;
}
