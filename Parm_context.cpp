/*
 * Parm_context.cpp
 *
 *  Created on: Jan 26, 2010
 *      Author: mikee
 */

#include "Parm_context.h"

Contexts::Contexts(MbRandom* rp, Tree* tp) {

	ranPtr = rp;
	treePtr = tp;
	numContexts = 16;
	a.resize(numContexts);
	contextProbs.resize(numContexts);
	contextNames.resize(numContexts);
	for (int i = 0; i < numContexts; i++)
			a[i] = 1.0;
	ranPtr->dirichletRv(a, contextProbs);
	setContextNames();
	print();
}

Contexts::Contexts(Contexts &c) {

	clone(c);
}

Contexts& Contexts::operator=(const Contexts &ctxts) {

	if (this != &ctxts)
			clone(ctxts);
	return *this;
}

double Contexts::change(void) {

	// update conditional likelihood and transition probability flags
	Topology* t = treePtr->getActiveTopology();
	t->updateAllCls(true);
	t->updateAllTis(true);
	t->flipAllActiveCls();
	t->flipAllActiveTis();

	// tuning parameters out of the move
	int numChangedContexts = 2;		// how many contexts to change in sequence!
	double alpha0 = 500.0;

	// pick some frequencies that will be changed
	std::vector<int> temp(numContexts);
	for (int i = 0; i < numContexts; i++)
		temp[i] = i;

	std::vector<int> contextIds(numChangedContexts);
	for (int i = 0; i < numChangedContexts; i++) {
		int whichContext = (int)(ranPtr->uniformRv() * (numContexts - i));
		contextIds[i] = temp[whichContext];
		temp[numContexts - i - 1] = contextIds[i];
	}

	// fill in a vector of the current Contexts
	std::vector<double> curF(numChangedContexts + 1);
	std::vector<double> newF(numChangedContexts + 1);
	for (int i = 0; i < numChangedContexts + 1; i++) {
		curF[i] = 0.0;
		newF[i] = 0.0;
	}
	for (int i = 0; i < numContexts; i++) {
		// is the context one of the context to have its frequency changed?
		bool isContextSpecial = false;
		int whichOffsetInContextIds;
		for (int j = 0; j < numChangedContexts; j++) {
			if (contextIds[j] == i) {
				isContextSpecial = true;
				whichOffsetInContextIds = j;
				break;
			}
		}

		// fill in the context frequency
		if (isContextSpecial == true)
			curF[whichOffsetInContextIds] = contextProbs[i];
		else
			curF[numChangedContexts] += contextProbs[i];
	}

	// fill in the Dirichlet prior parameter
	std::vector<double> forwardA(numChangedContexts + 1);
	std::vector<double> reverseA(numChangedContexts + 1);
	for (int i = 0; i < numChangedContexts + 1; i++)
		forwardA[i] = curF[i] * alpha0;

	// choose new context frequencies
	ranPtr->dirichletRv(forwardA, newF);
	for(int i = 0; i < numChangedContexts + 1; i++)
		reverseA[i] = newF[i] * alpha0;
	normalizeContextProbs(newF, 0.000001);

	// fill in the new context frequencies
	for (int i = 0; i < numContexts; i++) {
		// is the context one of the context to have its frequency changed?
		bool isContextSpecial = false;
		int whichOffsetInContextIds;
		for (int j = 0; j < numChangedContexts; j++) {
			if (contextIds[j] == i) {
				isContextSpecial = true;
				whichOffsetInContextIds = j;
				break;
			}
		}

		// fill in the context frequency
		if (isContextSpecial == true)
			contextProbs[i] = newF[whichOffsetInContextIds];
		else
			contextProbs[i] = newF[numChangedContexts] / curF[numChangedContexts];
	}

	double x = ranPtr->lnDirichletPdf(reverseA, curF) - ranPtr->lnDirichletPdf(forwardA, newF) + (numContexts - numChangedContexts) * log(newF[numChangedContexts] / curF[numChangedContexts]);
	return x;
}

double Contexts::lnProbability(void) {

	return ranPtr->lnDirichletPdf(a, contextProbs);
}

void Contexts::normalizeContextProbs(std::vector<double>& v, double threshhold) {

	double sum1 = 0.0;
	double sum2 = 0.0;
	for (int i = 0; i < (int)v.size(); i++) {
		if (v[i] >= threshhold)
			sum1 += v[i];
		else
			sum2 += threshhold;
	}
	double scaler = (1.0 - sum2) / sum1;
	for (int i = 0; i < (int)v.size(); i++) {
		if (v[i] < threshhold)
			v[i] = threshhold;
		else
			v[i] *= scaler;
	}
}

void Contexts::setContextNames() {

	// Indexed as base-4 numbers, where A = 0, C = 1, G = 2, T = 3.
	//   e.g. - CxT = 4(1) + 1(3) =  7
	//   e.g. - TxG = 4(3) + 1(2) = 14

	// immediate nucleotide neighbor context dependency
	if (numContexts == 16) {
		contextNames[ 0] = "AxA";
		contextNames[ 1] = "AxC";
		contextNames[ 2] = "AxG";
		contextNames[ 3] = "AxT";
		contextNames[ 4] = "CxA";
		contextNames[ 5] = "CxC";
		contextNames[ 6] = "CxG";
		contextNames[ 7] = "CxT";
		contextNames[ 8] = "GxA";
		contextNames[ 9] = "GxC";
		contextNames[10] = "GxG";
		contextNames[11] = "GxT";
		contextNames[12] = "TxA";
		contextNames[13] = "TxC";
		contextNames[14] = "TxG";
		contextNames[15] = "TxT";
	}
	else {
			std::cerr << "ERROR: Unknown context dependence schema!" << std::endl;
			exit(1);
	}
}

void Contexts::print(void) {

	std::cout << "Context Dependencies: ";
	for (int i = 0; i < numContexts; i++) {
		std::cout << contextNames[i] << ":";
		std::cout << std::fixed << std::setprecision(10) << contextProbs[i] << " ";
	}
	std::cout << std::endl;
}

void Contexts::clone(const Contexts &c) {

	ranPtr		= c.ranPtr;
	treePtr		= c.treePtr;
	numContexts	= c.numContexts;
	contextProbs.resize(c.contextProbs.size());
	a.resize(c.a.size());
	for (int i = 0; i < numContexts; i++) {
		a[i] = c.a[i];
		contextProbs[i] = c.contextProbs[i];
	}
}

ContextDependence::ContextDependence(MbRandom* rp, Tree* tp, std::string pn) : Parm(rp, pn) {

	contexts[0] = new Contexts(rp, tp);
	contexts[1] = new Contexts(*contexts[0]);
}

ContextDependence::~ContextDependence(void) {

	delete contexts[0];
	delete contexts[1];
}

double ContextDependence::lnPriorRatio(void) {

	return contexts[activeState]->lnProbability() - contexts[getInactiveState()]->lnProbability();
}

double ContextDependence::change(void) {

	numAttemptedChanges++;
	return contexts[activeState]->change();
}

void ContextDependence::print(void) {

	return contexts[activeState]->print();
}

void ContextDependence::keep(void) {

	*contexts[getInactiveState()] = *contexts[activeState];
}

void ContextDependence::restore(void) {

	*contexts[activeState] = *contexts[getInactiveState()];
}

std::string ContextDependence::getParameterStr(void) {

	std::vector<double> cp = contexts[activeState]->getContexts();
	std::string pStr = "";
	for(int i = 0; i < (int)cp.size(); i++) {
		char tempCh[50];
		sprintf(tempCh, "%1.3lf\t", cp[i]);
		std::string tempStr = tempCh;
		pStr += tempStr;
	}
	return pStr;
}

std::string ContextDependence::getParameterHeader(void) {

	std::vector<double> cp = contexts[activeState]->getContexts();
	std::string pStr = "";
	for (int i = 0; i < (int)cp.size(); i++) {
		char tempCh[50];
		sprintf(tempCh, "Pi(%d)\t", i + 1);
		std::string tempStr = tempCh;
		pStr += tempStr;
	}
	return pStr;
}
