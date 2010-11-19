/*
 * Path.cpp
 *
 *  Created on: Jan 28, 2010
 *      Author: mikee
 */

#include "Path.h"

Step::Step() {
	time = -1.0;
	state = 0;
	site = 0;
}

Step::Step(double t, int s, int x) {

	time = t;
	state = s;
	site = x;
}

bool 	Step::operator<(const Step &A) const {

	//std::cout << "operator<\n";
	//TODO: get merge to use the operator> properly
	if (time > A.time) return true;
	else if (time == A.time && site < A.site) return true;
	else return false;
}

bool 	Step::operator>(const Step &A) const {

	//std::cout << "operator>\n";
	if (time > A.time) return true;
	else return false;
}


Path::Path()
{

}

Path::~Path()
{

}

Path::Path(Path &p) {

	clone(p);
}

void Path::clearPathHistory(void)
{

	//pathHistory.clear();
	int i = pathHistory.size() - 1;

	while (pathHistory.empty() == false)
	{
		clearSitePathHistory(i);
		pathHistory.pop_back();
		i--;
	}

}

void Path::clearSitePathHistory(int x)
{
	while (pathHistory[x].empty() == false)
	{
		pathHistory[x].pop_back();
	}
}

void Path::clone(const Path &p) {

	pathHistory = p.pathHistory;
}

Path& Path::operator=(const Path &p) {

	if (this != &p)
		clone(p);
	return *this;
}


std::vector<Step>	Path::getSeqState() {

	int i = 0;
	std::vector<Step> seqState(pathHistory.size());
	for (std::vector<std::vector<Step> >::iterator it = pathHistory.begin(); it != pathHistory.end(); it++)
	{
		seqState[i] = (*it)[0];
		i++;
	}

	return seqState;
}

void	Path::setPathSize(int numCodonSites)
{
	pathHistory.resize(numCodonSites);
}


void	Path::pushStep(Step s, int x) {

	pathHistory[x].push_back(s);
}

void	Path::addNextStep(Step s, int x) {

	// inserts an element after the first position
	std::vector<Step>::iterator it = pathHistory[x].begin();
	it++;
	pathHistory[x].insert(it, s);
}


void Path::printPathHistory()
{
	double oldTime = -1.0;

	// create vector to hold events sorted by time
	std::vector<Step>	printVector;

	// add & sort all events to printVector
	for (std::vector<std::vector<Step> >::iterator it1 = pathHistory.begin(); it1 != pathHistory.end(); it1++)
	{
		for (std::vector<Step>::iterator it2 = (*it1).begin(); it2 != (*it1).end(); it2++)
		{
			printVector.push_back(*it2);
		}
	}
	sort(printVector.begin(), printVector.end());

	// print events from printVector
	for (std::vector<Step>::iterator it = printVector.begin(); it != printVector.end(); it++)
	{
		if ((*it).getTime() == oldTime)
			std::cout << std::setw(2) << (*it).getState() << " ";
		else {
			std::cout << std::endl;
			std::cout << "t" << std::setw(24) << std::setprecision(16) << (*it).getTime() << ": ";
			for (int i = 0; i < (*it).getSite(); i++)
				std::cout << "   ";
			std::cout << std::setw(2) << (*it).getState() << " ";
		}
		oldTime = (*it).getTime();
	}
	std::cout << "\n";
}

void Path::printSitePathHistory(int x) {

	double oldTime = -1.0;

	std::vector<Step> printVector = pathHistory[x];

	// print events from printVector
	for (std::vector<Step>::iterator it = printVector.begin(); it != printVector.end(); it++)
	{
		if ((*it).getTime() == oldTime)
			std::cout << std::setw(2) << (*it).getState() << " ";
		else {
			std::cout << std::endl;
			std::cout << "t" << std::setw(24) << std::setprecision(16) << (*it).getTime() << ": ";
			for (int i = 0; i < (*it).getSite(); i++)
				std::cout << "   ";
			std::cout << std::setw(2) << (*it).getState() << " ";
		}
		oldTime = (*it).getTime();
	}
	std::cout << "\n";
}

void Path::printSeqState() {

	double oldTime = (*(*(pathHistory.begin())).begin()).getTime();
	std::cout << "t" << std::setw(24) << std::setprecision(16) << oldTime << ": ";

	for (std::vector<std::vector<Step> >::iterator it1 = pathHistory.begin(); it1 != pathHistory.end(); it1++)
	{
		Step s = *(*it1).begin();
		std::cout << std::setw(2) << s.getState() << " ";
	}
	std::cout << std::endl;
}


History::History(void) {

	activeState = 0;
	numAttemptedChanges = 0;
	numAcceptedChanges = 0;
	paths[0] = new Path();
	paths[1] = new Path(*paths[0]);
}

History::~History(void) {

	delete paths[0];
	delete paths[1];
}

void History::initPathSize(int numCodonSites)
{
	paths[0]->setPathSize(numCodonSites);
	paths[1]->setPathSize(numCodonSites);
}

void History::keep(void) {

	*paths[getInactiveState()] = *paths[activeState];
}

void History::restore(void) {

	*paths[activeState] = *paths[getInactiveState()];
}

int History::getActiveState(void) {

	if (activeState == 1)	return 1;
	else					return 0;
}

int History::getInactiveState(void) {

	if (activeState == 1)	return 0;
	else					return 1;
}


/*

void Path::printSeqAtNode() {

	double oldTime = (*(eventHistory.begin())).getTime();

	std::cout << "t" << std::setw(24) << std::setprecision(16) << oldTime << ": ";

	for (std::list<Step>::iterator it = eventHistory.begin(); it != eventHistory.end(); ++it)
	{
		if ((*it).getTime() == oldTime)
			std::cout << std::setw(2) << (*it).getState() << " ";
		oldTime = (*it).getTime();
	}
	std::cout << std::endl;
}



void Path::printEventHistory()
{
	int j = 0;
	std::cout << "PATH HISTORY!" << std::endl;
	//for (int j = 0; (unsigned int)j < eventHistory.size(); j++)  // for each step
	for (std::list<Step>::iterator it = eventHistory.begin(); it != eventHistory.end(); ++it)
	{
		std::cout << "    (step#" << j+1 << ", site#" << (*it).getSite() << ", s" << (*it).getState() << ", t" << (*it).getTime() << ", 0x= " << (*it).getNextStep() << ")";
		std::cout << std::endl;
		j++;
	}

	std::cout << std::endl;
}


*/


/*
void	Path::pushStep(Step s) {

	eventHistory.push_back(s);
}
*/
/*
void	Path::pushEventHistory(std::list<Step> p)
{
	reverse(p.begin(), p.end());
	eventHistory.merge(p);
}
*/
