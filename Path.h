/*
 * Path.h
 *
 *  Created on: Jan 28, 2010
 *      Author: mikee
 */

#ifndef PATH_H_
#define PATH_H_

#include "MbRandom.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <list>
#include <vector>

class Step {

public:
						Step(void);
						Step(double t, int s, int x);
	bool				operator<(const Step &A) const;
	bool				operator>(const Step &A) const;
	void				setNextStep(Step* s);
	double				getTime(void)					{ return time; }
	int					getState(void)					{ return state; }
	int					getSite(void)					{ return site; }

private:
	double				time;		// time step occurred
	int					state;		// state acquired from step (e.g. ACT -> GCC; stepState = codonID(GCC))
	int					site;

};

class Path {

public:
										Path();
										~Path();
										Path(Path &p);
	std::vector<std::vector<Step> >		getPathHistory()					{ return pathHistory; }
	std::vector<Step> 					getSitePathHistory(int x)			{ return pathHistory[x]; }
	std::vector<Step>					getSeqState();
	Path								&operator=(const Path &d);
	double								change(void);
	void								setPathSize(int numCodonSites);
	void								pushStep(Step s, int x);
	void								addNextStep(Step s, int x);
	void								clearPathHistory(void);
	void								clearSitePathHistory(int x);
	void								printPathHistory(void);
	void								printSitePathHistory(int x);			// prints event history for a particular site
	void								printSeqState(void);


private:
	void								clone(const Path &p);
	std::vector< std::vector <Step> >	pathHistory;

};

class History {

public:
				History(void);
				~History(void);
				History(MbRandom* rp, std::string pn);
	Path*		getActivePath(void)							{return paths[activeState]; }
	Path*		getInactivePath(void)						{return paths[getInactiveState()]; }
	int			getActiveState(void);
	int			getInactiveState(void);
	double		lnPriorRatio(void);
	double		change(void);
	void		print(void);
	void		keep(void);
	void		restore(void);
	void		incrementNumAccepted(void)					{ numAcceptedChanges++; }
	int			getNumAcceptances(void)						{ return numAcceptedChanges; }
	int			getNumAttempts(void)						{ return numAttemptedChanges; }
	std::string	getParameterStr(void);
	std::string	getParameterHeader(void);
	void		initPathSize(int numCodonSites);

private:
	Path		*paths[2];
	bool		activeState;
	int			numAttemptedChanges;
	int			numAcceptedChanges;

};

#endif /* PATH_H_ */


// GRAVEYARD

/*
 * public:
	Step				getStep(int site, int step)			{ return pathHistory[site][step]; }
	double				getStepTime(int site, int step)		{ return pathHistory[site][step].getTime(); }
	int					getStepState(int site, int step)	{ return pathHistory[site][step].getState(); }
	int					getStepCount(int site)				{ return pathHistory[site].size(); }
	int					getSiteCount(void)					{ return pathHistory.size(); }
	void				clearPathHistory(void)				{ pathHistory.clear(); }
	void				clearStepHistory(int site)			{ pathHistory[site].clear(); }
	*/
	//std::vector<Step>	getPath(int site)					{ return pathHistory[site]; }
	//void				pushStep(int site, Step s);
	//void				pushStep(int site, double t, int s, int x);
	//void				pushPath(std::vector<Step> p);
	//void				popStep(int site)					{ pathHistory[site].pop_back(); }
	//void				popPath(void)						{ pathHistory.pop_back(); }

// private::
// std::vector < std::vector<Step> >	pathHistory;	// should now be obsolete
