/*
 * Alignment.h
 *
 *  Created on: Jan 18, 2010
 *      Author: mikee
 */

#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include "MbBitfield.h"
#include "Code.h"
#include "FileMgr.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

class MbBitfield;
class Code;
class FileMgr;

class Alignment {
public:
								Alignment(FileMgr *file);
								~Alignment(void);

	MbBitfield*					getCodonMatrixEntry(int t, int s)		{ return codonMatrix[t][s]; }
	int							getOrigMatrixEntry(int t, int s)		{ return originalMatrix[t][s]; }
	std::string					getNameForTaxon(int i)					{ return taxonNames[i]; }
	int							getContextID(int t, int pos);
	int							getCodonID(int t, int pos);
	int							getIndexForTaxon(std::string s);
	int   						getNumTaxa(void)						{ return numTaxa; }
	int							getNumSites(void)						{ return numSites; }
	int							getNumCodonSites(void)					{ return numCodonSites; }
	int							getNumSenseCodons(void);
	int							nucID(char nuc);
	Code*						getGeneticCode(void)					{ return geneticCode; }
	void						print(void);
	bool						translate(void);

private:
	void						getPossibleNucs(int nucCode, int nuc[]);
	bool						readAlignment(FileMgr *file);
	std::vector<std::string>	taxonNames;
	int							**originalMatrix;
	MbBitfield					***codonMatrix;
	Code						*geneticCode;
	bool						isTranslated;
	int							numTaxa;
	int							numSites;
	int							numCodonSites;



};

#endif /* ALIGNMENT_H_ */
