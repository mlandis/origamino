/*
 * Settings.h
 *
 *  Created on: Jan 15, 2010
 *      Author: mikee
 */

#ifndef SETTINGS_H_
#define SETTINGS_H_

#include "FileMgr.h"

#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

class FileMgr;

class Settings {
public:
				Settings(int, char**);
	std::string	getAlignmentFileName(void)			{ return alignmentFileName; }
	std::string	getTreeFileName(void)				{ return treeFileName; }
	std::string	getOutputFileName(void)				{ return outputFileName; }
	std::string	getSettingsFileName(void)			{ return settingsFileName; }
	std::string getCodonMatrixFileName(void)		{ return codonMatrixFileName; }
	std::string getCodonFreqFileName(void)         { return codonFreqFileName; }
	std::string getCoordinatesFileName(void)		{ return coordinatesFileName; }
	std::string getFreeEnergyFileName(void)			{ return freeEnergyFileName; }
	double		getBrlenLambda(void)				{ return brlenLambda; }
	double		getContactDistance(void)			{ return contactDistance; }
	int			getSiteChangesPerCycle(void)		{ return siteChangesPerCycle; }
	int			getChainLength(void)				{ return numCycles; }
	int			getPrintFrequency(void)				{ return printFrequency; }
	int			getSampleFrequency(void)			{ return sampleFrequency; }
	int			getSummarizeFrequency(void)			{ return summarizeFrequency; }
	int			getSVESamplesPerCycle(void)			{ return SVESamplesPerCycle; }
	int			getNumOmegaCats(void)				{ return numOmegaCats; }
	bool		getUseECM(void)						{ return useECM; }
	void		setAlignmentFileName(std::string s)	{ alignmentFileName = s; }
	void		setTreeFileName(std::string s)		{ treeFileName = s; }
	void		setOutputFileName(std::string s)	{ outputFileName = s; }
	void		setSettingsFileName(std::string s)	{ settingsFileName = s; }
	void		setCodonMatrixFileName(std::string s)	{ codonMatrixFileName = s; }
	void        setCodonFreqFileName(std::string s)    { codonFreqFileName = s; }
	void		setCoordinatesFileName(std::string s)	{ coordinatesFileName = s; }
	void		setFreeEnergyFileName(std::string s)	{ freeEnergyFileName = s; }
	void		setBrlenLambda(double x)			{ brlenLambda = x; }
	void		setContactDistance(double x)		{ contactDistance = x; }
	void		setSiteChangesPerCycle(int x)		{ siteChangesPerCycle = x; }
	void		setSVESamplesPerCycle(int x)		{ SVESamplesPerCycle = x; }

private:
	std::string	alignmentFileName;
	std::string	treeFileName;
	std::string	outputFileName;
	std::string	settingsFileName;
	std::string codonMatrixFileName;
	std::string codonFreqFileName;
	std::string coordinatesFileName;
	std::string freeEnergyFileName;
	double		brlenLambda;
	double		contactDistance;
	int			numOmegaCats;
	int			numCycles;
	int			printFrequency;
	int			sampleFrequency;
	int			summarizeFrequency;
	int			SVESamplesPerCycle;
	int			siteChangesPerCycle;
	bool		useECM;
};

#endif /* SETTINGS_H_ */
