/*
 * main.cpp
 *
 *  Created on: Jan 15, 2010
 *      Author: mikee
 */

#include "Mcmc.h"
#include "Alignment.h"
#include "FileMgr.h"
#include "Menu.h"
#include "MbRandom.h"
#include "ModelPath.h"
#include "ModelSeq.h"
#include "Parm_tree.h"
#include "Settings.h"

#include <iostream>

int main(int argc, char *argv[]) {

	// interpret settings
	Settings mySettings(argc, argv);
	Menu myMenu(&mySettings);

	// create file managers
	FileMgr alignFileMgr(mySettings.getAlignmentFileName());
	FileMgr treeFileMgr(mySettings.getTreeFileName());
	FileMgr outputFileMgr(mySettings.getOutputFileName());

	// read in the alignment
	Alignment myAlignment(&alignFileMgr);
	myAlignment.translate();

	//instantiate the random variable object
	MbRandom myRandom;
	//MbRandom myRandom(1);

	// read in the tree
	std::string treeStr = "";
	if (mySettings.getTreeFileName() != "")
		treeStr = treeFileMgr.readLineFromFile(1);
	std::cout << "treeStr: " << treeStr << "\n\n";

	Tree* tempTree;
	if (treeStr == "")
		tempTree = new Tree(&myRandom, "Tree", NULL, &myAlignment, (&mySettings)->getBrlenLambda());
	else
		tempTree = new Tree(&myRandom, "Tree", NULL, &myAlignment, treeStr, (&mySettings)->getBrlenLambda());
	std::cout << "\n";

	// set up the phylogenetic model
	ModelPath myModelPath(&myAlignment, &myRandom, &mySettings, tempTree);
	ModelSeq myModelSeq(&myAlignment, &myRandom, &mySettings, tempTree);

	// run the MCMC analysis
	Mcmc myMcmc( &myAlignment, &mySettings, &myRandom, &myModelPath, &myModelSeq, &outputFileMgr );

	std::cout << "origami complete!\n";
	return 0;
}
