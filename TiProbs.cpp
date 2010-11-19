#include <iomanip>
#include <iostream>
#include "MbTransitionMatrix.h"
#include "ModelPath.h"
#include "Parm_tree.h"
#include "TiProbs.h"



TiProbs::TiProbs(int nn, MbTransitionMatrix *tm, Tree* tp) {

	numNodes    	= nn;
	tiMatrixPtr 	= tm;
	treePtr 		= tp;
	numStates   	= tiMatrixPtr->getNumStates();
	
	for (int n = 0; n < 2; n++)
	{
		tis[n] = new MbMatrix<double>[numNodes];
		for (int i = 0; i < numNodes; i++)
			tis[n][i] = MbMatrix<double>(numStates, numStates);
	}
}

TiProbs::~TiProbs(void) {

	for (int n = 0; n < 2; n++)
		delete [] tis[n];
}

bool TiProbs::isMatrixGood(int space, int node, double thresshold) {

	MbMatrix<double> a;
	a = tis[space][node];
	bool isGood = true;
	for (int i = 0; i < numStates; i++)
	{
		double rowSum = 0.0;
		for (int j = 0; j < numStates; j++)
			rowSum += a[i][j];
		if ( rowSum > 1.0 + thresshold || rowSum < 1.0 - thresshold )
			isGood = false;
	}
	return isGood;
}

void TiProbs::print(void) {

#	if 1
	int numCols = numStates;
	numCols = 50;
	Topology* t = treePtr->getActiveTopology();
	for (int n=0; n<t->getNumNodes(); n++)
		{
		Node* p = t->getDownPassNode(n);
		if (p->getAnc() != NULL)
			{
			std::cout << "Matrix " << n << ":" << std::endl;
			std::cout << std::fixed << std::setprecision(4);
			for (int i=0; i<numStates; i++)
				{
				for (int j=0; j<numCols; j++)
					std::cout << tis[p->getActiveTi()][p->getIndex()][i][j] << " ";
				std::cout << std::endl;
				}
			}
		}
#	else
	for (int n=0; n<numNodes; n++)
	{
		std::cout << "Matrix " << n << ":" << std::endl;
		std::cout << std::fixed << std::setprecision(3);
		for (int i=0; i<numStates; i++)
		{
			for (int j=0; j<numStates; j++)
				std::cout << tis[0][n][i][j] << " ";
			std::cout << " -- ";
			for (int j=0; j<numStates; j++)
				std::cout << tis[1][n][i][j] << " ";
			std::cout << std::endl;
		}
	}
#	endif
}

void TiProbs::setTransitionProbs(void) {

	Topology* t = treePtr->getActiveTopology();
	MbMatrix<double> a;
	
	for (int n = 0; n < t->getNumNodes(); n++)
	{
		Node* p = t->getDownPassNode(n);
		if (p->getUpdateTi() == true)
		{
			a = getTiMatrix( p->getActiveTi(), p->getIndex() );
			a = tiMatrixPtr->tiProbs( p->getV(), a );
			p->setUpdateTi(false);
#			if 0
			if (isMatrixGood( p->getActiveTi(), p->getIndex(), 0.000001) == false)
			{
				std::cerr << "ERROR: Bad transition probability matrix!" << std::endl;
				exit(1);
			}
#			endif
		}
	}
}
