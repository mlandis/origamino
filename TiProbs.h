#ifndef TIPROBS_H
#define TIPROBS_H

#include "MbMatrix.h"
#include "Parm_scaling.h"

#include <vector>

class MbTransitionMatrix;
class SubstRate;

class TiProbs {

	public:
                            TiProbs(int nn, MbTransitionMatrix* tm, Tree* tp);
							~TiProbs(void);
		MbMatrix<double>&   getTiMatrix(int space, int node) { return tis[space][node]; }
					 void   setTransitionProbs(void);
					 void   print(void);

	private:
					 bool   isMatrixGood(int space, int node, double thresshold);
					  int   numNodes;
					  int   numStates;
	  MbTransitionMatrix*   tiMatrixPtr;
				  Tree*		treePtr;
		std::vector<double>	substVals;
	     MbMatrix<double>   *tis[2];
};


#endif
