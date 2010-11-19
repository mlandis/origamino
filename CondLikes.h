#ifndef CONDLIKES_H
#define CONDLIKES_H

#include <iomanip>
#include <iostream>
#include <vector>
#include "Alignment.h"
#include "MbBitfield.h"

class Alignment;
class CondLikes {

	public:
                            CondLikes(Alignment *ap);
							~CondLikes(void);
				  double*   getClPtr(int space, int node, int site) { return clsPtr[space][node][site]; }
				  double*   getClScalerPtr(int space, int node, int site) { return clScalerPtr[space][node][site]; }
				  double*   getClPtrUpL(int space, int node, int site) { return clsPtrUpL[space][node][site]; }
				  double*   getClScalerPtrUpL(int space, int node, int site) { return clScalerPtrUpL[space][node][site]; }
				  double*   getClPtrUpR(int space, int node, int site) { return clsPtrUpR[space][node][site]; }
				  double*   getClScalerPtrUpR(int space, int node, int site) { return clScalerPtrUpR[space][node][site]; }
					 void   print(void);

	private:
	            Alignment   *alignmentPtr;
				      int   numTaxa;
					  int   numSites;
					  int   numNodes;
					  int   numStates;
					  
	              double*   cls;
				double***   clsPtr[2];
				  double*   clScalers;
				double***   clScalerPtr[2];
				
	              double*   clsUpL;
				double***   clsPtrUpL[2];
				  double*   clScalersUpL;
				double***   clScalerPtrUpL[2];

	              double*   clsUpR;
				double***   clsPtrUpR[2];
				  double*   clScalersUpR;
				double***   clScalerPtrUpR[2];
};


#endif
