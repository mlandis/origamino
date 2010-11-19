/*
 * Code.h
 *
 *  Created on: Jan 22, 2010
 *      Author: mikee
 */

#ifndef CODE_H_
#define CODE_H_

#include <stdlib.h>
#include <iostream>

class Code {

	public:
					Code(void);
					Code(int whichCode);
		int			getAminoAcid(int s1, int s2, int s3) 		{ return codon[s1*16 + s2*4 + s3]; }
		char		getAminoAcidCode(int s1, int s2, int s3);
		int			getCodonNuc(int c, int i) 					{ return codonNucs[c][i]; }
		int			getTransCodon(int c) 						{ return transCodon[c]; }
		bool		isStopCodon(int s1, int s2, int s3);
		int			getNumSenseCodons(void) 					{ return nStates; }

	private:
		void		setCode(void);
		int			codon[64];                                                          /* holds amino acid identity for each of the 64 codons  */
		int			codonNucs[64][3];                                                   /* nucleotides for the coding codons                    */
		int			transCodon[64];                                                     /* which codons are coding (i.e., not stop codons)      */
		int			activeCode;                                                         /* which genetic code is in use                         */
		enum		activeCode { UNIVERSAL, VERTMT, MYCO, YEAST, CILIATES, METMT };     /* the possible genetic codes                           */
		int			nStates;

};

#endif /* CODE_H_ */
