/*
 * Code.cpp
 *
 *  Created on: Jan 22, 2010
 *      Author: mikee
 */

#include "Code.h"
Code::Code(void) {

	activeCode = UNIVERSAL;
	setCode();
}

Code::Code(int whichCode) {
	activeCode = whichCode;
	setCode();
}

char Code::getAminoAcidCode(int s1, int s2, int s3) {
	int x = codon[s1*16 + s2*4 + s3];

	if (x == 1)		return 'A';
	else if (x == 2)	return 'R';
	else if (x == 3)	return 'N';
	else if (x == 4)	return 'D';
	else if (x == 5)	return 'C';
	else if (x == 6)	return 'Q';
	else if (x == 7)	return 'E';
	else if (x == 8)	return 'G';
	else if (x == 9)	return 'H';
	else if (x == 10)	return 'I';
	else if (x == 11)	return 'L';
	else if (x == 12)	return 'K';
	else if (x == 13)	return 'M';
	else if (x == 14)	return 'F';
	else if (x == 15)	return 'P';
	else if (x == 16)	return 'S';
	else if (x == 17)	return 'T';
	else if (x == 18)	return 'W';
	else if (x == 19)	return 'Y';
	else if (x == 20)	return 'V';
	else if (x == 21)	return 'X';
	else {
		std::cerr << "Unknown code: " << x << std::endl;
		exit(1);
	}
}

bool Code::isStopCodon(int s1, int s2, int s3) {
	if ( codon[s1*16 + s2*4 + s3] == 21 )
		return true;
	else
		return false;
}

void Code::setCode(void) {
	codon[ 0] = 12; /* AAA Lys K */
	codon[ 1] =  3; /* AAC Asn N */
	codon[ 2] = 12; /* AAG Lys K */
	codon[ 3] =  3; /* AAT Asn N */
	codon[ 4] = 17; /* ACA Thr T */
	codon[ 5] = 17; /* ACC Thr T */
	codon[ 6] = 17; /* ACG Thr T */
	codon[ 7] = 17; /* ACT Thr T */
	codon[ 8] =  2; /* AGA Arg R */
	codon[ 9] = 16; /* AGC Ser S */
	codon[10] =  2; /* AGG Arg R */
	codon[11] = 16; /* AGT Ser S */
	codon[12] = 10; /* ATA Ile I */
	codon[13] = 10; /* ATC Ile I */
	codon[14] = 13; /* ATG Met M */
	codon[15] = 10; /* ATT Ile I */
	codon[16] =  6; /* CAA Gln Q */
	codon[17] =  9; /* CAC His H */
	codon[18] =  6; /* CAG Gln Q */
	codon[19] =  9; /* CAT His H */
	codon[20] = 15; /* CCA Pro P */
	codon[21] = 15; /* CCC Pro P */
	codon[22] = 15; /* CCG Pro P */
	codon[23] = 15; /* CCT Pro P */
	codon[24] =  2; /* CGA Arg R */
	codon[25] =  2; /* CGC Arg R */
	codon[26] =  2; /* CGG Arg R */
	codon[27] =  2; /* CGT Arg R */
	codon[28] = 11; /* CTA Leu L */
	codon[29] = 11; /* CTC Leu L */
	codon[30] = 11; /* CTG Leu L */
	codon[31] = 11; /* CTT Leu L */
	codon[32] =  7; /* GAA Glu E */
	codon[33] =  4; /* GAC Asp D */
	codon[34] =  7; /* GAG Glu E */
	codon[35] =  4; /* GAT Asp D */
	codon[36] =  1; /* GCA Ala A */
	codon[37] =  1; /* GCC Ala A */
	codon[38] =  1; /* GCG Ala A */
	codon[39] =  1; /* GCT Ala A */
	codon[40] =  8; /* GGA Gly G */
	codon[41] =  8; /* GGC Gly G */
	codon[42] =  8; /* GGG Gly G */
	codon[43] =  8; /* GGT Gly G */
	codon[44] = 20; /* GTA Val V */
	codon[45] = 20; /* GTC Val V */
	codon[46] = 20; /* GTG Val V */
	codon[47] = 20; /* GTT Val V */
	codon[48] = 21; /* TAA Stop  */
	codon[49] = 19; /* TAC Tyr Y */
	codon[50] = 21; /* TAG Stop  */
	codon[51] = 19; /* TAT Tyr Y */
	codon[52] = 16; /* TCA Ser S */
	codon[53] = 16; /* TCC Ser S */
	codon[54] = 16; /* TCG Ser S */
	codon[55] = 16; /* TCT Ser S */
	codon[56] = 21; /* TGA Stop  */
	codon[57] =  5; /* TGC Cys C */
	codon[58] = 18; /* TGG Trp W */
	codon[59] =  5; /* TGT Cys C */
	codon[60] = 11; /* TTA Leu L */
	codon[61] = 14; /* TTC Phe F */
	codon[62] = 11; /* TTG Leu L */
	codon[63] = 14; /* TTT Phe F */

	if ( activeCode == VERTMT ) {
		//printf ("   Vertebrate mitochondrial code\n");
		//printf ("      UGA: Ter -> Trp\n");
		//printf ("      AUA: Ile -> Met\n");
		//printf ("      AGA: Arg -> Ter\n");
		//printf ("      AGG: Arg -> Ter\n");
		codon[ 8] = 21; /* AGA Stop */
		codon[10] = 21; /* AGG Stop */
		codon[12] = 13; /* ATA Met */
		codon[56] = 18; /* TGA Trp */
	}
	else if ( activeCode == MYCO ) {
		//printf ("   Mycoplasma nuclear code\n");
		//printf ("      UGA: Ter -> Trp\n");
		codon[56] = 18; /* TGA Trp */
	}
	else if ( activeCode == YEAST ) {
		//printf ("   Yeast mitochondrial code\n");
		//printf ("      UGA: Ter -> Trp\n");
		//printf ("      AUA: Ile -> Met\n");
		//printf ("      CUA: Leu -> Thr\n");
		//printf ("      CUC: Leu -> Thr\n");
		//printf ("      CUG: Leu -> Thr\n");
		//printf ("      CUU: Leu -> Thr\n");
		codon[12] = 13; /* ATA Met */
		codon[28] = 17; /* CTA Thr */
		codon[29] = 17; /* CTC Thr */
		codon[30] = 17; /* CTG Thr */
		codon[31] = 17; /* CTT Thr */
		codon[56] = 18; /* TGA Trp */
	}
	else if ( activeCode == CILIATES ) {
		//printf ("   Ciliate nuclear code\n");
		//printf ("      UAA: Ter -> Gln\n");
		//printf ("      UAG: Ter -> Gln\n");
		codon[48] =  6; /* TAA Gln */
		codon[50] =  6; /* TAG Gln */
	}
	else if ( activeCode == METMT ) {
		//printf ("   Metazoan (except vertebrates) mitochondrial code\n");
		//printf ("      UGA: Ter -> Trp\n");
		//printf ("      AUA: Ile -> Met\n");
		//printf ("      AGA: Arg -> Ser\n");
		//printf ("      AGG: Arg -> Ser\n");
		codon[ 8] = 16; /* AGA Ser */
		codon[10] = 16; /* AGG Ser */
		codon[12] = 13; /* ATA Met */
		codon[56] = 18; /* TGA Trp */
	}

	nStates = 0;
	for (int i=0; i<64; i++) {
		if ( codon[i] != 21 )
			nStates++;
	}

	int s = 0;
	for (int s1=0; s1<4; s1++) {
		for (int s2=0; s2<4; s2++) {
			for (int s3=0; s3<4; s3++) {
				if ( codon[s1*16 + s2*4 + s3] != 21 ) {
					transCodon[s] = codon[s1*16 + s2*4 + s3];
					codonNucs[s][0] = s1;
					codonNucs[s][1] = s2;
					codonNucs[s][2] = s3;
					s++;
				}
			}
		}
	}
}
