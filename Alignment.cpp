/*
 * Alignment.cpp
 *
 *  Created on: Jan 18, 2010
 *      Author: mikee
 */

#include "Alignment.h"

Alignment::Alignment(FileMgr *file) {

	if (readAlignment(file) == false) {
		std::cerr << "ERROR: Problem reading alignment file" << std::endl;
		exit(1);
	}

	isTranslated = false;
	geneticCode = new Code;
}

Alignment::~Alignment(void) {

	delete [] originalMatrix[0];
	delete [] originalMatrix;
	if (isTranslated == true) {
		for (int i = 0; i < numTaxa; i++) {
			for (int j = 0; j < numCodonSites; j++) {
				delete codonMatrix[i][j];
			}
		}
		delete [] codonMatrix[0];
		delete [] codonMatrix;
	}

}

int Alignment::getIndexForTaxon(std::string s) {

	for (int i = 0; i < numTaxa; i++) {
		if (s == taxonNames[i]) {
			return i;
		}
	}
	std::cerr << "WARNING: Could not find taxon \"" << s << "\" in alignment" << std::endl;
	return -1;
}

int Alignment::getNumSenseCodons(void) {

	return geneticCode->getNumSenseCodons();
}

void Alignment::getPossibleNucs(int nucCode, int nuc[]) {
	if (nucCode == 1) {
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 0;
		nuc[3] = 0;
	}
	else if (nucCode == 2) {
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 0;
	}
	else if (nucCode == 3)	{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 0;
	}
	else if (nucCode == 4) {
		nuc[0] = 0;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 0;
	}
	else if (nucCode == 5) {
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 0;
	}
	else if (nucCode == 6) {
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 0;
	}
	else if (nucCode == 7) {
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 0;
	}
	else if (nucCode == 8) {
		nuc[0] = 0;
		nuc[1] = 0;
		nuc[2] = 0;
		nuc[3] = 1;
	}
	else if (nucCode == 9) {
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 0;
		nuc[3] = 1;
	}
	else if (nucCode == 10) {
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 1;
	}
	else if (nucCode == 11) {
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 1;
	}
	else if (nucCode == 12) {
		nuc[0] = 0;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 1;
	}
	else if (nucCode == 13) {
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 1;
	}
	else if (nucCode == 14) {
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 1;
	}
	else if (nucCode == 15) {
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 1;
	}
	else if (nucCode == 16) {
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 1;
	}
	else {
		std::cerr << "ERROR: Code::getPossibleNucs(), unexpected value: " << nucCode << std::endl;
	}
}

/*
int Alignment::getContextID(int t, int pos) {

	char lftPos = ' ';
	if (pos == 0) lftPos = ntSequences[t].at(0);
	else lftPos = ntSequences[t].at(pos - 1);

	char rhtPos = ' ';
	if (pos == numSites) rhtPos = numSites;
	else rhtPos = ntSequences[t].at(pos + 1);

	int contextID = 0;

	if ((lftPos == 'A' || lftPos == 'a') && (rhtPos == 'A' || rhtPos == 'a')) contextID =  0;
	else if ((lftPos == 'A' || lftPos == 'a') && (rhtPos == 'C' || rhtPos == 'c')) contextID =  1;
	else if ((lftPos == 'A' || lftPos == 'a') && (rhtPos == 'G' || rhtPos == 'g')) contextID =  2;
	else if ((lftPos == 'A' || lftPos == 'a') && (rhtPos == 'T' || rhtPos == 't')) contextID =  3;
	else if ((lftPos == 'C' || lftPos == 'c') && (rhtPos == 'A' || rhtPos == 'a')) contextID =  4;
	else if ((lftPos == 'C' || lftPos == 'c') && (rhtPos == 'C' || rhtPos == 'c')) contextID =  5;
	else if ((lftPos == 'C' || lftPos == 'c') && (rhtPos == 'G' || rhtPos == 'g')) contextID =  6;
	else if ((lftPos == 'C' || lftPos == 'c') && (rhtPos == 'T' || rhtPos == 't')) contextID =  7;
	else if ((lftPos == 'G' || lftPos == 'g') && (rhtPos == 'A' || rhtPos == 'a')) contextID =  8;
	else if ((lftPos == 'G' || lftPos == 'g') && (rhtPos == 'C' || rhtPos == 'c')) contextID =  9;
	else if ((lftPos == 'G' || lftPos == 'g') && (rhtPos == 'G' || rhtPos == 'g')) contextID = 10;
	else if ((lftPos == 'G' || lftPos == 'g') && (rhtPos == 'T' || rhtPos == 't')) contextID = 11;
	else if ((lftPos == 'T' || lftPos == 't') && (rhtPos == 'A' || rhtPos == 'a')) contextID = 12;
	else if ((lftPos == 'T' || lftPos == 't') && (rhtPos == 'C' || rhtPos == 'c')) contextID = 13;
	else if ((lftPos == 'T' || lftPos == 't') && (rhtPos == 'G' || rhtPos == 'g')) contextID = 14;
	else if ((lftPos == 'T' || lftPos == 't') && (rhtPos == 'T' || rhtPos == 't')) contextID = 15;
	else std::cerr << "ERROR: Context not found for " << lftPos << "_" << rhtPos << std::endl;

	return contextID;
}
*/

int Alignment::getCodonID(int t, int pos) {

	int codonPos = pos / 3;
	MbBitfield* codonBf = getCodonMatrixEntry(t, codonPos);

	for (int i = 0; i < codonBf->dim(); i++) {
		if (codonBf->isBitSet(i)) return i;
	}

	return -1;
}

int Alignment::nucID(char nuc) {

	char		n;

	if (nuc == 'U' || nuc == 'u')	n = 'T';
	else							n = nuc;
											  // TGCA
	if (n == 'A' || n == 'a')		return 1; // 0001
	else if (n == 'C' || n == 'c')	return 2; // 0010
	else if (n == 'G' || n == 'g')	return 4; // 0100
	else if (n == 'T' || n == 't')	return 8; // 1000
	else if (n == 'R' || n == 'r')	return 5;
	else if (n == 'Y' || n == 'y')	return 10;
	else if (n == 'M' || n == 'm')	return 3;
	else if (n == 'K' || n == 'k')	return 12;
	else if (n == 'S' || n == 's')	return 6;
	else if (n == 'W' || n == 'w')	return 9;
	else if (n == 'H' || n == 'h')	return 11;
	else if (n == 'B' || n == 'b')	return 14;
	else if (n == 'V' || n == 'v')	return 7;
	else if (n == 'D' || n == 'd')	return 13;
	else if (n == 'N' || n == 'n')	return 15;
	else if (n == '-')				return 15;
	else if (n == '?')				return 15;
	else							return -1;
}

void Alignment::print(void) {

	if (isTranslated == false) {
		for (int s = 0; s < numSites; s++) {
			std::cout << std::setw(4) << s << " --";
			for (int t = 0; t < numTaxa; t++) {
				std::cout << std::setw(3) << originalMatrix[t][s];
			}
			std::cout << std::endl;
		}
	}

	else {
		for (int s = 0; s < numCodonSites; s++) {
			std::cout << std::setw(4) << s << " --";
			for (int t = 0; t < numTaxa; t++) {
				for (int i = 0; i < codonMatrix[t][s]->dim(); i++)
				{
					if (codonMatrix[t][s]->isBitSet(i) == true)
						std::cout << "1";
					else
						std::cout << "0";
				}
				std::cout << " ";
			}
			std::cout << std::endl;

		}
	}
}

bool Alignment::readAlignment(FileMgr *file) {

	// open file
	std::ifstream seqStream;
	if (file->openFile(seqStream) == false) {
		std::cerr << "Cannot open file \"" + file->getFileName() + "\"" << std::endl;
		exit(1);
	}

	// intialize variables
	originalMatrix = NULL;
	numTaxa = 0;
	numSites = 0;

	// temp variables
	std::string linestring = "";
	std::string theSequence = "";
	int line = 0;
	int taxonNum = 0;
	bool excludeLine = false;
	bool charSetLine = false;

	// read alignment
	while (getline(seqStream, linestring).good()) {
		std::istringstream linestream(linestring);
		std::string word = "";
		std::string cmdString = "";
		int wordNum = 0;
		int siteNum = 0;
		int ch;
		excludeLine = false;
		charSetLine = false;

		do {
			word = "";
			linestream >> word;
			wordNum++;

			// initialize originalMatrix (untranslated)
			if (line == 0) {
				// read number of taxa/chars from the first line
				std::istringstream buf(word);
				int x;
				buf >> x;
				if (wordNum == 1)
					numTaxa = x;
				else
					numSites = x;

				// initialize originalMatrix according to numTaxa, numSites
				if (numTaxa > 0 && numSites > 0 && originalMatrix == NULL) {
					originalMatrix = new int*[numTaxa];
					originalMatrix[0] = new int[numTaxa * numSites];
					for (int i = 1; i < numTaxa; i++) {
						originalMatrix[i] = originalMatrix[i - 1] + numSites;
					}
					for (int i = 0; i < numTaxa; i++) {
						for (int j = 0; j < numSites; j++) {
							originalMatrix[i][j] = 0;
						}
					}
				}
			}

			// populate originalMatrix
			else {
				if (wordNum == 1) {
					// the first word of the alignment is the taxon name
					taxonNames.push_back(word);
					taxonNum++;
				}
				else {
					// add the sequence as codons (int)nucID to originalMatrix
					for (int i = 0; i < (int)word.length(); i++) {
						char site = word.at(i);
						originalMatrix[taxonNum - 1][siteNum++] = nucID(site);
					}
				}
			}
		} while ((ch = linestream.get()) != EOF);

		// advance to the next sequence
		line++;
	}

	file->closeFile(seqStream);
	return true;
}


bool Alignment::translate(void) {

	if (isTranslated == false) {

		// check the sequence is divisble by 3 (i.e. # nts / codon)
		if (numSites % 3 != 0) {
			std::cerr << "ERROR: Alignment is not evenly divisible by three" << std::endl;
			exit(1);
		}

		// translate
		numCodonSites = numSites / 3;
		codonMatrix = new MbBitfield**[numTaxa];

		// assign indices for codonMatrix based on numTaxa * numCodonSites
		codonMatrix[0] = new MbBitfield*[numTaxa * numCodonSites];
		for (int i = 1; i < numTaxa; i++) {
			codonMatrix[i] = codonMatrix[i - 1] + numCodonSites;
		}

		// assign a codon Bitfield for each codonMatrix element
		for (int i = 0; i < numTaxa; i++) {
			for (int j = 0; j < numCodonSites; j++) {
				codonMatrix[i][j] = new MbBitfield(geneticCode->getNumSenseCodons());
			}
		}

		// translate codon into a Bitfield stored in codonMatrix
		for (int i = 0; i < numTaxa; i++) {
			for (int j = 0; j < numSites; j += 3) {
				MbBitfield *bf = codonMatrix[i][j / 3];
				int nucCode1 = originalMatrix[i][j];
				int nucCode2 = originalMatrix[i][j + 1];
				int nucCode3 = originalMatrix[i][j + 2];
				int nuc1[4], nuc2[4], nuc3[4];
				getPossibleNucs(nucCode1, nuc1);
				getPossibleNucs(nucCode2, nuc2);
				getPossibleNucs(nucCode3, nuc3);

				int s = 0;
				for (int s1 = 0; s1 < 4; s1++) {
					for (int s2 = 0; s2 < 4; s2++) {
						for (int s3 = 0; s3 < 4; s3++) {
							if (nuc1[s1] == 1 && nuc2[s2] == 1 && nuc3[s3] == 1)
								bf->setBit(s);
							if (geneticCode->isStopCodon(s1, s2, s3) == false)
								s++;
						}
					}
				}
			}
		}
	}

	// translation complete
	isTranslated = true;
	return true;
}



