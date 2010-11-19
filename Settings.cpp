/*
 * Settings.cpp
 *
 *  Created on: Jan 15, 2010
 *      Author: mikee
 */

#include "Settings.h"

Settings::Settings(int argc, char* argv[]) {

	// set parameters to some default state
	alignmentFileName	= "";
	treeFileName		= "";
	outputFileName		= "";
	coordinatesFileName	= "";
	codonMatrixFileName = "";
	codonFreqFileName  = "";
	freeEnergyFileName  = "";
	settingsFileName	= "settingsmac.ini";
	numOmegaCats = 1;
	brlenLambda = 10.0;
	numCycles = 1000000;
	printFrequency = 1;
	sampleFrequency = 50;
	summarizeFrequency = 1000;
	contactDistance = 5.0;
	siteChangesPerCycle = 20;
	SVESamplesPerCycle = 5;
	useECM = false;

	// initialize method variables
	std::vector<std::string> cmdStringTokens;
	std::vector<std::string>::iterator it;
	std::string cmd = "";			// store individual command tokens
	std::string status = "none";	// temp string for argv arguments
	std::string lineString = "";	// temp string for settingsFile lines
	std::stringstream ss;			// convert to string
	int cmdStringFindPos = 0;	// used to find "=" in settings file

	// initialize cmdStringTokens vector
	for (int i = 1; i < argc; i++) {
		// save -sett_path value if found
		if (strcmp(argv[i],"-sett_path")) {
			settingsFileName = argv[i++];
		}
		// otherwise, store value in cmdStringTokens
	}

	// add settingsFile values to cmdStringTokens vector
	FileMgr settingsFileMgr(settingsFileName);
	std::ifstream settStream;
	if (settingsFileMgr.openFile(settStream) == false) {
		std::cerr << "Cannot open file \"" + settingsFileMgr.getFileName() + "\"" << std::endl;
		exit(1);
	}

	it = cmdStringTokens.begin();
	while(getline(settStream, lineString).good()) {
		std::istringstream linestream(lineString);
		cmdStringFindPos = lineString.find_first_of("=");
		cmd = "-" + lineString.substr(0, cmdStringFindPos); // get arg type
		cmdStringTokens.push_back(cmd);
		cmd = lineString.substr(cmdStringFindPos + 1, lineString.length() - cmdStringFindPos - 1); // get arg value
		cmdStringTokens.push_back(cmd);
	}

	// store argv value in cmdStringTokens vector
	for (int i = 1; i < argc; i++) {
		cmd = argv[i];
		cmdStringTokens.push_back(cmd);
	}

	/*
	// output to verify cmdStringTokens contents
	for (int i = 0; i < (int)cmdStringTokens.size(); i++) {
		std::cout << i << ": " << cmdStringTokens[i] << std::endl;
	}
	*/

	argc = cmdStringTokens.size() - 1;  // -1 adjusts for vector being indexed from 1

	// read values from arguments
	status = "none";
	if (argc > 1) {
		if (argc % 2 == 0) {
			std::cout << "Usage:" << std::endl;
			std::cout << "   -sett_path <PATH>              : Path to file containing settings" << std::endl;
			std::cout << "   -aln_path <PATH>               : Path to file containing the alignment" << std::endl;
			std::cout << "   -tree_path <PATH>              : Path to file containing the newick tree file" << std::endl;
			std::cout << "   -output_path <NAME>            : Output file name" << std::endl;
			std::cout << "   -matrix_path <PATH>            : Path to empirical codon substitution matrix\n";
			std::cout << "   -freq_path <PATH>              : Path to codon frequenceis\n";
			std::cout << "   -coord_path <PATH>             : Path to protein databank file coordinates\n";
			std::cout << "   -energy_path <PATH>            : Path to file containing free energy values\n";
			std::cout << "   -brlen_lambda <NUMBER>         : Exponential parameter for branch length prior" << std::endl;
			std::cout << "   -contact_dist <NUMBER>         : Distance threshold for AA contact mapping\n";
			std::cout << "   -chain_length <NUMBER>         : Number of MCMC cycles" << std::endl;
			std::cout << "   -print_freq <NUMBER>           : Frequency with which output is printed to screen" << std::endl;
			std::cout << "   -sample_freq <NUMBER>          : Frequency with which chain is sampled to a file" << std::endl;
			std::cout << "   -sum_freq <NUMBER>             : Frequency with which chain is summarized" << std::endl;
			std::cout << "   -sites_per_cycle <NUMBER>      : Number of site changes proposed per cycle\n";
			std::cout << "   -SVE_per_cycle <NUMBER>           : Number of SVE samples per cycle\n";
			std::cout << std::endl;
			std::cout << "Example:" << std::endl;
			std::cout << "   mmm -aln_path <input file> -output_name <output file>" << std::endl;
			exit(0);
			}

		/* read the command-line arguments */
		for (int i=0; i<argc; i++) {
			std::string cmd = cmdStringTokens[i];
			if (status == "none") {
				/* read the parameter specifier */
				if ( cmd == "-aln_path" )
					status = "aln_path";
				else if ( cmd == "-sett_path" )
					status = "sett_path";
				else if ( cmd == "-tree_path" )
					status = "tree_path";
				else if ( cmd == "-output_path" )
					status = "output_path";
				else if ( cmd == "-matrix_path" )
					status = "matrix_path";
				else if ( cmd == "-freq_path" )
					status = "freq_path";
				else if ( cmd == "-coord_path" )
					status = "coord_path";
				else if ( cmd == "-energy_path" )
					status = "energy_path";
				else if ( cmd == "-brlen_lambda" )
					status = "brlen_lambda";
				else if ( cmd == "-contact_dist" )
					status = "contact_dist";
				else if ( cmd == "-chain_length" )
					status = "chain_length";
				else if ( cmd == "-print_freq" )
					status = "print_freq";
				else if ( cmd == "-sample_freq" )
					status = "sample_freq";
				else if ( cmd == "-sum_freq" )
					status = "sum_freq";
				else if ( cmd == "-sites_per_cycle" )
					status = "sites_per_cycle";
				else if ( cmd == "-SVE_per_cycle" )
					status = "SVE_per_cycle";
				else {
					std::cerr << "Could not interpret option \"" << cmd << "\"." << std::endl;
					exit(1);
				}
			}
			else {
				/* read the parameter */
				if ( status == "aln_path" ) {
					alignmentFileName = cmdStringTokens[i];
				}
				else if ( status == "tree_path" ) {
					treeFileName = cmdStringTokens[i];
				}
				else if ( status == "output_path" ) {
					outputFileName = cmdStringTokens[i];
				}
				else if ( status == "matrix_path" ) {
					codonMatrixFileName = cmdStringTokens[i];
					useECM = true;
				}
				else if ( status == "freq_path" ) {
					codonFreqFileName = cmdStringTokens[i];
					useECM = true;
				}
				else if ( status == "coord_path" ) {
					coordinatesFileName = cmdStringTokens[i];
				}
				else if ( status == "energy_path" ) {
					freeEnergyFileName = cmdStringTokens[i];
				}
				else if ( status == "brlen_lambda" ) {
					ss << cmdStringTokens[i];
					brlenLambda = atof(cmdStringTokens[i].c_str());
				}
				else if ( status == "contact_dist" ) {
					ss << cmdStringTokens[i];
					contactDistance = atof(cmdStringTokens[i].c_str());
				}
				else if ( status == "chain_length" ) {
					ss << cmdStringTokens[i];
				//	ss >> numCycles;
					numCycles = atof(cmdStringTokens[i].c_str());
				}
				else if ( status == "print_freq" ) {
					ss << cmdStringTokens[i];
				//	ss >> printFrequency;
					printFrequency = atof(cmdStringTokens[i].c_str());
				}
				else if ( status == "sample_freq" )	{
					ss << cmdStringTokens[i];
				//	ss >> sampleFrequency;
					sampleFrequency = atof(cmdStringTokens[i].c_str());
				}
				else if ( status == "sum_freq" ) {
					ss << cmdStringTokens[i];
					summarizeFrequency = atoi(cmdStringTokens[i].c_str());
				}
				else if ( status == "sites_per_cycle" ) {
					ss << cmdStringTokens[i];
					siteChangesPerCycle = atoi(cmdStringTokens[i].c_str());
				}
				else if ( status == "SVE_per_cycle" ) {
					ss << cmdStringTokens[i];
					SVESamplesPerCycle = atoi(cmdStringTokens[i].c_str());
				}
				else {
					std::cerr << "Unknown status reading command line information" << std::endl;
					exit(1);
				}
				status = "none";
			}
		}

		std::cout << "Model settings:" << std::endl;
		std::cout << "   settingsFileName    = " << settingsFileName << std::endl;
		std::cout << "   alignmentFileName   = " << alignmentFileName << std::endl;
		std::cout << "   treeFileName        = " << treeFileName << std::endl;
		std::cout << "   outputFileName      = " << outputFileName << std::endl;
		std::cout << "   codonMatrixFileName = " << codonMatrixFileName << "\n";
		std::cout << "   codonFreqFileName   = " << codonFreqFileName << "\n";
		std::cout << "   coordinatesFileName = " << coordinatesFileName << "\n";
		std::cout << "   freeEnergyFileName  = " << freeEnergyFileName << "\n";
		std::cout << "   brlenLambda         = " << brlenLambda << std::endl;
		std::cout << "   contactDistance     = " << contactDistance << "\n";
		std::cout << "   numCycles           = " << numCycles  << std::endl;
		std::cout << "   printFrequency      = " << printFrequency << std::endl;
		std::cout << "   sampleFrequency     = " << sampleFrequency << std::endl;
		std::cout << "   summarizeFrequency  = " << summarizeFrequency << std::endl;
		std::cout << "   siteChangesPerCycle = " << siteChangesPerCycle << "\n";
		std::cout << "   SVESamplesPerCycle  = " << SVESamplesPerCycle << "\n";
	}
	else {
		std::cerr << "ERROR: You must provide command line arguments" << std::endl;
	}
}
