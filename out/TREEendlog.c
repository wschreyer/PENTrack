//
//	TREEendlog.c is a ROOT scipt that reads and edits a data files first line to refine the branch descriptor in order to
//	construct a ROOT tree named "[filename].root" and finally  restores the first line to its original. For each variable
//	will a branch be created in the tree contianing the data as double type by default.
//
//	User Instructions:
//	(1) Make sure the data files first line contains the names of the variables seperated by '[SPACE]' (e.g. "x y z")!
//	(2) Make sure that the data in the following lines are seperated by '[SPACE]', too (e.g. "0.42 0.42 0.42")!
//	(3) Run the script in ROOT like this: root [0] .x TREEendlog.c("[filename]");
//
//	Note:
//	(1) Any occurrence of errors may possibly lead to a loss of data.
//	(2) If the tree file (or a file with the same name) allready exists, it will by overwriten by default.
//	(3) The ROOT code allows the first line and in consequence the branch descriptor to be of 10000 characters at maximum.
//	(4) Alterable options in the script are marked with at least 8 '+' characters.
//
//	 ___________  __________  __________  __________    ____  ____
//	/___    ___/ /   __    / /   ______/ /   ______/   /___/ /   /
//	   /   /    /   /_/   / /   /____   /   /____     ____  /   _/
//	  /   /    /       __/ /   _____/  /   _____/    /   / /   /
//	 /   /    /   /\   \  /   /_____  /   /_____    /   / /   /
//	/___/    /___/  \___\/_________/ /_________/   /___/ /___/
//

#include <iostream>
#include <fstream>
#include <string>

void TREEendlog(TString filename)
{	gROOT->Reset(); // reset ROOT

	TString rootfilename(filename);
	rootfilename+=".root";
	TFile* file=TFile::Open(rootfilename, "RECREATE"); // recreating a new file wherein the tree will be saved.
	//++++++++ options: "CREATE" ~ create and open a new file if it does not already exist ++++++++++++++++++++++++++++++
	//+++++++++++++++++ "RECREATE" ~ create and overwrite if existing +++++++++++++++++++++++++++++++++++++++++++++++++++
	TTree* tree=new TTree("mytree", "mytree"); // creating a new TTree([treename],[treetitle])

	//======== Refining branch descriptor (~ "x/[type]:y:z") and the files first line (from "x y z" to "#[...]") ========
	std::cout << "Refining ...\n";

	fstream edfile; // new stream 'edfile' to the data file
	edfile.open(filename, ios_base::in | ios_base::out); // read and write access to the data file

	char line1[10000];
	edfile.getline(line1, 10000); // reading out the first line (~ "x y z") in 'line1' an array of chars
	int line1size = edfile.gcount() - 1; // sizes of first line (-1 for the line break '\n')

	bool b = false; // Boolean variable 'b' holding the information if the type is set in the branch descriptor
	std::string bdescriptor (""); // branch descriptor 'bdescriptor' as a empty string
	std::string type ("/D"); // type of all branches (D ~ Double_t)
	//++++++++ options: D ~ Double_t, F ~ Float_t, I ~ Int_t, L ~ Long63_t ++++++++++++++++++++++++++++++++++++++++++++++

	for(int i = 0; i < line1size; i++) // editing descriptor while looping over each character in 'line1'
	{	if(line1[i]==' ')
		{	
			if(!b)
			{	b = true;
				bdescriptor+ = type; // adding the string 'type' by the first call 
			}
			bdescriptor+ = ':'; // adding a ':' instead of a '[SPACE]'
		}
		else
		{	bdescriptor+ = line1[i]; // adding characters unlike '[SPACE]'
		}
	}
	//std::cout << bdescriptor << std::endl; // branch descriptor output (~ "x/[type]:y:z")

	char char1 = line1[0]; // backup the first character in 'char1'
	line1[0] = '#'; // editing first character to '#' to ignore the first line while building the tree

	char* p_line1 = &line1[0];
	edfile.seekp(0); // set position of output pointer
	edfile.write(p_line1, 1); // writing '#' as first charakter into the file

	edfile.close(); // closing the stream 'edfile'

	std::cout << "Refining done.\n";
	//======== Finished =================================================================================================

	TString t_bdescriptor = bdescriptor; // branch descriptor 't_bdescriptor' as TString
	tree->ReadFile(filename.Data(), t_bdescriptor.Data()); // building the tree using the data from the data file and the
	                                                       // self-made branch descriptor 't_bdescriptor'
	tree->Print(); // output of the tree overview
	tree->Write(); // saving the tree as "[rootfilename]" (equal "[filename].root")
	file->Close(); // closing the file

	//======== Restoring the files first line of file (from "x:y:z" to "x y z") =========================================
	std::cout << "Restoring ...\n";

	fstream edfile; // new stream 'edfile' to the data file
	edfile.open(filename, ios_base::in | ios_base::out); // read and write access to the data file

	p_line1 = &char1;
	edfile.seekp(0); // set position of output pointer
	edfile.write(p_line1, 1); // restoring the first character in the first line

	edfile.close(); // closing the stream 'edfile'

	std::cout << "Restoring done.\n";	
	//======== Finished =================================================================================================
}
