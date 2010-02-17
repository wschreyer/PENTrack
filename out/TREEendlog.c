//
//	TREEendlog.c is a ROOT (5.23/01) scipt that reads and edits a data files first line to refine the branch descriptor
//	in order to construct a ROOT tree named "[filename].root" and finally  restores the first line to its original. For
//	each variable will a branch be created in the tree contianing the data as double type by default.
//
//	User Instructions:
//	(1) Make sure the data files first line contains the names of the variables seperated by '[SPACE]' (e.g. "x y z")!
//	(2) Make sure that the data in the following lines are seperated by '[SPACE]', too (e.g. "0.42 0.42 0.42")!
//	(3) Make sure that the data are compatible to the type D / Double_t (a 64 bit floating point)!
//	(4) Run the script in ROOT like this: root [0] .x TREEendlog.c([jobnumber]);
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

void TREEendlog(Int_t jobnumber)
{	gROOT->Reset(); // reset ROOT
	TString filename;
	
	filename.Form("%08dend.out.root", jobnumber);
	cout << "Creating tree " << filename << endl;
	TFile* file=TFile::Open(filename.Data(), "RECREATE"); // recreating a new file wherein the tree will be saved.
	//++++++++ options: "CREATE" ~ create and open a new file if it does not already exist ++++++++++++++++++++++++++++++
	//+++++++++++++++++ "RECREATE" ~ create and overwrite if existing +++++++++++++++++++++++++++++++++++++++++++++++++++
	ifstream edfile; // new stream 'edfile' to the data file

	filename.Form("%08dend.out", jobnumber);
	edfile.open(filename, ios_base::in); // read and write access to the data file	
	if (edfile.good()){
		cout << "Reading data from " << filename << endl;
		TString bdescriptor; // branch descriptor 'bdescriptor' as a empty string
		bdescriptor.ReadLine(edfile); // read branch descriptor from file
		bdescriptor.ReplaceAll(" ",":"); // format branch descriptor for root ("x y z" -> "x:y:z")

		TNtupleD* tree=new TNtupleD("EndTree", "end log", bdescriptor.Data()); // creating a new TTree([treename],[treetitle])	
		Int_t n = tree->GetNvar();
		while (1){
			for (int i = 0; i < n; i++) edfile >> tree->GetArgs()[i];
			if (!edfile.good()) break;
			tree->Fill(tree->GetArgs());
		}
		
		edfile.close(); // closing the stream 'edfile'
		file->Write();
		tree->Print(); // output of the tree overview
		delete tree;
	}
	
	
	filename.Form("%08dtrack001.out",jobnumber);
	edfile.open(filename, ios_base::in); // read and write access to the data file	
	if (edfile.good()){
		cout << "Reading data from " << filename << endl;
		TString bdescriptor; // branch descriptor 'bdescriptor' as a empty string
		bdescriptor.ReadLine(edfile); // read branch descriptor from file
		bdescriptor.ReplaceAll(" ",":"); // format branch descriptor for root ("x y z" -> "x:y:z")

		TNtupleD* tree=new TNtupleD("TrackTree", "track log", bdescriptor.Data()); // creating a new TTree([treename],[treetitle])	
		Int_t n = tree->GetNvar();
		while (edfile.good()){
			for (int i = 0; i < n; i++) edfile >> tree->GetArgs()[i];
			if (!edfile.good()) break;
			tree->Fill(tree->GetArgs());
		}		
		edfile.close(); // closing the stream 'edfile'
		
		int f = 1;
		while (1){
			f++;
			filename.Form("%08dtrack%03d.out",jobnumber,f);
			edfile.open(filename, ios_base::in); // read and write access to the data file	
			if (edfile.good()){
				cout << "Reading data from " << filename << endl;
				bdescriptor.ReadLine(edfile);
				while (edfile.good()){
					for (int i = 0; i < n; i++) edfile >> tree->GetArgs()[i];
					if (!edfile.good()) break;
					tree->Fill(tree->GetArgs());
				}				
				edfile.close(); // closing the stream 'edfile'
			}
			else break;
		}
		file->Write();
		tree->Print(); // output of the tree overview
		delete tree;
	}
	
	
	
	file->Close(); // closing the file	
}
