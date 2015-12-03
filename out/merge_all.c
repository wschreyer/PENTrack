#include <map>
#include <iostream>
#include <fstream>

#include "TSystemDirectory.h"
#include "TFile.h"
#include "TNtupleD.h"
#include "TRegexp.h"

using namespace std;

/**
 * Merge all files in folder path with names matching the regular expression filematch into trees in a ROOT file.
 *
 * @param filematch Merge files whose names match this regular expression (default: read all)
 * @param path Merge files in this folder (default: current folder)
 */
void merge_all(const char *filematch = "[0-9]+[A-Za-z]+.out", const char *path = ".")
{
	TFile *outfile = new TFile("out.root","RECREATE"); // create ROOT file
	ifstream infile;
	TSystemDirectory dir(path,path); // open given folder
	TList *files = dir.GetListOfFiles(); // get all files in that folder
	TRegexp re(filematch); // create regular expression from given parameter

	outfile->cd(); // switch current directory to ROOT file
	
	Int_t n = 0;
	for (Int_t i = 0; ; i++){ // for loop incrementing index i
		TObject *f = files->At(i); // get file from folder with index i
		if (f){ // if next file was found
			TString filename = f->GetName(); // get filename
			TString fn = filename;
			fn.Resize(filename.Length() - 4); // shorten filename by extension ".out"
			
			ULong64_t jobnumber;
			char logtype[64];
			double data[1024];
			
			if ((re.Index(filename, &n) == 0) && (sscanf(fn.Data(), "%Ld%s", &jobnumber, logtype) == 2)){ // if filename matches regular expression and contains jobnumber
				infile.open(filename.Data()); // open file
				TNtupleD *tree = (TNtupleD*)outfile->Get(logtype); // get corresponding tree from file
				
				if (!tree) { // if tree does not yet exist
					TString bdescriptor;
					bdescriptor.ReadLine(infile); // read branch descriptor from file header
					bdescriptor.ReplaceAll(" ",":"); // format branch descriptor for root ("x y z" -> "x:y:z")

					tree = new TNtupleD(logtype, logtype, bdescriptor.Data()); // create new tree with name logtype from filename 
					printf("%ss have %i columns\n", logtype, tree->GetNvar());
				}
				else
					infile.ignore(9999, '\n'); // if tree already exists skip file header

				n = tree->GetNvar(); // get number of file columns
				cout << "Adding " << filename << '\n';
				while (1){
					for (int j = 0; j < n; j++) infile >> data[j]; // read values into data array
					if (!infile) break; // if something happened during reading: stop
					tree->Fill(data); // fill data into tree
				}
				infile.close(); // close file
				cout << "Entries: " << tree->GetEntries() << '\n';
			}
		}
		else break;
	}
	cout << "Done!" << endl;
	
	outfile->Write();
	outfile->Close(); // close ROOT file
	delete outfile;
}
