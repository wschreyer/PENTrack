//
//	TREEDRAWtrack.c is a ROOT (5.23/01) scipt that reads and edits a data files first line to refine the branch
//	descriptor in order to construct a ROOT tree named "[filename].root", restores the first line to its original and
//	finally three plots are drawn and saved to hard disc as a macros and PNGs. For each variable will a branch be created
//	in the tree contianing the data as double type by default. The Plots will show z versus r, z versus x and x-y-z.
//
//	User Instructions:
//	(1) Make sure the data files first line contains the names of the variables seperated by '[SPACE]' (e.g. "x y z")!
//	(2) Make sure that the data in the following lines are seperated by '[SPACE]', too (e.g. "0.42 0.42 0.42")!
//	(3) Make sure that the data are compatible to the type D / Double_t (a 64 bit floating point)!
//	(3) Make sure the data files first line contains at least "x", "y", "z" and "r", otherwise the script will fail to
//	    draw the plots!
//	(4) Run the script in ROOT like this: root [0] .x TREEDRAWtrack.c("[filename]");
//
//	Note:
//	(1) Any occurrence of errors may possibly lead to a loss of data.
//	(2) If the tree file (or a file with the same name) allready exists, it will by overwriten by default.
//	(3) The ROOT code allows the first line and in consequence the branch descriptor to be of 10000 characters at maximum.
//	(4) The created macros are named "[filename].root_[plotname].cxx" (e.g. "[filename].root_z-y-x.cxx").
//	(5) The created PNGs are named "[filename].root_[plotname]_[counter].png" (e.g. "[filename].root_z-y-x_0.png").
//	(5) If the macro or PNG files (or files with the same name) allready exist, they will by overwriten.
//	(6) Alterable options in the script are marked with at least 8 '+' characters.
//
//	 ___________  __________  __________  __________  ________    __________  ________    ____      ____    ____  ____
//	/___    ___/ /   __    / /   ______/ /   ______/ /   __   \  /   __    / /   __   \  /   /     /   /   /___/ /   /
//	   /   /    /   /_/   / /   /____   /   /____   /   / /   / /   /_/   / /   /_/   / /   / __  /   /   ____  /   _/
//	  /   /    /       __/ /   _____/  /   _____/  /   / /   / /       __/ /   __    / /   /_/ /_/   /   /   / /   /
//	 /   /    /   /\   \  /   /_____  /   /_____  /   /_/   / /   /\   \  /   / /   / /             /   /   / /   /
//	/___/    /___/  \___\/_________/ /_________/ /_________/ /___/  \___\/___/ /___/  \____________/   /___/ /___/
//

#include <iostream>
#include <fstream>
#include <string>

void TREEDRAWtrack(TString filename)
{	gROOT->Reset();

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

	edfile.close(); // closing stream 'edfile'

	std::cout << "Refining done.\n";
	//======== Finished =================================================================================================

	TString t_bdescriptor = bdescriptor; // branch descriptor 't_bdescriptor' as TString
	tree->ReadFile(filename.Data(), t_bdescriptor.Data()); // building the tree using the data from the data file and the
	                                                       // self-made branch descriptor 't_bdescriptor'	
	tree->Print(); // output of the tree overview
	tree->Write(); // saving the tree as "[rootfilename]" (equal "[filename].root")

	//======== Restoring the files first line of file (from "x:y:z" to "x y z") =========================================
	std::cout << "Restoring ...\n";

	fstream edfile; // new (editing) stream 'edfile' to the data file
	edfile.open(filename, ios_base::in | ios_base::out); // read and write access to the data file

	p_line1 = &char1;
	edfile.seekp(0); // set position of output pointer
	edfile.write(p_line1, 1); // restoring the first character in the first line

	edfile.close(); // closing stream 'edfile'

	std::cout << "Restoring done.\n";	
	//======== Finished =================================================================================================

	//======== Drawing and saving z versus r, z versus x and x-y-z  =====================================================
	std::cout << "Drawing ..." << std::endl;

	TCanvas *c1 = new TCanvas("c1", "z:r z:x data from " + rootfilename, 20, 20, 1200, 450); // creating a new
	// TCanvas([canvasname], [canvastitle], x pixel coordinate, y pixel coordinate, x pixle size, y pixel size)
	c1->Divide(2,1); // dividing 'c1' into 2*1 pads (numbered like text read)

	c1->cd(1); // select pad 1
	mytree->SetEstimate(mytree->GetEntries()); // setting the estimated lenght of V1, V2 and V3
	mytree->Draw("z:r", "", "goff"); // drawing "z:r" without graphical output
	g11 = new TGraph(mytree->GetEntries(), mytree->GetV2(), mytree->GetV1()); // generating graph and retrieving data
	                                                                          // from the draw command above
	g11->SetTitle("z:r");
	g11->SetMarkerColor(4);
	g11->SetLineColor(4);
	g11->Draw("AP"); //++++++++ options: "A" ~ axis, "P" ~ markers, "L" ~ a simple line between the points ++++++++++++++

	c1->Update(); // necessary command for setting the axis titles
	g11->GetHistogram()->SetXTitle("r [m]");
	g11->GetHistogram()->SetYTitle("z [m]");

	c1->cd(2); // select pad 2
	mytree->Draw("z:x", "", "goff"); // drawing "z:x" without graphical output
	g12 = new TGraph(mytree->GetEntries(), mytree->GetV2(), mytree->GetV1()); // generating graph and retrieving data
	                                                                          // from the draw command above
	g12->SetTitle("z:x");
	g12->SetMarkerColor(4);
	g12->SetLineColor(4);
	g12->Draw("AP"); //++++++++ options: "A" ~ axis, "P" ~ markers, "L" ~ a simple line between the points ++++++++++++++

	c1->Update(); // necessary command for setting the axis titles
	g12->GetHistogram()->SetXTitle("x [m]");
	g12->GetHistogram()->SetYTitle("z [m]");

	c1->SaveAs(rootfilename + "_z-r_z-x.cxx");     // saving canvas as macro
	c1->SaveAs(rootfilename + "_z-r_z-x_0.png");   // saving canvas as PNG
	c1_1->SaveAs(rootfilename + "_z-r_z-x_1.png"); // saving pad 1 as PNG
	c1_2->SaveAs(rootfilename + "_z-r_z-x_2.png"); // saving pad 2 as PNG
/*
	delete c1; // closing canvas window
	delete g11; // deleting graph
	delete g12; // deleting graph
*/
	TCanvas *c2 = new TCanvas("c2", "z:y:x data from " + rootfilename, 620, 200, 600, 450); // creating a new
	// TCanvas([canvasname], [canvastitle], x pixel coordinate, y pixel coordinate, x pixle size, y pixel size)

	c2->cd(); // select pad
	mytree->Draw("z:y:x", "", "goff"); // drawing "z:y:x" without graphical output
	g2 = new TGraph2D(mytree->GetEntries(), mytree->GetV3(), mytree->GetV2(), mytree->GetV1()); // generating graph and
	//retrieving data from the draw command above

	g2->SetTitle("z:y:x;x [m];y [m];z [m]"); // setting title (including axis titles seperated by ";")
	g2->SetMarkerColor(4);
	g2->SetLineColor(4);
	g2->Draw("P"); //++++++++ options: "P" ~ markers, "LINE" ~ a 3D polyline between the points +++++++++++++++++++++++++

	c2->SaveAs(rootfilename + "_z-y-x.cxx");   // saving canvas as macro
	c2->SaveAs(rootfilename + "_z-y-x_0.png"); // saving canvas as PNG
/*	
	delete c2; // closing canvas window
	delete g2; // deleting graph	
*/
	std::cout << "Drawing done." << std::endl;
	//======== Finished =================================================================================================

//	file->Close(); // closing the file
}
