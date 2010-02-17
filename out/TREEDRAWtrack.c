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
{	gROOT->Reset(); // reset ROOT

	TFile *file = TFile::Open(filename, "READ"); // opening the file wherein the tree is stored (read access only)
	TTree *mytree = (TTree*) file->Get("TrackTree"); // creating a pointer 'mytree' pointing on the TTree "mytree"
	
	//======== Drawing and saving z versus r, z versus x and x-y-z  =====================================================
	std::cout << "Drawing ..." << std::endl;
	
	gStyle->SetTitleFillColor(0);

	TCanvas *c1 = new TCanvas("c1", "z:r z:x data from " + filename, 20, 20, 1200, 450); // creating a new
	// TCanvas([canvasname], [canvastitle], x pixel coordinate, y pixel coordinate, x pixle size, y pixel size)
	c1->SetFillColor(0);
	c1->SetBorderMode(0);
	c1->Divide(2,1); // dividing 'c1' into 2*1 pads (numbered like text read)

	c1->cd(1); // select pad 1
	c1_1->SetFrameFillColor(0);
	c1_1->SetFrameBorderMode(0);
	mytree->SetEstimate(mytree->GetEntries()); // setting the estimated lenght of V1, V2 and V3
	mytree->Draw("z:r", "", "goff"); // drawing "z:r" without graphical output
	g11 = new TGraph(mytree->GetEntries(), mytree->GetV2(), mytree->GetV1()); // generating graph and retrieving data
	                                                                          // from the draw command above
	g11->SetTitle("z:r");
	g11->SetMarkerColor(4);
	g11->SetLineColor(4);

	c1->Update(); // necessary command for setting the axis titles
	g11->GetHistogram()->SetXTitle("r [m]");
	g11->GetHistogram()->SetYTitle("z [m]");

	g11->Draw("AP"); //++++++++ options: "A" ~ axis, "P" ~ markers, "L" ~ a simple line between the points ++++++++++++++

	c1->cd(2); // select pad 2
	c1_2->SetFrameFillColor(0);
	c1_2->SetFrameBorderMode(0);
	mytree->Draw("z:x", "", "goff"); // drawing "z:x" without graphical output
	g12 = new TGraph(mytree->GetEntries(), mytree->GetV2(), mytree->GetV1()); // generating graph and retrieving data
	                                                                          // from the draw command above
	g12->SetTitle("z:x");
	g12->SetMarkerColor(4);
	g12->SetLineColor(4);

	c1->Update(); // necessary command for setting the axis titles
	g12->GetHistogram()->SetXTitle("x [m]");
	g12->GetHistogram()->SetYTitle("z [m]");

	g12->Draw("AP"); //++++++++ options: "A" ~ axis, "P" ~ markers, "L" ~ a simple line between the points ++++++++++++++

	c1->SaveAs(filename + "_z-r_z-x.cxx");     // saving canvas as macro
	c1->SaveAs(filename + "_z-r_z-x_0.png");   // saving canvas as PNG
	c1_1->SaveAs(filename + "_z-r_z-x_1.png"); // saving pad 1 as PNG
	c1_2->SaveAs(filename + "_z-r_z-x_2.png"); // saving pad 2 as PNG
/*
	delete c1; // closing canvas window
	delete g11; // deleting graph
	delete g12; // deleting graph
*/
	TCanvas *c2 = new TCanvas("c2", "z:y:x data from " + filename, 620, 200, 600, 450); // creating a new
	// TCanvas([canvasname], [canvastitle], x pixel coordinate, y pixel coordinate, x pixle size, y pixel size)
	c2->SetFillColor(0);
	c2->SetBorderMode(0);

	c2->cd(); // select pad
	c2->SetFrameFillColor(0);
	c2->SetFrameBorderMode(0);
	mytree->Draw("z:y:x", "", "goff"); // drawing "z:y:x" without graphical output
	g2 = new TGraph2D(mytree->GetEntries(), mytree->GetV3(), mytree->GetV2(), mytree->GetV1()); // generating graph and
	//retrieving data from the draw command above

	g2->SetTitle("z:y:x;x [m];y [m];z [m]"); // setting title (including axis titles seperated by ";")
	g2->SetMarkerColor(4);
	g2->SetLineColor(4);
	g2->Draw("P"); //++++++++ options: "P" ~ markers, "LINE" ~ a 3D polyline between the points +++++++++++++++++++++++++

	c2->SaveAs(filename + "_z-y-x.cxx");   // saving canvas as macro
	c2->SaveAs(filename + "_z-y-x_0.png"); // saving canvas as PNG
/*	
	delete c2; // closing canvas window
	delete g2; // deleting graph	
*/
	std::cout << "Drawing done." << std::endl;
	//======== Finished =================================================================================================

//	file->Close(); // closing the file
}
