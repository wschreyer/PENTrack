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

void TREEDRAWdecaytrack(TString filename)
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
	
	TLine *lineP = new TLine();   // line to depict prot det
	lineP->SetLineWidth(3);
	lineP->SetLineColor(kGray+2);

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
	
	Double_t xmax = +0.7, ymax = +0.7, zmax = +1.2; // dummy ~ range maximum in 'g2', 'g2n', 'g2p' and 'g2e'
	Double_t xmin = -0.7, ymin = -0.7, zmin = -0.5; // dummy ~ range minuium in 'g2', 'g2n', 'g2p' and 'g2e'
	
	gStyle->SetTitleFillColor(0);

	
	
	
	TCanvas *c3 = new TCanvas("c3", "r:z data from " + rootfilename, 20, 20, 600, 800); // creating a new
	// TCanvas([canvasname], [canvastitle], x pixel coordinate, y pixel coordinate, x pixle size, y pixel size)
	c3->SetFillColor(0);
	c3->SetBorderMode(0);

	c3->cd(); // select pad
	c3->SetFrameFillColor(0);
	c3->SetFrameBorderMode(0);
	
	mytree->SetEstimate(mytree->GetEntries()); // setting the estimated lenght of V1, V2 and V3
	mytree->Draw("z:r", "", "goff"); // drawing "z:r" without graphical output
	g11 = new TGraph(mytree->GetEntries(), mytree->GetV2(), mytree->GetV1()); // generating graph and retrieving data
	                                                                          // from the draw command above
	g11->SetTitle("z:r");
	g11->SetMarkerColor(0);
	g11->SetLineColor(0);

	c3->Update(); // necessary command for setting the axis titles
	g11->GetHistogram()->SetXTitle("r [m]");
	g11->GetHistogram()->SetYTitle("z [m]");

	g11->Draw("AP"); //++++++++ options: "A" ~ axis, "P" ~ markers, "L" ~ a simple line between the points ++++++++++++++

	//-------- Generating the graphs for neutrons, protons and electrons ------------------------------------------------
	{	mytree->Draw("z:r", "protneut==1", "goff"); // drawing "z:x" (if "protneut==1") without graphical output
		lineP->DrawLine(0.18,1.15,0.36,1.15);	
		
		g11n = new TGraph(mytree->GetEntries("protneut==1"), mytree->GetV2(), mytree->GetV1());
		// generating graph and retrieving data from the draw command above
		g11n->SetMarkerColor(kGreen+1);
		g11n->SetMarkerStyle(6);
		g11n->SetLineColor(kGreen+1);

		mytree->Draw("z:r", "protneut==2", "goff"); // drawing "z:x" (if "protneut==2") without graphical output
		g11p = new TGraph(mytree->GetEntries("protneut==2"), mytree->GetV2(), mytree->GetV1());
		// generating graph and retrieving data from the draw command above
		g11p->SetMarkerColor(kBlue+1);
		g11p->SetMarkerStyle(6);
		g11p->SetLineColor(kBlue+1);
		
		mytree->Draw("z:r", "protneut==6", "goff"); // drawing "z:x" (if "protneut==6") without graphical output
		g11e = new TGraph(mytree->GetEntries("protneut==6"), mytree->GetV2(), mytree->GetV1());
		// generating graph and retrieving data from the draw command above
		g11e->SetMarkerColor(kRed+1);
		g11e->SetMarkerStyle(6);
		g11e->SetLineColor(kRed+1);
		
		
	}
	//-------- Finished -------------------------------------------------------------------------------------------------

	
	g11n->Draw("PSAME");
	g11p->Draw("PSAME");
	g11e->Draw("PSAME");
	
	TCanvas *c1 = new TCanvas("c1", "z:r z:x data from " + rootfilename, 20, 20, 1200, 900); // creating a new
	// TCanvas([canvasname], [canvastitle], x pixel coordinate, y pixel coordinate, x pixle size, y pixel size)
	c1->SetFillColor(0);
	c1->SetBorderMode(0);
	c1->Divide(2,2); // dividing 'c1' into 2*2 pads (numbered like text read)

	c1->cd(1); // select pad 1
	c1_1->SetFrameFillColor(0);
	c1_1->SetFrameBorderMode(0);
	

	c1->cd(2); // select pad 2
	c1_2->SetFrameFillColor(0);
	c1_2->SetFrameBorderMode(0);

	mytree->Draw("z:x", "", "goff"); // drawing "z:x" without graphical output
	g12 = new TGraph(mytree->GetEntries(), mytree->GetV2(), mytree->GetV1()); // generating graph and retrieving data
	                                                                          // from the draw command above

	g12->SetTitle("z:x");
	g12->SetMarkerColor(0);
	g12->SetLineColor(0);

	c1->Update(); // necessary command for setting the axis titles
	g12->GetHistogram()->SetXTitle("x [m]");
	g12->GetHistogram()->SetYTitle("z [m]");

	g12->Draw("AP"); //++++++++ options: "A" ~ axis, "P" ~ markers, "L" ~ a simple line between the points ++++++++++++++

	//-------- Generating the graphs for neutrons, protons and electrons ------------------------------------------------
	{	mytree->Draw("z:x", "protneut==1", "goff"); // drawing "z:x" (if "protneut==1") without graphical output
		g12n = new TGraph(mytree->GetEntries("protneut==1"), mytree->GetV2(), mytree->GetV1());
		// generating graph and retrieving data from the draw command above
		g12n->SetMarkerColor(kGreen+1);
		g12n->SetLineColor(kGreen+1);

		mytree->Draw("z:x", "protneut==2", "goff"); // drawing "z:x" (if "protneut==2") without graphical output
		g12p = new TGraph(mytree->GetEntries("protneut==2"), mytree->GetV2(), mytree->GetV1());
		// generating graph and retrieving data from the draw command above
		g12p->SetMarkerColor(kBlue+1);
		g12p->SetLineColor(kBlue+1);

		mytree->Draw("z:x", "protneut==6", "goff"); // drawing "z:x" (if "protneut==6") without graphical output
		g12e = new TGraph(mytree->GetEntries("protneut==6"), mytree->GetV2(), mytree->GetV1());
		// generating graph and retrieving data from the draw command above
		g12e->SetMarkerColor(kRed+1);
		g12e->SetLineColor(kRed+1);
	}
	//-------- Finished -------------------------------------------------------------------------------------------------

	g12n->Draw("PSAME");
	g12p->Draw("PSAME");
	g12e->Draw("PSAME");
//
	c1->cd(4); // select pad 2
	c1_4->SetFrameFillColor(0);
	c1_4->SetFrameBorderMode(0);

	mytree->Draw("y:x", "", "goff"); // drawing "y:x" without graphical output
	g14 = new TGraph(mytree->GetEntries(), mytree->GetV2(), mytree->GetV1()); // generating graph and retrieving data
	                                                                          // from the draw command above

	g14->SetTitle("y:x");
	g14->SetMarkerColor(0);
	g14->SetLineColor(0);

	c1->Update(); // necessary command for setting the axis titles
	g14->GetHistogram()->SetXTitle("x [m]");
	g14->GetHistogram()->SetYTitle("y [m]");

	g14->Draw("AP"); //++++++++ options: "A" ~ axis, "P" ~ markers, "L" ~ a simple line between the points ++++++++++++++

	//-------- Generating the graphs for neutrons, protons and electrons ------------------------------------------------
	{	mytree->Draw("y:x", "protneut==1", "goff"); // drawing "z:x" (if "protneut==1") without graphical output
		g14n = new TGraph(mytree->GetEntries("protneut==1"), mytree->GetV2(), mytree->GetV1());
		// generating graph and retrieving data from the draw command above
		g14n->SetMarkerColor(kGreen+1);
		g14n->SetLineColor(kGreen+1);

		mytree->Draw("y:x", "protneut==2", "goff"); // drawing "z:x" (if "protneut==2") without graphical output
		g14p = new TGraph(mytree->GetEntries("protneut==2"), mytree->GetV2(), mytree->GetV1());
		// generating graph and retrieving data from the draw command above
		g14p->SetMarkerColor(kBlue+1);
		g14p->SetLineColor(kBlue+1);

		mytree->Draw("y:x", "protneut==6", "goff"); // drawing "z:x" (if "protneut==6") without graphical output
		g14e = new TGraph(mytree->GetEntries("protneut==6"), mytree->GetV2(), mytree->GetV1());
		// generating graph and retrieving data from the draw command above
		g14e->SetMarkerColor(kRed+1);
		g14e->SetLineColor(kRed+1);
	}
	//-------- Finished -------------------------------------------------------------------------------------------------

	g14n->Draw("PSAME");
	g14p->Draw("PSAME");
	g14e->Draw("PSAME");
//
	c1->SaveAs(rootfilename + "_z-r_z-x.cxx");     // saving canvas as macro
	c1->SaveAs(rootfilename + "_z-r_z-x_0.png");   // saving canvas as PNG
	c1_1->SaveAs(rootfilename + "_z-r_z-x_1.png"); // saving pad 1 as PNG
	c1_2->SaveAs(rootfilename + "_z-r_z-x_2.png"); // saving pad 2 as PNG
/*
	delete c1; // closing canvas window
	delete g11;  // deleting graph
	delete g11n; // deleting graph
	delete g11p; // deleting graph
	delete g11e; // deleting graph
	delete g12;  // deleting graph
	delete g12n; // deleting graph
	delete g12p; // deleting graph
	delete g12e; // deleting graph
	delete g14;  // deleting graph
	delete g14n; // deleting graph
	delete g14p; // deleting graph
	delete g14e; // deleting graph
*/
	TCanvas *c2 = new TCanvas("c2", "z:y:x data from " + rootfilename, 20, 20, 600, 800); // creating a new
	// TCanvas([canvasname], [canvastitle], x pixel coordinate, y pixel coordinate, x pixle size, y pixel size)
	c2->SetFillColor(0);
	c2->SetBorderMode(0);

	c2->cd(); // select pad
	c2->SetFrameFillColor(0);
	c2->SetFrameBorderMode(0);

	mytree->Draw("z:y:x", "", "goff"); // drawing "z:y:x" without graphical output
	g2 = new TGraph2D(mytree->GetEntries(), mytree->GetV3(), mytree->GetV2(), mytree->GetV1()); // generating graph and
	//retrieving data from the draw command above
	g2->SetPoint((mytree->GetEntries() + 1), xmax, ymax, zmax);

	g2->SetTitle("z:y:x;x [m];y [m];z [m]"); // setting title (including axis titles seperated by ";")
	g2->SetMarkerColor(0);
	g2->SetLineColor(0);
	g2->Draw("P"); //++++++++ options: "P" ~ markers, "LINE" ~ a 3D polyline between the points +++++++++++++++++++++++++

	//-------- Generating the graphs for neutrons, protons and electrons ------------------------------------------------
	{	mytree->Draw("z:y:x", "protneut==1", "goff"); // drawing "z:y:x" (if "protneut==1") without graphical output
		g2n = new TGraph2D(mytree->GetEntries("protneut==1"), mytree->GetV3(), mytree->GetV2(), mytree->GetV1());
		// generating graph and retrieving data from the draw command above
		g2n->SetPoint((mytree->GetEntries("protneut==1") + 1), xmin, ymin, zmin);
		g2n->SetPoint((mytree->GetEntries("protneut==1") + 2), xmax, ymax, zmax);
		g2n->SetMarkerColor(kGreen+1);
		g2n->SetLineColor(kGreen+1);

		mytree->Draw("z:y:x", "protneut==2", "goff"); // drawing "z:y:x" (if "protneut==2") without graphical output
		g2p = new TGraph2D(mytree->GetEntries("protneut==2"), mytree->GetV3(), mytree->GetV2(), mytree->GetV1());
		// generating graph and retrieving data from the draw command above
		g2p->SetPoint((mytree->GetEntries("protneut==2") + 1), xmin, ymin, zmin);
		g2p->SetPoint((mytree->GetEntries("protneut==2") + 2), xmax, ymax, zmax);
		g2p->SetMarkerColor(kBlue+1);
		g2p->SetLineColor(kBlue+1);

		mytree->Draw("z:y:x", "protneut==6", "goff"); // drawing "z:y:x" (if "protneut==6") without graphical output
		g2e = new TGraph2D(mytree->GetEntries("protneut==6"), mytree->GetV3(), mytree->GetV2(), mytree->GetV1());
		// generating graph and retrieving data from the draw command above
		g2e->SetPoint((mytree->GetEntries("protneut==6") + 1), xmin, ymin, zmin);
		g2e->SetPoint((mytree->GetEntries("protneut==6") + 2), xmax, ymax, zmax);
		g2e->SetMarkerColor(kRed+1);
		g2e->SetLineColor(kRed+1);
	}
	//-------- Finished -------------------------------------------------------------------------------------------------

	g2n->Draw("PSAME");
	g2p->Draw("PSAME");
	g2e->Draw("PSAME");

	c2->SaveAs(rootfilename + "_z-y-x.cxx");   // saving canvas as macro
	c2->SaveAs(rootfilename + "_z-y-x_0.png"); // saving canvas as PNG
/*	
	delete c2; // closing canvas window
	delete g2;  // deleting graph
	delete g2n; // deleting graph
	delete g2p; // deleting graph
	delete g2e; // deleting graph
*/
	std::cout << "Drawing done." << std::endl;
	//======== Finished =================================================================================================

//	file->Close(); // closing the file
}

