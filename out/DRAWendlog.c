//
//	DRAWendlog.c is a ROOT (5.23/01) scipt that draws a few plots from a ROOT tree named "[filename].root" and saves them
//	as macros and PNGs. It is highly recommended to run TREEendlog.c before with the same data file.
//	The followong Plots will be created:
//	(1) 'kennz' versus 'zstart' versus 'rstart'
//	(2) histogram of 't'
//	(3) histogram of 'H'
//	(4) histogram of 'rend'
//	(5) histogram of 'gammaend'
//	(6) histogram of 'NeutEnergie'
//	(7) histogram of 'NeutEnergie'
//	(8) 'NeutEnergie' versus 't'
//	(9) histograms of 'EFeldSkal' and 'BFeldSkal'
//	(10) histogram of 'kennz'
//	(11) 'zstart' versus 'rstart'
// (12) 'kennz' versus 'zend' versus 'rend'
//
//	User Instructions:
//	(1) Make sure the ROOT tree exists and is named "mytree"!
//	(2) Make sure the ROOT tree contains at least the branches named in the list of plots above!
//	(3) Make sure the branches contain leaves compatible to the type D / Double_t (a 64 bit floating point)!
//	(4) Run the script in ROOT like this: root [0] .x DRAWendlog.c("[filename]");
//
//	Note:
//	(1) Any occurrence of errors may possibly lead to a loss of data.
//	(2) The created macros are named "[filename].root_[plotname].cxx" (e.g. "[filename].root_kennz.cxx").
//	(3) The created PNGs are named "[filename].root_[plotname]_[counter].cxx" (e.g. "[filename].root_kennz.png").
//	(4) If the macro or PNG files (or files with the same name) allready exist, they will by overwriten.
//	(5) Alterable options in the script are marked with at least 8 '+' characters.
//	(6) At the beginning of each draw routine you can find the assignment of a few dummies, which offer an easy way to
//	    modify the routine without altering the hole code.
//
//	     ________    __________  ________    ____      ____    ____  ____
//	    /   __   \  /   __    / /   __   \  /   /     /   /   /___/ /   /
//	   /   / /   / /   /_/   / /   /_/   / /   / __  /   /   ____  /   _/
//	  /   / /   / /       __/ /   __    / /   /_/ /_/   /   /   / /   /
//	 /   /_/   / /   /\   \  /   / /   / /             /   /   / /   /
//	/_________/ /___/  \___\/___/ /___/  \____________/   /___/ /___/
//

#include <iostream>
#include <cmath>
#include <string>

void DRAWendlog(TString filename)
{	gROOT->Reset();

	TString vnamex; // dummy
	TString vnamey; // dummy
	TString vnamec; // dummy
	Int_t nbins;    // dummy
	Double_t xmin;  // dummy
	Double_t xmax;  // dummy

	TString rootfilename(filename);
	rootfilename+ = ".root";
	TFile *file = TFile::Open(rootfilename, "READ"); // opening the file wherein the tree is stored (read access only)
	TTree *mytree = (TTree*) file->Get("mytree"); // creating a pointer 'mytree' pointing on the TTree "mytree"
	//file->cd();
	
	Int_t particle = 0; //++++++++ options: 0 ~ unspecified, 1 ~ neutron, 2 ~ protron, 6 ~ electron +++++++++++++++++++++
	//-------- Generating particle type selection criterion 'pnamec' ----------------------------------------------------
	TString pnamec = ""; // particle type selection criterion
	{	if(particle == 0) // choose predominant particle type or 1 (~ neutron)
		{	particle = 1;
			if(mytree->GetEntries("protneut==6") > mytree->GetEntries("protneut==1")) particle = 6;
			if(mytree->GetEntries("protneut==2") > mytree->GetEntries("protneut==1")) particle = 2;
		}	
		pnamec+ = "(protneut==";
		pnamec+ = particle;
		pnamec+ = ")";
	}
	//-------- Finished -------------------------------------------------------------------------------------------------

	const Long_t nentries = mytree->GetEntries(pnamec); // total number of entries in 'mytree'
	
	const Int_t ncolors = gROOT->GetListOfColors()->LastIndex(); // number of the last idicated colour
	//-------- Generating a gray scale palette with 20 shades -----------------------------------------------------------
	Int_t grayscale[20];
	{	for(int i = 0; i < 20; i++) // filling the gray scale with 20 shades of gray
		{	TColor *color = new TColor(ncolors + i, 1 - Float_t(i)/20, 1 - Float_t(i)/20, 1 - Float_t(i)/20, "mygray" + i);
			grayscale[i] = ncolors + i;
		}
	}
	//-------- Finished -------------------------------------------------------------------------------------------------
	
	std::cout << "Drawing ..." << std::endl;

	gStyle->SetTitleFillColor(0);
	gStyle->SetStatColor(0);

	//======== (1) Drawing and saving 'kennz' versus 'zstart' versus 'rstart' ===========================================
	vnamex = "rstart"; // dummy ~ x-variable
	vnamey = "zstart"; // dummy ~ y-variable
	vnamec = "kennz";  // dummy ~ z-variable

	TCanvas *c1 = new TCanvas("c1", vnamec + ":" + vnamey + ":" + vnamex + " data from " + rootfilename, 20, 20, 400, 800);
	// creating a new TCanvas([canvasname], [canvastitle], x, y pixel coordinate, x, y pixel size)
	c1->SetFillColor(0);
	c1->SetBorderMode(0);
	c1->cd(); // select pad
	c1->SetFrameFillColor(0);
	c1->SetFrameBorderMode(0);
	mytree->SetEstimate(nentries); // setting the estimated lenght of V1, V2 and V3
	mytree->Draw(vnamec + ":" + vnamey + ":" + vnamex, pnamec, "goff"); // drawing "[vnamec]:[vnamey]:[vnamex]" without
	                                                                    // graphical output
	g1 = new TGraph2D(nentries, mytree->GetV3(), mytree->GetV2(), mytree->GetV1()); // generating graph and retrieving
	                                                                                // data from the draw command above
	g1->SetTitle(vnamec + ":" + vnamey + ":" + vnamex + ";" + vnamex + " [m];" + vnamey + " [m];" + vnamec);
	// setting title (including axis titles seperated by ";")
	g1->SetMarkerStyle(21);
	g1->SetMarkerSize(0.4);
	gStyle->SetPalette(1, 0); //++++++++ options: 1, 0 ~ coloured palette +++++++++++++++++++++++++++++++++++++++++++++++
	                          //+++++++++++++++++ 20, grayscale ~ gray scale palette with 20 shades +++++++++++++++++++++
	g1->Draw("PCOL"); //++++++++ options: "P" ~ markers, "COL" ~ coloured z-values, "Z" ~ colour palette ++++++++++++++++

	vnamec = rootfilename + "_" + vnamec + "-" + vnamey + "-" + vnamex;
	c1->SaveAs(vnamec + ".cxx");   // saving canvas as macro
	c1->SaveAs(vnamec + "_0.png"); // saving canvas as PNG
//	delete c1; // closing canvas window
//	delete g1; // deleting graph
	//======== (1) Finished =============================================================================================

	//======== (2) Drawing and saving the histogram of 't' ==============================================================
	vnamex = "dt";      // dummy ~ x-variable
	vnamey = "counts"; // dummy ~ y-axis title

	TCanvas *c2 = new TCanvas("c2", vnamex + " data from " + rootfilename, 40, 40, 800, 600);
	// creating a new TCanvas([canvasname], [canvastitle], x, y pixel coordinate, x, y pixel size)
	c2->SetFillColor(0);
	c2->SetBorderMode(0);
	c2->cd(); // select pad
	c2->SetFrameFillColor(0);
	c2->SetFrameBorderMode(0);
	mytree->Draw(vnamex + ">>h2", pnamec, "goff"); // drawing "[vnamex]" and storing the histogram 'h2'
	TH1D *h2 = (TH1D*) gDirectory->Get("h2"); // retrieving histogram 'h2' from the draw command above

	h2->GetXaxis()->SetTitle(vnamex + " [s]");
	//h2->GetXaxis()->SetRangeUser(0, ??);
	h2->GetYaxis()->SetTitle(vnamey);
	c2->SetLogy(1); // defining a logarithmical y-axis
	h2->SetLineWidth(2);
	h2->SetFillColor(17);
	h2->Draw();	
	
	vnamec = rootfilename + "_" + vnamex;
	c2->SaveAs(vnamec + ".cxx");   // saving canvas as macro
	c2->SaveAs(vnamec + "_0.png"); // saving canvas as PNG
//	delete c2; // closing canvas window
//	delete h2; // deleting histogram
	//======== (2) Finished =============================================================================================

	//======== (3) Drawing and saving the histogram of 'H' ==============================================================
	vnamex = "H";        // dummy ~ x-variable
	vnamey = "counts";   // dummy ~ y-axis title
	vnamec = "kennz==6"; // dummy ~ selection criterion

	TCanvas *c3 = new TCanvas("c3", vnamex + " {" + vnamec + "} data from " + rootfilename, 60, 60, 800, 600);
	// creating a new TCanvas([canvasname], [canvastitle], x, y pixel coordinate, x, y pixel size)
	c3->SetFillColor(0);
	c3->SetBorderMode(0);
	c3->cd(); // select pad
	c3->SetFrameFillColor(0);
	c3->SetFrameBorderMode(0);
	mytree->Draw(vnamex + ">>h3", "(" + vnamec + ")&&" + pnamec); // drawing "[vnamex]" (if "[vnamec]")
	                                                              // and storing the histogram 'h3'
	TH1D *h3 = (TH1D*) gDirectory->Get("h3"); // retrieving histogram 'h3' from the draw command above

	h3->GetXaxis()->SetTitle(vnamex + " [neV]");
	//h3->GetXaxis()->SetRangeUser(0, ??);
	h3->GetYaxis()->SetTitle(vnamey);
	h3->SetLineWidth(2);
	h3->SetFillColor(17);
	h3->Draw();
	
	vnamec = rootfilename + "_" + vnamex + "_{" + vnamec + "}";
	c3->SaveAs(vnamec + ".cxx");   // saving canvas as macro
	c3->SaveAs(vnamec + "_0.png"); // saving canvas as PNG
//	delete c3; // closing canvas window
//	delete h3; // deleting histogram
	//======== (3) Finished =============================================================================================

	//======== (4) Drawing and saving the histogram of 'rend' ===========================================================
	vnamex = "rend";     // dummy ~ x-variable
	vnamey = "counts";   // dummy ~ y-axis title
	vnamec = "kennz==6"; // dummy ~ selection criterion
	nbins = 75;          // dummy ~ number of bins
	xmin = 0;            // dummy ~ x-minimum
	xmax = 0.4;          // dummy ~ x-maximum

	TCanvas *c4 = new TCanvas("c4", vnamex + " {" + vnamec + "} data from " + rootfilename, 80, 80, 800, 600);
	// creating a new TCanvas([canvasname], [canvastitle], x, y pixel coordinate, x, y pixel size)
	c4->SetFillColor(0);
	c4->SetBorderMode(0);
	c4->cd();
	c4->SetFrameFillColor(0);
	c4->SetFrameBorderMode(0);
	TH1D *h4 = new TH1D("h4", vnamex + " {" + vnamec + "}", nbins, xmin, xmax * 1.0001);
	// creating a new TH1D([histogramname], [histogramtitle], number of bins, lower edge of the first bin, excluded!
	// upper edge of the last bin)
	mytree->Draw(vnamex + ">>+h4", "(" + vnamec + ")&&" + pnamec); // drawing "[vnamex]" (if "[vnamec]")
	                                                               // and storing the result in 'h4'
	h4->GetXaxis()->SetTitle(vnamex + " [m]");
	h4->GetXaxis()->SetRangeUser(xmin, xmax);
	h4->GetYaxis()->SetTitle(vnamey);
	h4->SetLineWidth(2);
	h4->SetFillColor(17);
	h4->Draw();

	vnamec = rootfilename + "_" + vnamex + "_{" + vnamec + "}";
	c4->SaveAs(vnamec + ".cxx");   // saving canvas as macro
	c4->SaveAs(vnamec + "_0.png"); // saving canvas as PNG
//	delete c4; // closing canvas window
//	delete h4; // deleting histogram
	//======== (4) Finished =============================================================================================

	//======== (5) Drawing and saving the histogram of 'gammaend' =======================================================
	vnamex = "gammaend"; // dummy ~ x-variable
	vnamey = "counts";   // dummy ~ y-axis title
	vnamec = "kennz==6"; // dummy ~ 
	nbins = 75;          // dummy ~ number of bins
	xmin = 0;            // dummy ~ x-minimum
	xmax = 90;           // dummy ~ x-maximum

	TCanvas *c5 = new TCanvas("c5", vnamex + " {" + vnamec + "} data from " + rootfilename, 100, 100, 800, 600);
	// creating a new TCanvas([canvasname], [canvastitle], x, y pixel coordinate, x, y pixel size)
	c5->SetFillColor(0);
	c5->SetBorderMode(0);
	c5->cd(); // select pad
	c5->SetFrameFillColor(0);
	c5->SetFrameBorderMode(0);
	TH1D *h5 = new TH1D("h5", vnamex + " {" + vnamec + "}", nbins, xmin, xmax * 1.0001);
	// creating a new TH1D([histogramname], [histogramtitle], number of bins, lower edge of the first bin, excluded!
	// upper edge of the last bin)
	mytree->Draw(vnamex + ">>+h5", "(" + vnamec + ")&&" + pnamec); // drawing "[vnamex]" (if "[vnamec]")
	                                                               // and storing the result in 'h5'
	h5->GetXaxis()->SetTitle(vnamex + " [°]");
	h5->GetXaxis()->SetRangeUser(xmin, xmax);
	h5->GetYaxis()->SetTitle(vnamey);
	h5->SetLineWidth(2);
	h5->SetFillColor(17);
	h5->Draw();

	vnamec = rootfilename + "_" + vnamex + "_{" + vnamec + "}";
	c5->SaveAs(vnamec + ".cxx");   // saving canvas as macro
	c5->SaveAs(vnamec + "_0.png"); // saving canvas as PNG
//	delete c5; // closing canvas window
//	delete h5; // deleting histogram
	//======== (5) Finished =============================================================================================

	//======== (6) Drawing and saving the histogram of 'NeutEnergie' ====================================================
	vnamex = "NeutEnergie"; // dummy ~ x-variable
	vnamey = "counts";      // dummy ~ y-axis title

	TCanvas *c6 = new TCanvas("c6", vnamex + " data from " + rootfilename, 120, 120, 800, 600);
	// creating a new TCanvas([canvasname], [canvastitle], x, y pixel coordinate, x, y pixel size)
	c6->SetFillColor(0);
	c6->SetBorderMode(0);
	c6->cd(); // select pad
	c6->SetFrameFillColor(0);
	c6->SetFrameBorderMode(0);
	mytree->Draw(vnamex + ">>h6", pnamec); // drawing "[vnamex]" and storing the histogram 'h6'
	TH1D *h6 = (TH1D*) gDirectory->Get("h6"); // retrieving histogram 'h6' from the draw command above

	h6->GetXaxis()->SetTitle(vnamex + " [neV]");
	//h6->GetXaxis()->SetRangeUser(0, ??);
	h6->GetYaxis()->SetTitle(vnamey);
	h6->SetLineWidth(2);
	h6->SetFillColor(17);
	h6->Draw();

	vnamec = rootfilename + "_" + vnamex;
	c6->SaveAs(vnamec + ".cxx");   // saving canvas as macro
	c6->SaveAs(vnamec + "_0.png"); // saving canvas as PNG
//	delete c6; // closing canvas window
//	delete h6; // deleting histogram
	//======== (6) Finished =============================================================================================

	//======== (7) Drawing and saving the histogram of 'H' ==============================================================
	vnamex = "H";      // dummy ~ x-variable
	vnamey = "counts"; // dummy ~ y-axis title

	TCanvas *c7 = new TCanvas("c7", vnamex + " data from " + rootfilename, 140, 140, 800, 600);
	// creating a new TCanvas([canvasname], [canvastitle], x, y pixel coordinate, x, y pixel size)
	c7->SetFillColor(0);
	c7->SetBorderMode(0);
	c7->cd(); // select pad
	c7->SetFrameFillColor(0);
	c7->SetFrameBorderMode(0);
	mytree->Draw(vnamex + ">>h7", pnamec); // drawing "[vnamex]" and storing the histogram 'h7'
	TH1D *h7 = (TH1D*) gDirectory->Get("h7"); // retrieving histogram 'h7' from the draw command above

	h7->GetXaxis()->SetTitle(vnamex + " [neV]");
	//h7->GetXaxis()->SetRangeUser(0, ??);
	h7->GetYaxis()->SetTitle(vnamey);
	h7->SetLineWidth(2);
	h7->SetFillColor(17);
	h7->Draw();

	vnamec = rootfilename + "_" + vnamex;
	c7->SaveAs(vnamec + ".cxx");   // saving canvas as macro
	c7->SaveAs(vnamec + "_0.png"); // saving canvas as PNG
//	delete c7; // closing canvas window
//	delete h7; // deleting histogram
	//======== (7) Finished =============================================================================================

	//======== (8) Drawing and saving 'NeutEnergie' versus 't' ==========================================================
	vnamex = "dt";           // dummy ~ x-variable
	vnamey = "NeutEnergie"; // dummy ~ y-variable
	vnamec = "kennz==6";    // dummy ~ selection criterion

	TCanvas *c8 = new TCanvas("c8", vnamey + ":" + vnamex + " {" + vnamec + "} data from " + rootfilename, 140, 140, 800, 600);
	// creating a new TCanvas([canvasname], [canvastitle], x, y pixel coordinate, x, y pixel size)
	c8->SetFillColor(0);
	c8->SetBorderMode(0);
	c8->cd(); // select pad
	c8->SetFrameFillColor(0);
	c8->SetFrameBorderMode(0);
	mytree->SetEstimate(mytree->GetEntries("(" + vnamec + ")&&" + pnamec)); // setting the estimated lenght of V1 and V2
	mytree->Draw(vnamey + ":" + vnamex, "(" + vnamec + ")&&" + pnamec, "goff");
	// drawing "[vnamey]:[vnamey]" (if "[vnamec]") without graphical output

	g8 = new TGraph(mytree->GetEntries("(" + vnamec + ")&&" + pnamec), mytree->GetV2(), mytree->GetV1());
	// generating graph and retrieving data from the draw command above
	
	g8->SetTitle(vnamey + ":" + vnamex + " {" + vnamec + "}");
	g8->SetMarkerStyle(21);
	g8->SetMarkerSize(0.4);
	g8->SetMarkerColor(4);

	c8->Update(); // necessary command for setting the axis titles
	g8->GetHistogram()->SetXTitle(vnamex + " [s]");
	g8->GetHistogram()->SetYTitle(vnamey + " [neV]");
	
	g8->Draw("AP"); //++++++++ options: "A" ~ axis, "P" ~ markers, "L" ~ a simple line between the points +++++++++++++++

	vnamec = rootfilename + "_" + vnamey + "-" + vnamex + "_{" + vnamec + "}";
	c8->SaveAs(vnamec + ".cxx");   // saving canvas as macro
	c8->SaveAs(vnamec + "_0.png"); // saving canvas as PNG
//	delete c8; // closing canvas window
//	delete g8; // deleting graph
	//======== (8) Finished =============================================================================================

	//======== (9) Drawing and saving the histogram of 'EFeldSkal' and 'BFeldSkal' ======================================
	vnamex = "EFeldSkal"; // dummy ~ x-variable in 'h91', 'h93' and 'h94'
	vnamey = "BFeldSkal"; // dummy ~ x-variable in 'h92' and 'h93', y-variable in 'h94'
	vnamec = "kennz==6";  // dummy ~ selection criterion
	nbins = 10;           // dummy ~ number of bins
	xmin = 0;             // dummy ~ x-minimum in 'h91', 'h92', 'h93' and 'h94', y-minimum in 'h94'
	xmax = 1;             // dummy ~ x-maximum in 'h91', 'h92', 'h93' and 'h94', y-maximum in 'h94'

	TCanvas *c9 = new TCanvas("c9", vnamex + " " + vnamey + " {" + vnamec + "} data from " + rootfilename, 160, 160, 800, 600);
	// creating a new TCanvas([canvasname], [canvastitle], x, y pixel coordinate, x, y pixel size)
	c9->SetFillColor(0);
	c9->SetBorderMode(0);
	c9->Divide(2, 2); // dividing 'c9' into 2*2 pads (numbered like text read)
	c9->cd(1); // select pad 1	
	c9_1->SetFrameFillColor(0);
	c9_1->SetFrameBorderMode(0);
	TH1D *h91 = new TH1D("h91", vnamex + " {" + vnamec + "}", nbins, xmin, xmax * 1.0001);
	// creating a new TH1D([histogramname], [histogramtitle], number of bins, lower edge of the first bin, excluded!
	// upper edge of the last bin)
	mytree->Draw(vnamex + ">>+h91", "(" + vnamec + ")&&" + pnamec); // drawing "[vnamex]" (if "[vnamec]")
	                                                                // and storing the result in 'h91'
	h91->SetStats(0);
	h91->GetXaxis()->SetTitle(vnamex);
	h91->GetXaxis()->SetRangeUser(xmin, xmax);
	h91->GetYaxis()->SetTitle("percentage [%]");
	//h91->GetYaxis()->SetRangeUser(0, 100);
	h91->SetLineWidth(2);
	h91->SetFillColor(4);
	h91->SetFillStyle(3018);
	for(int i = 1; i < (nbins + 1); i++) // resize to percent of total entries
	{	h91->SetBinContent(i, (h91->GetBinContent(i))/nentries*100); // resize the content of bin 'i' to percent
	}
	h91->Draw();

	c9->cd(2); // select pad 2
	c9_2->SetFrameFillColor(0);
	c9_2->SetFrameBorderMode(0);
	TH1D *h92 = new TH1D("h92", vnamey + " {" + vnamec + "}", nbins, xmin, xmax * 1.0001);
	// creating a new TH1D([histogramname], [histogramtitle], number of bins, lower edge of the first bin, excluded!
	// upper edge of the last bin)
	mytree->Draw(vnamey + ">>+h92", "(" + vnamec + ")&&" + pnamec); // drawing "[vnamey]" (if "[vnamec]")
	                                                                // and storing the result in 'h92'
	h92->SetStats(0);
	h92->GetXaxis()->SetTitle(vnamey);
	h92->GetXaxis()->SetRangeUser(0, 1);
	h92->GetYaxis()->SetTitle("percentage [%]");
	//h91->GetYaxis()->SetRangeUser(0, 100);
	h92->SetLineWidth(2);
	h92->SetFillColor(2);
	h92->SetFillStyle(3017);
	for(int i = 1; i < (nbins + 1); i++) // resize to percent of total entries
	{	h92->SetBinContent(i, (h92->GetBinContent(i))/nentries*100); // resize the content of bin 'i' to percent
	}
	h92->Draw();

	c9->cd(3); // select pad 3
	c9_3->SetFrameFillColor(0);
	c9_3->SetFrameBorderMode(0);
	THStack *h93 = new THStack("h93", vnamex + ", " + vnamey + " {" + vnamec + "}");
	// creating a new THStack([histogramname], [histogramtitle])
	h93->Add(h91); // adding 'h91' to 'h93'
	h93->Add(h92); // adding 'h92' to 'h93'
	h93->Draw();
	h93->GetXaxis()->SetTitle(vnamex + ", " + vnamey);
	h93->GetXaxis()->SetRangeUser(xmin, xmax);
	h93->GetYaxis()->SetTitle("percentage [%]");
	//h93->GetYaxis()->SetRangeUser(0, 100);
	h93	->Draw("nostack"); //++++++++ options: "" ~ stack, "nostack" ~ no stack (last added one in the front) +++++++++++

	c9->cd(4); // select pad 4
	c9_4->SetFrameFillColor(0);
	c9_4->SetFrameBorderMode(0);
	TH2D *h94 = new TH2D("h94", vnamey + ":" + vnamex + " {" + vnamec + "}", nbins, xmin, xmax * 1.0001, nbins, xmin, xmax * 1.0001);
	// creating a new TH2D([histogramname], [histogramtitle], number of x-bins, lower edge of the first x-bin, excluded!
	// upper edge of the last x-bin, number of y-bins, lower edge y-bin, excluded! upper edge y-bin)
	mytree->Draw(vnamey + ":" + vnamex + ">>+h94", "(" + vnamec + ")&&" + pnamec);
	// drawing "[vnamey]:[vnamex]" (if "[vnamec]") and storing the result in 'h94'

	h94->GetXaxis()->SetTitle(vnamex);
	h94->GetXaxis()->SetTitleOffset(2.0);
	h94->GetXaxis()->SetRangeUser(xmin, xmax);
	h94->GetYaxis()->SetTitle(vnamey);
	h94->GetYaxis()->SetTitleOffset(2.0);
	h94->GetYaxis()->SetRangeUser(xmin, xmax);
	h94->GetZaxis()->SetTitle("percentage [%]");
	//h94->GetZaxis()->SetRangeUser(0, 100);
	h94->SetFillColor(17);
	for(int i = 1; i < (nbins + 1); i++) // resize to percent of total entries
	{	for(int j = 1; j < (nbins + 1); j++)
		{	h94->SetCellContent(i, j, (h94->GetBinContent(i, j))/nentries*100); // resize the content of bin 'i,j' to percent
		}
	}
	h94->Draw("LEGO1"); //++++++++ options: "LEGO1" ~ lego plot with hidden surface removal +++++++++++++++++++++++++++++
	                    //+++++++++++++++++ "LEGO2" ~ lego plot using colours to show the cell contents +++++++++++++++++

	vnamec = rootfilename + "_" + vnamex + "_" + vnamey;
	c9->SaveAs(vnamec + ".cxx");     // saving canvas as macro
	c9->SaveAs(vnamec + "_0.png");   // saving canvas as PNG
	c9_1->SaveAs(vnamec + "_1.png"); // saving pad 1 as PNG
	c9_2->SaveAs(vnamec + "_2.png"); // saving pad 2 as PNG
	c9_3->SaveAs(vnamec + "_3.png"); // saving pad 3 as PNG
	c9_4->SaveAs(vnamec + "_4.png"); // saving pad 4 as PNG
//	delete c9; // closing canvas window
//	delete h91; // deleting histograms
//	delete h92; 
//	delete h93;
//	delete h94;
	//======== (9) Finished =============================================================================================

	//======== (10) Drawing and saving the histogram of 'kennz' =========================================================
	vnamec = "kennz"; // dummy ~ x-variable
	nbins = 12;       // dummy ~ number of bins
	xmin = 0;         // dummy ~ x-minimum
	xmax = 12;        // dummy ~ x-maximum

	TCanvas *c10 = new TCanvas("c10", vnamec + " data from " + rootfilename, 180, 180, 800, 600);
	// creating a new TCanvas([canvasname], [canvastitle], x, y pixel coordinate, x, y pixel size)
	c10->SetFillColor(0);
	c10->SetBorderMode(0);
	c10->cd(); // select pad
	c10->SetFrameFillColor(0);
	c10->SetFrameBorderMode(0);
	TH1D *h10 = new TH1D("h10", vnamec, nbins, xmin, xmax * 1.0001); // creating a new TH1D([histogramname],
	// [histogramtitle], number of bins, lower edge of the first bin, excluded! upper edge of the last bin)
	mytree->Draw(vnamec + ">>+h10", pnamec); // drawing "[vnamec]" and storing the result in 'h10'

	h10->SetStats(0);
	h10->GetXaxis()->SetTitle(vnamec);
	h10->GetXaxis()->SetRangeUser(0, 12);
	h10->GetYaxis()->SetTitle("percentage [%]");
	h10->SetLineWidth(2);
	h10->SetFillColor(17);
	for(int i = 1; i < 13; i++) // resize to percent of total entries
	{	h10->SetBinContent(i, (h10->GetBinContent(i))/nentries*100); // resize the content of bin 'i' to percent
	}
	h10->Draw();

	vnamec = rootfilename + "_" + vnamec;
	c10->SaveAs(vnamec + ".cxx");   // saving canvas as macro
	c10->SaveAs(vnamec + "_0.png"); // saving canvas as PNG
//	delete c10; // closing canvas window
//	delete h10; // deleting histogram
	//======== (10) Finished ============================================================================================

	//======== (11) Drawing and saving 'zstart' versus 'rstart' =========================================================
	vnamex = "rstart";   // dummy ~ x-variable
	vnamey = "zstart";   // dummy ~ y-variable
	vnamec = "kennz==6"; // dummy ~ selection criterion
	nbins = 75;          // dummy ~ number of bins
	xmin = 0.1;          // dummy ~ x-minimum
	xmax = 0.5;          // dummy ~ x-maximum

	TCanvas *c11 = new TCanvas("c11", vnamey + ":" + vnamex + " {" + vnamec + "} data from " + rootfilename, 200, 200, 500, 800);
	// creating a new TCanvas([canvasname], [canvastitle], x, y pixel coordinate, x, y pixel size)
	c11->SetFillColor(0);
	c11->SetBorderMode(0);
	TH2D *h110 = new TH2D("h110", "reference histogram", nbins, xmin, xmax * 1.0001, nbins, 0, 1.2 * 1.0001);
	// creating a new TH2D([histogramname], [histogramtitle], number of x-bins, lower edge of the first x-bin, excluded!
	// upper edge of the last x-bin, number of y-bins, lower edge y-bin, excluded! upper edge y-bin)
	mytree->Draw(vnamey + ":" + vnamex + ">>+h110", pnamec, "goff"); // drawing "[vnamey]:[vnamex]"
	                                                                 // and storing the result in 'h110'
	c11->cd(); // select pad
	c11->SetFrameFillColor(0);
	c11->SetFrameBorderMode(0);
	TH2D *h111 = new TH2D("h111", vnamey + ":" + vnamex + " {" + vnamec + "}", nbins, xmin, xmax * 1.0001, nbins, 0, 1.2 * 1.0001);
	TH2D *h112 = new TH2D("h112", vnamey + ":" + vnamex + " {" + vnamec + "}", nbins, xmin, xmax * 1.0001, nbins, 0, 1.2 * 1.0001);
	// creating two new TH2D([histogramname], [histogramtitle], number of x-bins, lower edge of the first x-bin, excluded!
	// upper edge of the last x-bin, number of y-bins, lower edge y-bin, excluded! upper edge y-bin)
	mytree->Draw(vnamey + ":" + vnamex + ">>+h111", "(" + vnamec + ")&&" + pnamec);
	// drawing "[vnamey]:[vnamex]" (if "[vnamec]") and storing the result in 'h111'

	h111->SetStats(0);
	h111->GetXaxis()->SetTitle(vnamex + " [m]");
	h111->GetXaxis()->SetRangeUser(xmin, xmax);
	h111->GetYaxis()->SetTitle(vnamey + " [m]");
	h111->GetYaxis()->SetRangeUser(0, 1.2);
	h111->GetZaxis()->SetTitle("percentage per bin [%]");
	h111->GetZaxis()->SetTitleColor(0);
	h111->GetZaxis()->SetTitleOffset(-0.35);
	for(int i = 1; i < (nbins + 1); i++) // resize to percent of bin entries without selection criterion
	{	for(int j = 1; j < (nbins + 1); j++)
		{	if(h111->GetBinContent(i, j) > 0)
			{	h111->SetCellContent(i, j, (h111->GetBinContent(i, j))/(h110->GetBinContent(i, j))*100);
				// resize the content of bin 'i,j' to percent
			}
			h112->SetCellContent(i, j, h111->GetBinContent(i, j));
		}
	}
	gStyle->SetPalette(20, grayscale);
	h111->Draw("CONT4Z"); //++++++++ options: "CONT" ~ contour plot, "CONT4" ~ analog, "Z" ~ colour palette +++++++++++++

	vnamec = rootfilename + "_" + vnamey + "-" + vnamex + "_{" + vnamec + "}";	
	c11->SaveAs(vnamec + "_0g.png"); // saving canvas as PNG ("g" ~ gray scale palette)

	h112->SetStats(0);
	h112->GetXaxis()->SetTitle(vnamex + " [m]");
	h112->GetXaxis()->SetRangeUser(xmin, xmax);
	h112->GetYaxis()->SetTitle(vnamey + " [m]");
	h112->GetYaxis()->SetRangeUser(0, 1.2);
	h112->GetZaxis()->SetTitle("percentage per bin [%]");
	h112->GetZaxis()->SetTitleColor(1);
	h112->GetZaxis()->SetTitleOffset(-0.35);
	gStyle->SetPalette(1, 0);
	h112->Draw("CONT4Z"); //++++++++ options: "CONT" ~ contour plot, "CONT4" ~ analog, "Z" ~ colour palette +++++++++++++

	c11->SaveAs(vnamec + ".cxx");    // saving canvas as macro
	c11->SaveAs(vnamec + "_0c.png"); // saving canvas as PNG ("c" ~ coloured palette)
//	delete c11; // closing canvas window
//	delete h111; // deleting histograms
//	delete h112;
	//======== (11) Finished ============================================================================================
	
	//======== (12) Drawing and saving 'kennz' versus 'zend' versus 'rend' ===========================================
	vnamex = "rend"; // dummy ~ x-variable
	vnamey = "zend"; // dummy ~ y-variable
	vnamec = "kennz";  // dummy ~ z-variable

	TCanvas *c12 = new TCanvas("c12", vnamec + ":" + vnamey + ":" + vnamex + " data from " + rootfilename, 500, 20, 400, 800);
	// creating a new TCanvas([canvasname], [canvastitle], x, y pixel coordinate, x, y pixel size)
	c12->SetFillColor(0);
	c12->SetBorderMode(0);
	c12->cd(); // select pad
	c12->SetFrameFillColor(0);
	c12->SetFrameBorderMode(0);
	mytree->SetEstimate(nentries); // setting the estimated lenght of V1, V2 and V3
	mytree->Draw(vnamec + ":" + vnamey + ":" + vnamex, pnamec, "goff"); // drawing "[vnamec]:[vnamey]:[vnamex]" without
	                                                                    // graphical output
	g12 = new TGraph2D(nentries, mytree->GetV3(), mytree->GetV2(), mytree->GetV1()); // generating graph and retrieving
	                                                                                // data from the draw command above
	g12->SetTitle(vnamec + ":" + vnamey + ":" + vnamex + ";" + vnamex + " [m];" + vnamey + " [m];" + vnamec);
	// setting title (including axis titles seperated by ";")
	g12->SetMarkerStyle(21);
	g12->SetMarkerSize(0.4);
	gStyle->SetPalette(1, 0); //++++++++ options: 1, 0 ~ coloured palette +++++++++++++++++++++++++++++++++++++++++++++++
	                          //+++++++++++++++++ 20, grayscale ~ gray scale palette with 20 shades +++++++++++++++++++++
	g12->Draw("PCOL"); //++++++++ options: "P" ~ markers, "COL" ~ coloured z-values, "Z" ~ colour palette ++++++++++++++++

	vnamec = rootfilename + "_" + vnamec + "-" + vnamey + "-" + vnamex;
	c12->SaveAs(vnamec + ".cxx");   // saving canvas as macro
	c12->SaveAs(vnamec + "_0.png"); // saving canvas as PNG
//	delete c12; // closing canvas window
//	delete g1; // deleting graph
	//======== (1) Finished =============================================================================================

	std::cout << "Drawing done." << std::endl;

	//-------- Removing the 20 gray shades from the list of colours -----------------------------------------------------
	{	for(int i = 0; i < 20; i++) // removing the last 20 colours
		{	gROOT->GetListOfColors()->RemoveLast();
		}
	}
	//-------- Finished -------------------------------------------------------------------------------------------------

//	file->Close();  // closing the file
}
