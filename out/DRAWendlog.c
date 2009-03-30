//
//	DRAWendlog.c is a ROOT scipt that draws a few plots from a ROOT tree named "[filename].root" and saves them as macros.
//	It is highly recommended to run TREEendlog.c before with the same data file.
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
//
//	User Instructions:
//	(1) Make sure the ROOT tree exists and is named "mytree"!
//	(2) Make sure the ROOT tree contains at least the branches named in the list of plots above!
//	(3) Make sure the branches contain doubles!
//	(4) Run the script in ROOT like this: root [0] .x DRAWendlog.c("[filename]");
//
//	Note:
//	(1) Any occurrence of errors may possibly lead to a loss of data.
//	(2) The created macros are named "[filename].root_[vnamec].C" (e.g. "[filename].root_kennz.C").
//	(3) If the macro files (or files with the same name) allready exist, they will by overwriten.
//	(4) Alterable options in the script are marked with at least 8 '+' characters.
//	(5) At the beginning of each draw routine you can find the assignment of a few dummies, which offer an easy way to
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

	const Long_t nentries = mytree->GetEntries(); // total number of entries in 'mytree'

	std::cout << "Drawing ..." << std::endl;

	//======== (1) Drawing and saving the histogram of 'kennz' versus 'zstart' versus 'rstart' ==========================
	vnamex = "rstart"; // dummy ~ x-variable
	vnamey = "zstart"; // dummy ~ y-variable
	vnamec = "kennz";  // dummy ~ z-variable

	TCanvas *c1 = new TCanvas("c1", vnamec + ":" + vnamey + ":" + vnamex + " data from " + rootfilename, 20, 20, 800, 600);
	// creating a new TCanvas([canvasname], [canvastitle], x, y pixel coordinate, x, y pixel size)
	c1->cd(); // select pad
	mytree->SetEstimate(nentries); // setting the estimated lenght of V1, V2 and V3
	mytree->Draw(vnamec + ":" + vnamey + ":" + vnamex, "", "goff"); // drawing "[vnamec]:[vnamey]:[vnamex]" without
	                                                                // graphical output
	g1 = new TGraph2D(nentries, mytree->GetV3(), mytree->GetV2(), mytree->GetV1()); // generating graph and retrieving
	                                                                                // data from the draw command above
	gStyle->SetPalette(1);
	//c1->UseCurrentStyle();
	//g1->UseCurrentStyle();
	g1->SetTitle("kennz:zstart:rstart;rstart [m];zstart [m];kennz"); // setting title (including axis titles seperated by ";")
	g1->SetMarkerStyle(21);
	g1->SetMarkerSize(0.4);
	g1->Draw("PCOL"); //++++++++ options: "P" ~ markers, "COL" ~ coloured z-values, "Z" ~ colour palette ++++++++++++++++

	vnamec = rootfilename + "_" + vnamec + "-" + vnamey + "-" + vnamex + ".C";
	c1->SaveAs(vnamec); // saving canvas
//	delete c1; // closing canvas window
//	delete g1; // deleting graph
	//======== (1) Finished =============================================================================================

	//======== (2) Drawing and saving the histogram of 't' ==============================================================
	vnamex = "t";      // dummy ~ x-variable
	vnamey = "counts"; // dummy ~ y-axis title

	TCanvas *c2 = new TCanvas("c2", vnamex + " data from " + rootfilename, 40, 40, 800, 600);
	// creating a new TCanvas([canvasname], [canvastitle], x, y pixel coordinate, x, y pixel size)
	c2->cd(); // select pad
	mytree->Draw(vnamex + ">>h2", "", "goff"); // drawing "[vnamex]" and storing the histogram 'h2'
	TH1D *h2 = (TH1D*) gDirectory->Get("h2"); // retrieving histogram 'h2' from the draw command above

	h2->GetXaxis()->SetTitle(vnamex + " [s]");
	//h2->GetXaxis()->SetRangeUser(0, ??);
	h2->GetYaxis()->SetTitle(vnamey);
	c2->SetLogy(1); // defining a logarithmical y-axis
	h2->SetLineWidth(2);
	h2->SetFillColor(17);
	h2->Draw();	
	
	vnamec = rootfilename + "_" + vnamex + ".C";
	c2->SaveAs(vnamec); // saving canvas
//	delete c2; // closing canvas window
//	delete h2; // deleting histogram
	//======== (2) Finished =============================================================================================

	//======== (3) Drawing and saving the histogram of 'H' ==============================================================
	vnamex = "H";        // dummy ~ x-variable
	vnamey = "counts";   // dummy ~ y-axis title
	vnamec = "kennz==6"; // dummy ~ selection criterion

	TCanvas *c3 = new TCanvas("c3", vnamex + " {" + vnamec + "} data from " + rootfilename, 60, 60, 800, 600);
	// creating a new TCanvas([canvasname], [canvastitle], x, y pixel coordinate, x, y pixel size)
	c3->cd(); // select pad
	mytree->Draw(vnamex + ">>h3", vnamec); // drawing "[vnamex]" (if "[vnamec]") and storing the histogram 'h3'
	TH1D *h3 = (TH1D*) gDirectory->Get("h3"); // retrieving histogram 'h3' from the draw command above

	h3->GetXaxis()->SetTitle(vnamex + " [neV]");
	//h3->GetXaxis()->SetRangeUser(0, ??);
	h3->GetYaxis()->SetTitle(vnamey);
	h3->SetLineWidth(2);
	h3->SetFillColor(17);
	h3->Draw();
	
	vnamec = rootfilename + "_" + vnamex + "_{" + vnamec + "}.C";
	c3->SaveAs(vnamec); // saving canvas
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
	c4->cd();
	TH1D *h4 = new TH1D("h4", vnamex + " {" + vnamec + "}", nbins, xmin, xmax * 1.0001);
	// creating a new TH1D([histogramname], [histogramtitle], number of bins, lower edge of the first bin, excluded!
	// upper edge of the last bin)
	mytree->Draw(vnamex + ">>+h4", vnamec); // drawing "[vnamex]" (if "[vnamec]") and storing the result in 'h4'

	h4->GetXaxis()->SetTitle(vnamex + " [m]");
	h4->GetXaxis()->SetRangeUser(xmin, xmax);
	h4->GetYaxis()->SetTitle(vnamey);
	h4->SetLineWidth(2);
	h4->SetFillColor(17);
	h4->Draw();

	vnamec = rootfilename + "_" + vnamex + "_{" + vnamec + "}.C";
	c4->SaveAs(vnamec); // saving canvas
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
	c5->cd(); // select pad
	TH1D *h5 = new TH1D("h5", vnamex + " {" + vnamec + "}", nbins, xmin, xmax * 1.0001);
	// creating a new TH1D([histogramname], [histogramtitle], number of bins, lower edge of the first bin, excluded!
	// upper edge of the last bin)
	mytree->Draw(vnamex + ">>+h5", vnamec); // drawing "[vnamex]" (if "[vnamec]") and storing the result in 'h5'

	h5->GetXaxis()->SetTitle(vnamex + " [°]");
	h5->GetXaxis()->SetRangeUser(xmin, xmax);
	h5->GetYaxis()->SetTitle(vnamey);
	h5->SetLineWidth(2);
	h5->SetFillColor(17);
	h5->Draw();

	vnamec = rootfilename + "_" + vnamex + "_{" + vnamec + "}.C";
	c5->SaveAs(vnamec); // saving canvas
//	delete c5; // closing canvas window
//	delete h5; // deleting histogram
	//======== (5) Finished =============================================================================================

	//======== (6) Drawing and saving the histogram of 'NeutEnergie' ====================================================
	vnamex = "NeutEnergie"; // dummy ~ x-variable
	vnamey = "counts";      // dummy ~ y-axis title

	TCanvas *c6 = new TCanvas("c6", vnamex + " data from " + rootfilename, 120, 120, 800, 600);
	// creating a new TCanvas([canvasname], [canvastitle], x, y pixel coordinate, x, y pixel size)
	c6->cd(); // select pad
	mytree->Draw(vnamex + ">>h6"); // drawing "[vnamex]" and storing the histogram 'h6'
	TH1D *h6 = (TH1D*) gDirectory->Get("h6"); // retrieving histogram 'h6' from the draw command above

	h6->GetXaxis()->SetTitle(vnamex + " [neV]");
	//h6->GetXaxis()->SetRangeUser(0, ??);
	h6->GetYaxis()->SetTitle(vnamey);
	h6->SetLineWidth(2);
	h6->SetFillColor(17);
	h6->Draw();

	vnamec = rootfilename + "_" + vnamex + ".C";
	c6->SaveAs(vnamec); // saving canvas
//	delete c6; // closing canvas window
//	delete h6; // deleting histogram
	//======== (6) Finished =============================================================================================

	//======== (7) Drawing and saving the histogram of 'H' ==============================================================
	vnamex = "H";      // dummy ~ x-variable
	vnamey = "counts"; // dummy ~ y-axis title

	TCanvas *c7 = new TCanvas("c7", vnamex + " data from " + rootfilename, 140, 140, 800, 600);
	// creating a new TCanvas([canvasname], [canvastitle], x, y pixel coordinate, x, y pixel size)
	c7->cd(); // select pad
	mytree->Draw(vnamex + ">>h7"); // drawing "[vnamex]" and storing the histogram 'h7'
	TH1D *h7 = (TH1D*) gDirectory->Get("h7"); // retrieving histogram 'h7' from the draw command above

	h7->GetXaxis()->SetTitle(vnamex + " [neV]");
	//h7->GetXaxis()->SetRangeUser(0, ??);
	h7->GetYaxis()->SetTitle(vnamey);
	h7->SetLineWidth(2);
	h7->SetFillColor(17);
	h7->Draw();

	vnamec = rootfilename + "_" + vnamex + ".C";
	c7->SaveAs(vnamec); // saving canvas
//	delete c7; // closing canvas window
//	delete h7; // deleting histogram
	//======== (7) Finished =============================================================================================

	//======== (8) Drawing and saving 'NeutEnergie' versus 't' ==========================================================
	vnamex = "t";           // dummy ~ x-variable
	vnamey = "NeutEnergie"; // dummy ~ y-variable
	vnamec = "kennz==6";    // dummy ~ selection criterion

	TCanvas *c8 = new TCanvas("c8", vnamey + ":" + vnamex + " {" + vnamec + "} data from " + rootfilename, 140, 140, 800, 600);
	// creating a new TCanvas([canvasname], [canvastitle], x, y pixel coordinate, x, y pixel size)
	c8->cd(); // select pad
	mytree->SetEstimate(mytree->GetEntries(vnamec)); // setting the estimated lenght of V1 and V2
	mytree->Draw(vnamey + ":" + vnamex, vnamec, "goff"); // drawing "[vnamey]:[vnamey]" (if "[vnamec]")
	                                                     // without graphical output
	g8 = new TGraph(mytree->GetEntries(vnamec), mytree->GetV2(), mytree->GetV1()); // generating graph and retrieving
	                                                                               // data from the draw command above
	g8->SetTitle(vnamey + ":" + vnamex + " {" + vnamec + "}");
	g8->SetMarkerStyle(21);
	g8->SetMarkerSize(0.4);
	g8->SetMarkerColor(4);
	g8->Draw("AP"); //++++++++ options: "A" ~ axis, "P" ~ markers, "L" ~ a simple line between the points +++++++++++++++

	c8->Update(); // necessary command for setting the axis titles
	g8->GetHistogram()->SetXTitle(vnamex + " [s]");
	g8->GetHistogram()->SetYTitle(vnamey + " [neV]");

	vnamec = rootfilename + "_" + vnamey + "-" + vnamex + "_{" + vnamec + "}.C";
	c8->SaveAs(vnamec); // saving canvas
//	delete c8; // closing canvas window
//	delete g8; // deleting graph
	//======== (8) Finished =============================================================================================

	//======== (9) Drawing and saving the histogram of 'EFeldSkal' and 'BFeldSkal' ======================================
	vnamex = "EFeldSkal"; // dummy ~ x-variable in h91, h93 and h94
	vnamey = "BFeldSkal"; // dummy ~ x-variable in h92 and h93, y-variable in h94
	vnamec = "kennz==6";  // dummy ~ selection criterion
	nbins = 10;           // dummy ~ number of bins
	xmin = 0;             // dummy ~ x-minimum in h91, h92, h93 and h94, y-minimum in h94
	xmax = 1;             // dummy ~ x-maximum in h91, h92, h93 and h94, y-maximum in h94

	TCanvas *c9 = new TCanvas("c9", vnamex + " " + vnamey + " {" + vnamec + "} data from " + rootfilename, 160, 160, 800, 600);
	// creating a new TCanvas([canvasname], [canvastitle], x, y pixel coordinate, x, y pixel size)
	c9->Divide(2,2); // dividing 'c9' into 2*2 pads (numbered like text read)
	c9->cd(1); // select pad 1
	TH1D *h91 = new TH1D("h91",  vnamex + " {" + vnamec + "}", nbins, xmin, xmax * 1.0001);
	// creating a new TH1D([histogramname], [histogramtitle], number of bins, lower edge of the first bin, excluded!
	// upper edge of the last bin)
	mytree->Draw(vnamex + ">>+h91", vnamec); // drawing "[vnamex]" (if "[vnamec]") and storing the result in 'h91'

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
	TH1D *h92 = new TH1D("h92", vnamey + " {" + vnamec + "}", nbins, xmin, xmax * 1.0001);
	// creating a new TH1D([histogramname], [histogramtitle], number of bins, lower edge of the first bin, excluded!
	// upper edge of the last bin)
	mytree->Draw(vnamey + ">>+h92", vnamec); // drawing "[vnamey]" (if "[vnamec]") and storing the result in 'h92'

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
	TH2D *h94 = new TH2D("h94", vnamey + ":" + vnamex + " {" + vnamec + "}", nbins, xmin, xmax * 1.0001, nbins, xmin, xmax * 1.0001);
	// creating a new TH2D([histogramname], [histogramtitle], number of x-bins, lower edge of the first x-bin, excluded!
	// upper edge of the last x-bin, number of y-bins, lower edge y-bin, excluded! upper edge y-bin)
	mytree->Draw(vnamey + ":" + vnamex + ">>+h94", vnamec); // drawing "[vnamey]:[vnamex]" (if "[vnamec]") and storing
	                                                        // the result in 'h94'
	h94->GetXaxis()->SetTitle(vnamex);
	h94->GetXaxis()->SetTitleOffset(2.0);
	h94->GetXaxis()->SetRangeUser(xmin, xmax);
	h94->GetYaxis()->SetTitle(vnamey);
	h94->GetYaxis()->SetTitleOffset(2.0);
	h94->GetYaxis()->SetRangeUser(xmin, xmax);
	h94->GetZaxis()->SetTitle("percentage [%]");
	//h94->GetZaxis()->SetTitleOffset(2.0);
	//h94->GetZaxis()->SetRangeUser(0, 100);
	h94->SetFillColor(17);
	for(int i = 1; i < (nbins + 1); i++) // resize to percent of total entries
	{	for(int j = 1; j < (nbins + 1); j++)
		{	h94->SetCellContent(i, j, (h94->GetBinContent(i, j))/nentries*100); // resize the content of bin 'i,j' to percent
		}
	}
	h94->Draw("LEGO1"); //++++++++ options: "LEGO1" ~ lego plot with hidden surface removal +++++++++++++++++++++++++++++

	vnamec = rootfilename + "_" + vnamex + "_" + vnamey + ".C";
	c9->SaveAs(vnamec); // saving canvas
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
	c10->cd(); // select pad
	TH1D *h10 = new TH1D("h10", vnamec, nbins, xmin, xmax * 1.0001); // creating a new TH1D([histogramname],
	// [histogramtitle], number of bins, lower edge of the first bin, excluded! upper edge of the last bin)
	mytree->Draw(vnamec + ">>+h10"); // drawing "[vnamec]" and storing the result in 'h10'

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

	vnamec = rootfilename + "_" + vnamec + ".C";
	c9->SaveAs(vnamec); // saving canvas
//	delete c9; // closing canvas window
//	delete h10; // deleting histogram
	//======== (10) Finished ============================================================

	std::cout << "Drawing done." << std::endl;

//	file->Close();  // closing the file
}
