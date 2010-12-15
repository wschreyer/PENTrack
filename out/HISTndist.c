// Create ndist histogram of all ndist files in current directory

void HISTndist(){
	TSystemDirectory dir(".",".");
	TList *files = dir.GetListOfFiles();
	TFile *outfile = new TFile("ndist.root","RECREATE");
	TH2D *ndist = new TH2D("ndist", "Neutron Distribution", 301, 0.001, 0.601, 1201, -0.499, 1.901);
	ifstream infile;
	Int_t xi, yi, hits;
	Double_t x, y, prob;
	for (Int_t i = 0; ; i++){
		TObject *f;
		if (f = files->At(i)){
			TString filename = f->GetName();
			if (filename.EndsWith("ndist.out") && TString(filename(0,8)).IsDigit()){
				cout << "Adding " << filename << endl;
				infile.open(filename.Data());
				infile.ignore(9999, '\n');
				while (infile){
					infile >> xi >> yi >> x >> y >> prob >> hits;
					if (infile && (hits != 0)) ndist->Fill(x,y,prob);
				}
				infile.close();
				cout << "Entries: " << ndist->GetEntries() << endl;
			}
		}
		else break;
	}
	cout << "Done!" << endl;
	ndist->Write();
	outfile->Close();
	ndist->Draw();
}

