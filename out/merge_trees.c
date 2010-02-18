// merge trees from different jobs ("xxxxxxxxend.out.root") into one tree ("end.out.root")

void merge_trees(){
	TSystemDirectory dir(".",".");
	TList *files = dir.GetListOfFiles();
	TChain endlogs("EndTree", "end log");
	TChain tracks("TrackTree", "track log");
	for (Int_t i = 0; ; i++){
		TObject *f;
		if (f = files->At(i)){
			TString filename = f->GetName();
			if (filename.EndsWith("end.out.root") && TString(filename(0,8)).IsDigit()){
				cout << "Adding " << filename << endl;
				endlogs.Add(filename.Data());
				tracks.Add(filename.Data());
			}
		}
		else break;
	}
	cout << "Merging trees, read entries: " << endlogs.GetEntries()+tracks.GetEntries() << endl;
	endlogs.Merge("end.out.root");
	tracks.Merge("end.out.root");
	cout << "Done, written entries: " << endlogs.GetEntries()+tracks.GetEntries() << endl;
}