void CentralityPlot(const char* inputFile = "o2sim.root") {

    TString inFileName(inputFile);
    TFile *fIn = new TFile(inFileName, "read");
    TTree *tree = (TTree*)fIn->Get("o2sim");
    Int_t nEntries = (Int_t)tree->GetEntries();

    vector<o2::fv0::Hit>* hitArrFV0 = nullptr;
    tree->SetBranchAddress("FV0Hit", &hitArrFV0);

    vector<double> hitList;

    TH1D *hCharge = new TH1D("hCharge", "hCharge", 100, 0., 1.);

    int maxVal = 0;

    int nOutput = nEntries/20;
    if (nOutput<1) nOutput = 1;
    for (int iev=0; iev<nEntries; iev++) {
        if (iev % nOutput == 0)
            cout << 100*iev/nEntries << " % finished" << endl;

        tree->GetEntry(iev);

        Int_t arrSize = hitArrFV0->size();

        if (arrSize > maxVal) maxVal = arrSize;

        hitList.push_back(arrSize);
    }

    for (int i=0; i<hitList.size(); i++) {
        hCharge->Fill((double)hitList[i]/(double)maxVal);
    }

    TCanvas *c1 = new TCanvas("c1", "c1");
    hCharge->Draw();
}
