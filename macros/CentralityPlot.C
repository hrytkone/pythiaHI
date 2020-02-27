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

    TF1 *fNBD = new TF1("fNBD", "([0]*TMath::Gamma(x+[1])/(TMath::Gamma(x+1)*TMath::Gamma([1])))*(TMath::Power([2]/[1], x)/TMath::Power([2]/[1] + 1, x + [1]))", 0.1, 1.5);
    fNBD->SetParameter(0, 1.0);
    fNBD->SetParameter(1, 1.0);
    fNBD->SetParameter(2, 1.0);

    TF1 *fNBD2 = new TF1("fNBD2","[2]*ROOT::Math::negative_binomial_pdf(x,[0],[1])", 0.0, 1.0);
    fNBD2->SetParameter(0, 1.0);
    fNBD2->SetParameter(1, 1.0);
    fNBD2->SetParameter(2, 1.0);

    hCharge->Fit("fNBD2", "R");

    TCanvas *c1 = new TCanvas("c1", "c1");
    //hCharge->Draw();
    fNBD2->Draw();
}
