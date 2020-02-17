void PlotResults(const char* inputFile = "analysisdata.root") {

    TString inFileName(inputFile);
    TFile *fIn = new TFile(inFileName, "read");

    TH1D *hVobs;
    fIn->GetObject("hVobs", hVobs);

    TH1D *hRsub;
    fIn->GetObject("hRsub", rSub);
}
