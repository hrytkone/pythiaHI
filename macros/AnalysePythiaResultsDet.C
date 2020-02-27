#include <boost/math/special_functions/bessel.hpp>

#define N 5

double GetEta(double px, double py, double pz);
double GetPhi(double x, double y);
void CalculateQvector(TComplex &Qvec, double &norm, int n, double phi, double pt, bool bUseWeight);
double GetEventPlane(TComplex Qvec, int n);
double CalculateVnobs(vector<double> phiList, TComplex &Qvec, int n);
double GetVnObs(TComplex Qvec, double phi, int n);

double R1(double khi);
double func(double *x, double *p);
double RIter(double x0, double R0, double err);
double CalculateRerror(double khi, double khiErr);

void AnalysePythiaResultsDet(const char* inputFile = "pythia.root") {
    TString inFileName(inputFile);
    TFile *fIn = new TFile(inFileName, "read");

    TString outFileName("analysisdata.root");
    TFile *fOut = new TFile(outFileName, "recreate");

    TNtuple *ntuple;
    fIn->GetObject("pythiaEvents", ntuple);

    // Event variables
    Float_t eventid = 0;

    // Track variables
    Float_t particleid = 0;
    Float_t charge = 0;
    Float_t px = 0.0, py = 0.0, pz = 0.0, x = 0.0, y = 0.0, z = 0.0;

    ntuple->SetBranchAddress("particleId",&particleid);
    ntuple->SetBranchAddress("eventId",&eventid);
    ntuple->SetBranchAddress("charge",&charge);
    ntuple->SetBranchAddress("px",&px);
    ntuple->SetBranchAddress("py",&py);
    ntuple->SetBranchAddress("pz",&pz);
    ntuple->SetBranchAddress("x",&x);
    ntuple->SetBranchAddress("y",&y);
    ntuple->SetBranchAddress("z",&z);

    TH1D *hVobs[N];
    TH1D *hRsub[N];

    for (int i=0; i<N; i++) {
        hVobs[i] = new TH1D(Form("hV%dobs", i+1), Form("hV%dobs", i+1), 100, -2.0, 2.0);
        hRsub[i] = new TH1D(Form("hRsub%d", i+1), Form("hRsub%d", i+1), 100, -1.0, 1.0);
    }

    TComplex Qvec[N], QvecA[N], QvecB[N];
    double norm[N], normA[N], normB[N];
    vector<double> phiList;

    double vobs = 0.0, epA = 0.0, epB = 0.0, rsub = 0.0;

    Int_t previd = 0, nevents = 0;
    Int_t nentries = (Int_t)ntuple->GetEntries();

    int noutput = nentries/20;
    if (noutput<1) noutput = 1;
    for ( Int_t i=0; i<nentries; i++ ) {
    //for ( Int_t i=0; i<10; i++ ) {
        if (i % noutput == 0)
            cout << 100*(double)i/(double)nentries << " % finished" << endl;

        Int_t entry = ntuple->GetEntry(i);
        if (eventid!=previd) {
            previd = eventid;

            for (Int_t j=0; j<N; j++) {
                vobs = CalculateVnobs(phiList, Qvec[j], j+1);
                vobs /= phiList.size();
                hVobs[j]->Fill(vobs);

                epA = GetEventPlane(QvecA[j], j+1);
                epB = GetEventPlane(QvecB[j], j+1);
                rsub = TMath::Cos((j+1)*(epA - epB));

                hRsub[j]->Fill(rsub);

                norm[j] = 0; normA[j] = 0; normB[j] = 0;
                Qvec[j] = TComplex(0, 0); QvecA[j] = TComplex(0, 0); QvecB[j] = TComplex(0, 0);
            }

            phiList.clear();
        }

        double phi = GetPhi(px, py);
        double eta = GetEta(px, py, pz);

        for (Int_t j=0; j<N; j++) {
            if ((eta > 2.2) && (eta < 5.06)) {
                CalculateQvector(Qvec[j], norm[j], j+1, phi, 0, 0);
                phiList.push_back(phi);
            }
            if ((eta > 3.8) && (eta < 5.4)) {
                CalculateQvector(QvecA[j], normA[j], j+1, phi, 0, 0);
            }
            if ((eta > -3.3) && (eta < -2.2)) {
                CalculateQvector(QvecB[j], normB[j], j+1, phi, 0, 0);
            }
        }
    }

    for (Int_t j=0; j<N; j++) {
        vobs = CalculateVnobs(phiList, Qvec[j], j+1);
        vobs /= phiList.size();
        hVobs[j]->Fill(vobs);

        epA = GetEventPlane(QvecA[j], j+1);
        epB = GetEventPlane(QvecB[j], j+1);
        rsub = TMath::Cos((j+1)*(epA - epB));
        hRsub[j]->Fill(rsub);
    }

    for (Int_t i=0; i<N; i++) {

        double rinit = TMath::Sqrt(hRsub[i]->GetMean());
        double khi = RIter(0.5, rinit, 0.0001);
        double res = R1(TMath::Sqrt(2)*khi);

        cout << "R" << i+1 << " : " << res;
        cout << "   v" << i+1 << " : " << hVobs[i]->GetMean()/res << endl;

        hVobs[i]->Write(Form("hV%dobs", i+1));
        hRsub[i]->Write(Form("hRsub%d", i+1));
    }
}

//______________________________________________________________________________
double GetEta(double px, double py, double pz) {
    double p = TMath::Sqrt(px*px + py*py + pz*pz);
    return TMath::ATanH(pz/p);
}


double GetPhi(double x, double y) {
    return TMath::ATan2(y, x);
}

void CalculateQvector(TComplex &Qvec, double &norm, int n, double phi, double pt, bool bUseWeight) {

    double w = 1.0;
    if (bUseWeight) w = pt;
    norm += w*w;

    Qvec += TComplex(w*TMath::Cos(n*phi), w*TMath::Sin(n*phi));
}

double GetEventPlane(TComplex Qvec, int n) {
    return TMath::ATan2(Qvec.Im(), Qvec.Re())/n;
}

double CalculateVnobs(vector<double> phiList, TComplex &Qvec, int n) {

    double vobs = 0.0;
    double w = 1.0;
    for (int i=0; i<phiList.size(); i++) {

        //if (bUseWeight) w = pt;
        double phi = phiList[i];
        Qvec -= TComplex(w*TMath::Cos(n*phi), w*TMath::Sin(n*phi));
        vobs += GetVnObs(Qvec, phi, n);
        Qvec += TComplex(w*TMath::Cos(n*phi), w*TMath::Sin(n*phi));
    }

    return vobs;
}

double GetVnObs(TComplex Qvec, double phi, int n) {
    return TMath::Cos(n*(phi - GetEventPlane(Qvec, n)));
}

double R1(double khi) {
    return (TMath::Sqrt(TMath::Pi())/2)*khi*TMath::Exp(-khi*khi/2)*(TMath::BesselI0(khi*khi/2) + TMath::BesselI1(khi*khi/2));
}

double func(double *x, double *p) {
    double khi = x[0];
    double Rk = p[0];
    return (TMath::Sqrt(TMath::Pi())/2)*khi*TMath::Exp(-khi*khi/2)*(TMath::BesselI0(khi*khi/2) + TMath::BesselI1(khi*khi/2)) - Rk;
}

double RIter(double x0, double R0, double err) {
    double x = 0;
    TF1 *fRes = new TF1("fRes", func, 0, 50.0, 1);
    fRes->SetParameter(0, R0);
    while (TMath::Abs(R1(x) - R0) > err) {
        x = x0 - fRes->Eval(x0)/fRes->Derivative(x0);
        x0 = x;
    }
    return x;
}

// Virheen yleisellä etenemisellä R(khi):n lausekkeesta
double CalculateRerror(double khi, double khiErr) {
    double bessel0 = -(khi*khi-2)*TMath::BesselI0(khi*khi/2);
    double bessel1 = 2*TMath::BesselI1(khi*khi/2);
    double bessel2 = khi*khi*TMath::BesselI(2,khi*khi/2);
    return TMath::Sqrt(TMath::Pi()/4)*TMath::Exp(-khi*khi)*(bessel0 + bessel1 + bessel2)*khiErr;
}
