#include <boost/math/special_functions/bessel.hpp>

#include "DataFormatsFV0/Hit.h"
#include "DataFormatsFT0/HitType.h"
#include "FT0Base/Geometry.h"

#define N 5
#define FT0A_CH_N 96
#define FT0C_CH_N 112

double GetFT0Eta(int chID, double z);
double GetFV0Eta(int chID);
double GetFV0Phi(int chID);
double GetFT0Phi(int chID);
double GetFT0ChannelCenterX(int chID);
double GetFT0ChannelCenterY(int chID);

void CalculateQvector(TComplex &Qvec, double &norm, int n, double phi, double pt, bool bUseWeight);
double GetEventPlane(TComplex Qvec, int n);
double CalculateResolution(double sumAB, double sumAC, double sumBC);
double CalculateVnobs(vector<double> phiList, TComplex &Qvec, int n);
double GetVnObs(TComplex Qvec, double phi, int n);

double R1(double khi);
double func(double *x, double *p);
double RIter(double x0, double R0, double err);
double CalculateRerror(double khi, double khiErr);

const double xFT0A[FT0A_CH_N] = {-13.2745, -13.2745, -10.3255, -10.3255, -7.3745, -7.3745,
                    -4.4255, -4.4255, -1.4745, -1.4745, 1.4745, 1.4745, 4.4255,
                    4.4255, 7.3745, 7.3745, 10.3255, 10.3255, 13.2745, 13.2745,
                    -13.2745, -13.2745, -10.3255, -10.3255, -7.3745, -7.3745,
                    -4.4255, -4.4255, -1.4745, -1.4745, 1.4745, 1.4745, 4.4255,
                    4.4255, 7.3745, 7.3745, 10.3255, 10.3255, 13.2745, 13.2745,
                    -14.2745, -14.2745, -11.3255, -11.3255, -8.3745, -8.3745,
                    -5.4255, -5.4255, 5.4255, 5.4255, 8.3745, 8.3745, 11.3255,
                    11.3255, 14.2745, 14.2745, -13.2745, -13.2745, -10.3255,
                    -10.3255, -7.3745, -7.3745, -4.4255, -4.4255, -1.4745,
                    -1.4745, 1.4745, 1.4745, 4.4255, 4.4255, 7.3745, 7.3745,
                    10.3255, 10.3255, 13.2745, 13.2745, -13.2745, -13.2745,
                    -10.3255, -10.3255, -7.3745, -7.3745, -4.4255, -4.4255,
                    -1.4745, -1.4745, 1.4745, 1.4745, 4.4255, 4.4255, 7.3745,
                    7.3745, 10.3255, 10.3255, 13.2745, 13.2745};

const double yFT0A[FT0A_CH_N] = {10.4255, 13.3745, 13.3745, 10.4255, 10.4255, 13.3745,
                    13.3745, 10.4255, 11.4255, 14.3745, 14.3745, 11.4255,
                    10.4255, 13.3745, 13.3745, 10.4255, 10.4255, 13.3745,
                    13.3745, 10.4255, 4.5255, 7.4745, 7.4745, 4.5255,
                    4.5255, 7.4745, 7.4745, 4.5255, 5.5255, 8.4745, 8.4745,
                    5.5255, 4.5255, 7.4745, 7.4745, 4.5255, 4.5255, 7.4745,
                    7.4745, 4.5255, -1.4745, 1.4745, 1.4745, -1.4745, -1.4745,
                    1.4745, 1.4745, -1.4745, -1.4745, 1.4745, 1.4745, -1.4745,
                    -1.4745, 1.4745, 1.4745, -1.4745, -7.4745, -4.5255, -4.5255,
                    -7.4745, -7.4745, -4.5255, -4.5255, -7.4745, -8.4745,
                    -5.5255, -5.5255, -8.4745, -7.4745, -4.5255, -4.5255,
                    -7.4745, -7.4745, -4.5255, -4.5255, -7.4745, -13.3745,
                    -10.4255, -10.4255, -13.3745, -13.3745, -10.4255, -10.4255,
                    -13.3745, -14.3745, -11.4255, -11.4255, -14.3745, -13.3745,
                    -10.4255, -10.4255, -13.3745, -13.3745, -10.4255, -10.4255,
                    -13.3745};

const double xFT0C[FT0C_CH_N] = {-10.30056, -10.30056, -7.35156, -7.35156, -4.42159,
                    -4.42159, -1.47259, -1.47259, 1.47259, 1.47259, 4.42159,
                    4.42159, 7.35156, 7.35156, 10.30056, 10.30056, -16.1339,
                    -16.1339, -13.1849, -13.1849, -10.30056, -10.30056,
                    -7.35156, -7.35156, -4.42159, -4.42159, -1.47259, -1.47259,
                    1.47259, 1.47259, 4.42159, 4.42159, 7.35156, 7.35156,
                    10.30056, 10.30056, 13.1849, 13.1849, 16.1339, 16.1339,
                    -16.1339, -16.1339, -13.1849, -13.1849, -10.30056,
                    -10.30056, -7.35156, -7.35156, 7.35156, 7.35156, 10.30056,
                    10.30056, 13.1849, 13.1849, 16.1339, 16.1339, -16.1339,
                    -16.1339, -13.1849, -13.1849, -10.30056, -10.30056,
                    -7.35156, -7.35156, 7.35156, 7.35156, 10.30056, 10.30056,
                    13.1849, 13.1849, 16.1339, 16.1339, -16.1339, -16.1339,
                    -13.1849, -13.1849, -10.30056, -10.30056, -7.35156,
                    -7.35156, -4.42159, -4.42159, -1.47259, -1.47259, 1.47259,
                    1.47259, 4.42159, 4.42159, 7.35156, 7.35156, 10.30056,
                    10.30056, 13.1849, 13.1849, 16.1339, 16.1339, -10.30056,
                    -10.30056, -7.35156, -7.35156, -4.42159, -4.42159,
                    -1.47259, -1.47259, 1.47259, 1.47259, 4.42159, 4.42159,
                    7.35156, 7.35156, 10.30056, 10.30056};

const double yFT0C[FT0C_CH_N] = {13.1849, 16.1339, 16.1339, 13.1849, 13.1849, 16.1339,
                    16.1339, 13.1849, 13.1849, 16.1339, 16.1339, 13.1849,
                    13.1849, 16.1339, 16.1339, 13.1849, 7.35156, 10.30056,
                    10.30056, 7.35156, 7.35156, 10.30056, 10.30056, 7.35156,
                    7.35156, 10.30056, 10.30056, 7.35156, 7.35156, 10.30056,
                    10.30056, 7.35156, 7.35156, 10.30056, 10.30056, 7.35156,
                    7.35156, 10.30056, 10.30056, 7.35156, 1.47259, 4.42159,
                    4.42159, 1.47259, 1.47259, 4.42159, 4.42159, 1.47259,
                    1.47259, 4.42159, 4.42159, 1.47259, 1.47259, 4.42159,
                    4.42159, 1.47259, -4.42159, -1.47259, -1.47259, -4.42159,
                    -4.42159, -1.47259, -1.47259, -4.42159, -4.42159, -1.47259,
                    -1.47259, -4.42159, -4.42159, -1.47259, -1.47259, -4.42159,
                    -10.30056, -7.35156, -7.35156, -10.30056, -10.30056,
                    -7.35156, -7.35156, -10.30056, -10.30056, -7.35156,
                    -7.35156, -10.30056, -10.30056, -7.35156, -7.35156,
                    -10.30056, -10.30056, -7.35156, -7.35156, -10.30056,
                    -10.30056, -7.35156, -7.35156, -10.30056, -16.1339,
                    -13.1849, -13.1849, -16.1339, -16.1339, -13.1849, -13.1849,
                    -16.1339, -16.1339, -13.1849, -13.1849, -16.1339, -16.1339,
                    -13.1849, -13.1849, -16.1339};


void AnalyseUsing3Sub(const char* inputFile = "o2sim.root") {
    TString inFileName(inputFile);
    TFile *fIn = new TFile(inFileName, "read");

    TString outFileName("o2analysisdata.root");
    TFile *fOut = new TFile(outFileName, "recreate");

    TTree *tree = (TTree*)fIn->Get("o2sim");

    vector<o2::ft0::HitType>* hitArrFT0 = nullptr;
    tree->SetBranchAddress("FT0Hit", &hitArrFT0);
    TLeaf* ft0ChID = tree->GetLeaf("FT0Hit.mDetectorID");

    vector<o2::fv0::Hit>* hitArrFV0 = nullptr;
    tree->SetBranchAddress("FV0Hit", &hitArrFV0);
    TLeaf* fv0ChID = tree->GetLeaf("FV0Hit.mDetectorID");

    //TH1D *hVobs[N];
    //TH1D *hRsub[N];
    TH1D *hSumAB[N];
    TH1D *hSumAC[N];
    TH1D *hSumBC[N];

    for (int i=0; i<N; i++) {
        //hVobs[i] = new TH1D(Form("hV%dobs", i+1), Form("hV%dobs", i+1), 100, -2.0, 2.0);
        //hRsub[i] = new TH1D(Form("hRsub%d", i+1), Form("hRsub%d", i+1), 100, -1.0, 1.0);
        hSumAB[i] = new TH1D(Form("hSumAB%d", i+1), Form("hSumAB%d", i+1), 100, -1.0, 1.0);
        hSumAC[i] = new TH1D(Form("hSumAC%d", i+1), Form("hSumAC%d", i+1), 100, -1.0, 1.0);
        hSumBC[i] = new TH1D(Form("hSumBC%d", i+1), Form("hSumBC%d", i+1), 100, -1.0, 1.0);
    }

    TComplex QvecA[N], QvecB[N], QvecC[N];
    double normA[N], normB[N], normC[N];
    vector<double> phiList;

    double vobs = 0.0, epA = 0.0, epB = 0.0, epC =0.0;
    double sumAB[N], sumAC[N], sumBC[N];
    int n = 0;

    Int_t nentries = (Int_t)tree->GetEntries();
    Int_t arrSize;

    int noutput = nentries/20;
    if (noutput<1) noutput = 1;
    for ( Int_t i=0; i<nentries; i++ ) {
        if ( i % noutput == 0 )
            cout << 100*(double)i/(double)nentries << " % finished" << endl;

        Int_t entry = tree->GetEntry(i);

        for ( Int_t j=0; j<N; j++ ) {
            n = j+1;

            normA[j] = 0; normB[j] = 0; normC[j] = 0;
            QvecA[j] = TComplex(0, 0); QvecB[j] = TComplex(0, 0); QvecC[j] = TComplex(0, 0);

            arrSize = hitArrFV0->size();
            for ( Int_t k=0; k<arrSize; k++ ) {

                const auto& fv0Hit = (*hitArrFV0)[k];
                double phi = GetFV0Phi(fv0ChID->GetValue(k));

                CalculateQvector(QvecA[j], normA[j], n, phi, 0, 0);
            }

            arrSize = hitArrFT0->size();
            for ( Int_t k=0; k<arrSize; k++ ) {

                const auto& ft0Hit = (*hitArrFT0)[k];
                double phi = GetFT0Phi(ft0ChID->GetValue(k));

                double z = ft0Hit.GetZ();

                if (z>0) {
                    CalculateQvector(QvecB[j], normB[j], n, phi, 0, 0);
                } else {
                    CalculateQvector(QvecC[j], normC[j], n, phi, 0, 0);
                }
            }

            epA = GetEventPlane(QvecA[j], n);
            epB = GetEventPlane(QvecB[j], n);
            epC = GetEventPlane(QvecC[j], n);

            hSumAB[j]->Fill(TMath::Cos(n*(epA - epB)));
            hSumAC[j]->Fill(TMath::Cos(n*(epA - epC)));
            hSumBC[j]->Fill(TMath::Cos(n*(epB - epC)));

            sumAB[j] += TMath::Cos(n*(epA - epB));
            sumAC[j] += TMath::Cos(n*(epA - epC));
            sumBC[j] += TMath::Cos(n*(epB - epC));
        }
    }

    for ( Int_t i=0; i<N; i++ ) {

        sumAB[i] /= nentries;
        sumAC[i] /= nentries;
        sumBC[i] /= nentries;

        double res = CalculateResolution(sumAB[i], sumAC[i], sumBC[i]);

        cout << "R" << i+1 << " : " << res << "\n";

        hSumAB[i]->Write(Form("hSumAB%d", i+1));
        hSumAC[i]->Write(Form("hSumAC%d", i+1));
        hSumBC[i]->Write(Form("hSumBC%d", i+1));
    }
}

//______________________________________________________________________
double GetFT0Eta(int chID, double z) {
    double x = GetFT0ChannelCenterX(chID);
    double y = GetFT0ChannelCenterY(chID);
    double theta = TMath::ATan2(TMath::Sqrt(x*x + y*y), z);
    return -TMath::Log(TMath::Tan(theta/2.0));
}

//______________________________________________________________________
double GetFV0Eta(int chID) {

    double eta[5] = {4.765, 4.185, 3.655, 3.11, 2.505};

    for (int i=1; i<6; i++) {
        if (chID < i*8) return eta[i-1];
    }

    return 0.0;
}

//______________________________________________________________________________
double GetFV0Phi(int chID) {
    return TMath::Pi()/8.0 + (chID%8)*TMath::Pi()/4.0;
}

//______________________________________________________________________
double GetFT0Phi(int chID) {
    double x = GetFT0ChannelCenterX(chID);
    double y = GetFT0ChannelCenterY(chID);
    return TMath::ATan2(y, x);
}

//______________________________________________________________________
// FT0-A: Channels 0-95
// FT0-C: Channels 96-207
double GetFT0ChannelCenterX(int chID) {
    if (chID<FT0A_CH_N) {
        return xFT0A[chID];
    } else {
        return xFT0A[chID-FT0A_CH_N-1];
    }
}

//______________________________________________________________________
double GetFT0ChannelCenterY(int chID) {
    if (chID<FT0A_CH_N) {
        return yFT0A[chID];
    } else {
        return yFT0A[chID-FT0A_CH_N-1];
    }
}

//______________________________________________________________________
void CalculateQvector(TComplex &Qvec, double &norm, int n, double phi, double pt, bool bUseWeight) {

    double w = 1.0;
    if (bUseWeight) w = pt;
    norm += w*w;

    Qvec += TComplex(w*TMath::Cos(n*phi), w*TMath::Sin(n*phi));
}

//______________________________________________________________________
double GetEventPlane(TComplex Qvec, int n) {
    return TMath::ATan2(Qvec.Im(), Qvec.Re())/n;
}

//______________________________________________________________________
double CalculateResolution(double sumAB, double sumAC, double sumBC) {
    return TMath::Sqrt((sumAB*sumAC)/sumBC);
}

//______________________________________________________________________
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

//______________________________________________________________________
double GetVnObs(TComplex Qvec, double phi, int n) {
    return TMath::Cos(n*(phi - GetEventPlane(Qvec, n)));
}

//______________________________________________________________________
double R1(double khi) {
    return (TMath::Sqrt(TMath::Pi())/2)*khi*TMath::Exp(-khi*khi/2)*(TMath::BesselI0(khi*khi/2) + TMath::BesselI1(khi*khi/2));
}

//______________________________________________________________________
double func(double *x, double *p) {
    double khi = x[0];
    double Rk = p[0];
    return (TMath::Sqrt(TMath::Pi())/2)*khi*TMath::Exp(-khi*khi/2)*(TMath::BesselI0(khi*khi/2) + TMath::BesselI1(khi*khi/2)) - Rk;
}

//______________________________________________________________________
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

//______________________________________________________________________
double CalculateRerror(double khi, double khiErr) {
    double bessel0 = -(khi*khi-2)*TMath::BesselI0(khi*khi/2);
    double bessel1 = 2*TMath::BesselI1(khi*khi/2);
    double bessel2 = khi*khi*TMath::BesselI(2,khi*khi/2);
    return TMath::Sqrt(TMath::Pi()/4)*TMath::Exp(-khi*khi)*(bessel0 + bessel1 + bessel2)*khiErr;
}
