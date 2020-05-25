#include <vector>
#include <cstring>

#include "Pythia8/Pythia.h"
#include "Pythia8/HeavyIons.h"
#include "Pythia8/HIUserHooks.h"

#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TNtuple.h"
#include "TH1D.h"

using namespace Pythia8;

double GetPhi0(double phi, double *vn, double *psi);
double AnisotropicPhiDist(double *x, double *p);
double GetAnisotropicPhi(double phi0, double phiInit, double err, double *vn, double *psi, TF1 *fPhiDist);
double BisectionMethod(double phi0, double err, double *vn, double *psi, TF1 *fPhiDist);
int GetCentralityBin(double b);

#define NCENTBINS 10
#define NCOEF 5

int main(int argc, char *argv[]) {

    if (argc==1) {
        cout << "Usage : ./pythiaHI 'number of events'[=100] 'output file name'[=output.root] 'seed'[=0] 'save impact parameter'[=0] 'use centrality dependent vn'[=0] 'give bMin'[=0] 'give bMax'[=20.0]" << endl;
        return 0;
    }

    int nEvents = argc > 1 ? atol(argv[1]) : 100;
    TString outFileName = argc > 2 ? argv[2] : "output.root";
    int seed = argc > 3 ? atol(argv[3]) : 0;
    bool bSaveb = argc > 4 ? atol(argv[3]) : 0;
    bool bCentDep = argc > 5 ? atol(argv[5]) : 0;
    double bMin = argc > 6 ? atof(argv[6]) : 0.0;
    double bMax = argc > 7 ? atof(argv[7]) : 20.0;

    TFile *fOut = new TFile(outFileName, "RECREATE");

    Pythia pythia;

    // Initialise pythia
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed = %d", seed));
    pythia.readString("Beams:idA = 1000822080");
    pythia.readString("Beams:idB = 1000822080"); // The lead ion.
    //pythia.readString("ParticleDecays:limitTau0 = On");
    //pythia.readString("ParticleDecays:tau0Max = 1.0");

    // Settings 1
    //pythia.readString("Beams:eCM = 2760.0");
    //pythia.readString("HeavyIon:SigFitErr = 0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0");
    //pythia.readString("HeavyIon:SigFitDefPar = 17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0");
    //pythia.readString("HeavyIon:SigFitNGen = 20");

    // Settings 2
    pythia.readString("Beams:eCM = 5520.0");
    pythia.readString("HeavyIon:SigFitNGen 0");
    pythia.readString("HeavyIon:SigFitDefPar 14.82,1.82,0.25,0.0,0.0,0.0,0.0,0.0");

    pythia.init();

    std::vector<int> moms;

    // Flow related stuff

    //double vn[NCOEF] = {0., 0.15, 0.08, 0.03, 0.01};
    double vn[NCOEF] = {0., 0.1, 0.0, 0.0, 0.0};

    if (bCentDep==0) {
        cout << "vn : [ ";
        for (int i=0; i<NCOEF; i++) {
            cout << vn[i] << " ";
            if (i==NCOEF-1) cout << "]\n";
        }
    }

    // Data from arXiv:1804.02944 (10.17182/hepdata.83737)
    // These are just rough estimates
    double centvn[NCOEF][NCENTBINS] = {
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.035, 0.063, 0.082, 0.093, 0.097, 0.095, 0.086, 0.073, 0.060, 0.050},
        {0.021, 0.026, 0.028, 0.029, 0.029, 0.027, 0.022, 0.020, 0.018, 0.015},
        {0.010, 0.012, 0.013, 0.014, 0.015, 0.015, 0.014, 0.013, 0.012, 0.010},
        {0.0045, 0.0053, 0.0062, 0.0070, 0.0070, 0.0065, 0.0060, 0.0055, 0.0050, 0.0040}
    };

    double psi[5] = {TMath::Pi(), TMath::Pi(), TMath::Pi(), TMath::Pi(), TMath::Pi()};
    //double psi[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    TRandom3 *rand = new TRandom3(seed);
    TF1 *fPhiDist = new TF1("fPhiDist", AnisotropicPhiDist, -TMath::Pi(), TMath::Pi(), 11);

    // For testing purposes
    TH1D *hPhi = new TH1D("hPhi", "hPhi", 50, -TMath::Pi(), TMath::Pi());
    TH1D *hV2 = new TH1D("hV2", "hV2", 100, -0.2, 0.4);
    TH1D *hPhiDiff = new TH1D("hPhiDiff", "hPhiDiff", 321, -TMath::Pi()/2.0, TMath::Pi()/2.0);

    // NTuple to save events
    TNtuple *ntuple;
    if (bSaveb) {
        ntuple = new TNtuple("pythiaEvents", "data from Pythia8 with afterburner", "eventId:particleId:px:py:pz:x:y:z:isHadron:charge:b");
    } else {
        ntuple = new TNtuple("pythiaEvents", "data from Pythia8 with afterburner", "eventId:particleId:px:py:pz:x:y:z:isHadron:charge");
    }

    // Track variables
    Int_t pid = 0;
    Float_t px = 0.0, py = 0.0, pz = 0.0, m = 0.0, x = 0.0, y = 0.0, z = 0.0, t = 0.0;
    Double_t charge = 0.0;
    bool isHadron = false;

    // Loop over events.
    for ( int iEvent = 0; iEvent < nEvents; ++iEvent ) {

        if ( !pythia.next() ) continue;

        for (int j=0; j<5; j++) {
            psi[j] = rand->Uniform(-TMath::Pi(), TMath::Pi());
        }

        //pythia.event.list(false, true);

        double b = pythia.info.hiinfo->b();

	    if (b<bMin || b>bMax) {
            iEvent--;
            continue;
	    }

        int ibin = 0;
        if (bCentDep) {
            ibin = GetCentralityBin(b);
            if (ibin==-1) ibin = NCENTBINS-1;
            for (int i=0; i<NCOEF; i++) {
                vn[i] = centvn[i][ibin];
            }
        }

        double v2 = 0.0;
        int mult = 0;
        for (int iPart = 0; iPart < pythia.event.size(); iPart++) {

            if (pythia.event[iPart].isFinal() && (pythia.event[iPart].name()!="NucRem")) {

                charge = pythia.event[iPart].charge();
                if (charge==0) continue;

                int partIndex = iPart;
                while (partIndex!=0) {
                    if (!pythia.event.at(pythia.event.at(partIndex).mother1()).isHadron()) break; // Check that the mother is hadron
                    partIndex = pythia.event.at(partIndex).mother1();
                }

                double phi0 = pythia.event.at(partIndex).phi();
                double phi = GetAnisotropicPhi(phi0, phi0, 0.0001, vn, psi, fPhiDist);

                double phiDiff = phi - phi0;
                pythia.event[iPart].rot(0, phiDiff);

                pid = pythia.event[iPart].id();

                px = pythia.event[iPart].px();
                py = pythia.event[iPart].py();
                pz = pythia.event[iPart].pz();

                x = pythia.event[iPart].xProd();
                y = pythia.event[iPart].yProd();
                z = pythia.event[iPart].zProd();

                isHadron = pythia.event[iPart].isHadron();

                if (bSaveb) {
                    ntuple->Fill(iEvent, pid, px, py, pz, x, y, z, isHadron, charge, b);
                } else {
                    ntuple->Fill(iEvent, pid, px, py, pz, x, y, z, isHadron, charge);
                }

                hPhi->Fill(pythia.event[iPart].phi());
                hPhiDiff->Fill(phiDiff);
                v2 += TMath::Cos(2.0 * (pythia.event[iPart].phi() - psi[1]));
                mult++;
            }
        }

        hV2->Fill(v2/mult);

        cout << "\nEvent " << iEvent + 1
             << " done, n=" << pythia.event.size() << endl;
        cout << "Impact parameter : " << pythia.info.hiinfo->b() << endl;
        cout << "vn : [ ";
        for (int i=0; i<NCOEF; i++) {
            cout << vn[i] << " ";
        }
        cout << "]\n";
    }

    ntuple->Write("", TObject::kOverwrite);
    hPhi->Write("hPhi");
    hV2->Write("hV2");
    hPhiDiff->Write("hPhiDiff");
    fOut->Close();

    pythia.stat();

    return 0;
}

//______________________________________________________________________________
//              afterburner functions
//______________________________________________________________________________

double GetPhi0(double phi, double *vn, double *psi) {
    return phi + 2.0*vn[0]*TMath::Sin(phi-psi[0]) + vn[1]*TMath::Sin(2.0*(phi-psi[1])) + (2.0/3.0)*vn[2]*TMath::Sin(3.0*(phi-psi[2])) + (1.0/2.0)*vn[3]*TMath::Sin(4.0*(phi-psi[3])) + (2.0/5.0)*vn[4]*TMath::Sin(5.0*(phi-psi[4]));
}

double AnisotropicPhiDist(double *x, double *p) {
    double phi = x[0];
    double phi0 = p[0];
    double v1 = p[1];
    double v2 = p[2];
    double v3 = p[3];
    double v4 = p[4];
    double v5 = p[5];
    double psi1 = p[6];
    double psi2 = p[7];
    double psi3 = p[8];
    double psi4 = p[9];
    double psi5 = p[10];
    return phi - phi0 + 2.0*v1*TMath::Sin(phi-psi1) + v2*TMath::Sin(2.0*(phi-psi2)) + (2.0/3.0)*v3*TMath::Sin(3.0*(phi-psi3)) + (1.0/2.0)*v4*TMath::Sin(4.0*(phi-psi4)) + (2.0/5.0)*v5*TMath::Sin(5.0*(phi-psi5));
}


double GetAnisotropicPhi(double phi0, double phiInit, double err, double *vn, double *psi, TF1 *fPhiDist) {

    double phi = 0;
    fPhiDist->SetParameters(phi0, vn[0], vn[1], vn[2], vn[3], vn[4], psi[0], psi[1], psi[2], psi[3], psi[4]);

    int step = 0;
    while (TMath::Abs(GetPhi0(phi, vn, psi) - phi0) > err) {
        phi = phiInit - fPhiDist->Eval(phiInit)/fPhiDist->Derivative(phiInit);
        phiInit = phi;

        step++;

        if (step==100) {
            cout << "Newton not converging, switch to bisection method" << endl;
            phi = BisectionMethod(phi0, err, vn, psi, fPhiDist);
            break;
        }
    }

    if (phi>TMath::Pi()) { phi -= 2*TMath::Pi(); }
    else if (phi<-TMath::Pi()) { phi += 2*TMath::Pi(); }

    return phi;
}

double BisectionMethod(double phi0, double err, double *vn, double *psi, TF1 *fPhiDist) {

    fPhiDist->SetParameters(phi0, vn[0], vn[1], vn[2], vn[3], vn[4], psi[0], psi[1], psi[2], psi[3], psi[4]);
    double a = -TMath::Pi(), b = TMath::Pi(), c = 0;
    while (TMath::Abs(b-a) > err) {
        c = (a+b)/2.0;
        if (fPhiDist->Eval(c) < 0) {
            a = c;
        } else if (fPhiDist->Eval(c) > 0) {
            b = c;
        } else {
            break;
        }
    }

    return c;
}

int GetCentralityBin(double b) {
    double bins[NCENTBINS+1] = {0.0, 4.94, 6.98, 8.55, 9.88, 11.04, 12.09,
                                13.05, 13.97, 14.96, 20.0};

    for (int i=0; i<NCENTBINS-1; i++) {
        if ((b>bins[i-1]) && (b<=bins[i])) return i;
    }

    return -1;
}
