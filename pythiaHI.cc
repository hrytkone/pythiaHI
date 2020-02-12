#include <vector>

#include "Pythia8/Pythia.h"
#include "Pythia8/HeavyIons.h"

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

int main() {

    int seed = 313;

    TString outFileName("pythia-output.root");
    TFile *fOut = new TFile(outFileName, "RECREATE");

    Pythia pythia;

    // Initialise pythia
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed = %d", seed));
    pythia.readString("Beams:idA = 1000822080");
    pythia.readString("Beams:idB = 1000822080"); // The lead ion.
    pythia.readString("Beams:eCM = 2760.0");
    pythia.readString("Beams:frameType = 1");
    pythia.readString("HeavyIon:SigFitErr = 0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0");
    pythia.readString("HeavyIon:SigFitDefPar = 17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0");
    pythia.readString("HeavyIon:SigFitNGen = 20");
    //pythia.readString("ParticleDecays:limitTau0 = On");
    //pythia.readString("ParticleDecays:tau0Max = 10.0");
    pythia.readString("HeavyIon:bWidth = 0.0");
    pythia.init();

    std::vector<int> moms;

    // Flow related stuff
    double vn[5] = {0., 0.15, 0.08, 0.03, 0.01};
    double psi[5] = {TMath::Pi(), TMath::Pi(), TMath::Pi(), TMath::Pi(), TMath::Pi()};
    TRandom3 *rand = new TRandom3(seed);
    TF1 *fPhiDist = new TF1("fPhiDist", AnisotropicPhiDist, -TMath::Pi(), TMath::Pi(), 11);

    TH1D *hPhi0 = new TH1D("hPhi0", "hPhi0", 100, -TMath::Pi(), TMath::Pi());
    TH1D *hPhi = new TH1D("hPhi", "hPhi", 100, -TMath::Pi(), TMath::Pi());

    // NTuple to save events
    TNtuple *ntuple = new TNtuple("pythiaEvents", "data from Pythia8 with afterburner", "eventId:particleId:px:py:pz:x:y:z:isHadron");

    // Track variables
    Int_t pid = 0;
    Float_t px = 0.0, py = 0.0, pz = 0.0, m = 0.0, x = 0.0, y = 0.0, z = 0.0, t = 0.0;
    Double_t charge = 0.0;
    bool isHadron = false;

    // Loop over events.
    int nEvents = 10;
    for ( int iEvent = 0; iEvent < nEvents; ++iEvent ) {

        if ( !pythia.next() ) continue;
        cout << "\nEvent " << iEvent + 1 << endl;

        for (int j=0; j<5; j++) {
            psi[j] = rand->Uniform(-TMath::Pi(), TMath::Pi());
        }

        //pythia.event.list();

        for (int iPart = 0; iPart < pythia.event.size(); iPart++) {

            if (pythia.event[iPart].isFinal()) {

                moms = pythia.event[iPart].motherList();
                double phi0 = pythia.event[iPart].phi();

                // Change particle angle according to flow.
                // Check that if the particle is from long-lived unstable
                // particle.
                if (moms.size()!=0) {
                    bool bTau0 = 0;
                    for (int i=0; i<moms.size(); i++) {
                        if (pythia.event.at(moms[i]).tau0()>10.0) bTau0 = 1;
                    }

                    if (bTau0) {

                    }
                }

                double phi = GetAnisotropicPhi(phi0, 1., 0.001, vn, psi, fPhiDist);
                double phiDiff = phi - phi0;
                pythia.event[iPart].rot(0, phiDiff);

                px = pythia.event[iPart].px();
                py = pythia.event[iPart].py();
                pz = pythia.event[iPart].pz();

                x = pythia.event[iPart].xProd();
                y = pythia.event[iPart].yProd();
                z = pythia.event[iPart].zProd();

                charge = pythia.event[iPart].charge();
                isHadron = pythia.event[iPart].isHadron();

                ntuple->Fill(iEvent, pid, px, py, pz, x, y, z, isHadron, charge);

                hPhi0->Fill(phi0);
                hPhi->Fill(phi);



            }
        }
    }

    ntuple->Write("", TObject::kOverwrite);
    hPhi0->Write("hPhi0");
    hPhi->Write("hPhi");
    fOut->Close();

    pythia.stat();

    return 0;
}

//______________________________________________________________________________
//              afterburner functions
//______________________________________________________________________________

double GetPhi0(double phi, double *vn, double *psi) {
    return phi - 2.0*vn[0]*TMath::Sin(phi-psi[0]) + vn[1]*TMath::Sin(2.0*(phi-psi[1])) + (2.0/3.0)*vn[2]*TMath::Sin(3.0*(phi-psi[2])) + (1.0/2.0)*vn[3]*TMath::Sin(4.0*(phi-psi[3])) + (2.0/5.0)*vn[4]*TMath::Sin(5.0*(phi-psi[4]));
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

    while (TMath::Abs(GetPhi0(phi, vn, psi) - phi0) > err) {
        phi = phiInit - fPhiDist->Eval(phiInit)/fPhiDist->Derivative(phiInit);
        phiInit = phi;
    }

    if (phi>TMath::Pi())phi -= 2*TMath::Pi();
    if (phi<-TMath::Pi()) phi += 2*TMath::Pi();

    return phi;
}
