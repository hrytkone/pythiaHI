/*
 * Iteration using Newton's method
 */

#include "TMath.h"
#include <boost/math/special_functions/bessel.hpp>

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
