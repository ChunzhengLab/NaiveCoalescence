#include "TRandom3.h"
#include "Particle.h"

void Hardon::MomentumBoost(float rho0, float rho2) {
  float etas = 0.;
  float phib = position_vector.Phi();
  // float rhob = pow(RadSys(rads,phis,Rx),1)*(rho0+rho2*cos(2.*phib)+rho3*cos(3.*phib)+rho4*cos(4.*phib));
  float rhob = position_vector.Mod() * (rho0 + rho2*cos(2*phib));
  TLorentzVector u;
  u.SetXYZT(sinh(rhob)*cos(phib),sinh(rhob)*sin(phib),cosh(rhob)*sinh(etas),cosh(rhob)*cosh(etas));
  TVector3 boostvector= u.BoostVector();
  p.Boost(boostvector);
}

void Hardon::SampleRawMomentum(float mass) {
  float pT = 1;
  float phi = gRandom->Uniform(-TMath::Pi(),TMath::Pi());
  p.SetPxPyPzE(pT*cos(phi),pT*sin(phi),0,pT+mass);
}
