#include <fstream>
#include <iostream>
#include <chrono>  // std::chrono::system_clock
#include <random>  // std::default_random_engine
#include <cmath>
#include <vector>
#include <algorithm>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TPolyLine.h"
#include "TRandom3.h"
#include "TString.h"
#include "TTree.h"
#include "TVector2.h"

#include "Particle.h"

using namespace std;
const float mass_pion = 0.13957039;
const float mass_proton = 0.93827208816;
const bool shuffle_mode = false;

float GetDistence(TVector2 p0, TVector2 p1) {
  float dis = (p0 - p1).Mod();
  return dis;
}
float GetDistence(TVector2 p0, TVector2 p1, TVector2 p2) {
  float dis = 1./3.* ((p0-p1).Mod() + (p1-p2).Mod() + (p2-p0).Mod());
  return dis;
}

Parton GetLBCFriendForThisParton(Parton parton, float sigmaLBC ,float bn_friend) {
  Parton partonLBC;
  partonLBC.SetBaryonNumber(bn_friend);
  float r = gRandom->Uniform(0,1);
  float phi = parton.GetPositionVector().Phi();
  float dphi = gRandom->Gaus(0, sigmaLBC);
  phi += dphi;
  partonLBC.SetPositionVector(sqrt(r) * cos(phi), sqrt(r) * sin(phi));
  partonLBC.SetSerialNumber(-parton.GetSerialNumber());
  return partonLBC;
}

const int MaxNumTracks = 100000;
typedef struct {
  unsigned int nevent;
  unsigned int multi;
  int fBaryonNumber[MaxNumTracks];
  float fHadronX[MaxNumTracks];
  float fHadronY[MaxNumTracks];
  float fHadronPx[MaxNumTracks];
  float fHadronPy[MaxNumTracks];
  float fDistance[MaxNumTracks];
  float fParton0BaryonNumber[MaxNumTracks];
  float fParton0X[MaxNumTracks];
  float fParton0Y[MaxNumTracks];
  float fParton1BaryonNumber[MaxNumTracks];
  float fParton1X[MaxNumTracks];
  float fParton1Y[MaxNumTracks];
  float fParton2BaryonNumber[MaxNumTracks];
  float fParton2X[MaxNumTracks];
  float fParton2Y[MaxNumTracks];
} Cell_t;

int main(int argc, char** argv) {
  unsigned int nEvents = 10000;
  unsigned int nPartons = 30;
  unsigned int nAntiPartons = 30;
  float rM = 1.5;
  float rho0 = 0.;
  float rho2 = 0.;
  float frac_qq = 1;
  float frac_qqbar = 1;
  float sigmaLBC = 0.1;
  
  char* FileInput = 0;
  char* FileOutput = 0;
  if (argc != 3 && argc != 1) return 0;
  if (argc == 3) {
    FileInput = argv[1];
    FileOutput = argv[2];
  }

  ifstream filein;
  filein.open(FileInput);
  if (!filein) return (0);

  while (filein.good()) {
    filein >> nEvents >> nPartons >> nAntiPartons >> rM >> rho0 >> rho2 >> frac_qq >> frac_qqbar >> sigmaLBC;
  }
  filein.close();
  cout << "nEvents      = " << nEvents << endl;
  cout << "nPartons     = " << nPartons << endl;
  cout << "nAntiPartons = " << nAntiPartons << endl;
  cout << "rM           = " << rM << endl;
  cout << "rho0         = " << rho0 << endl;
  cout << "rho2         = " << rho2 << endl;
  cout << "frac_qq      = " << rho2 << endl;
  cout << "frac_qqbar   = " << rho2 << endl;
  cout << "sigmaLBC     = " << sigmaLBC << endl;
  sigmaLBC *= TMath::Pi();
  
  // Generate events
  TTree* tree = new TTree("CoalData", "CoalData DST Tree");
  Cell_t cell;
  tree->Branch("Event", &cell.nevent, "nevent/I:multi/I");
  tree->Branch("BaryonNum", cell.fBaryonNumber, "BaryonNum[multi]/I");
  tree->Branch("HadronX", cell.fHadronX, "HadronX[multi]/F");
  tree->Branch("HadronY", cell.fHadronY, "HadronY[multi]/F");
  tree->Branch("HadronPx", cell.fHadronPx, "HadronPx[multi]/F");
  tree->Branch("HadronPy", cell.fHadronPy, "HadronPy[multi]/F");
  tree->Branch("Distance", cell.fDistance, "Distance[multi]/F");
  tree->Branch("Parton0BaryonNum", cell.fParton0BaryonNumber, "Parton0BaryonNumber[multi]/F");
  tree->Branch("Parton0X", cell.fParton0X, "Parton0X[multi]/F");
  tree->Branch("Parton0Y", cell.fParton0Y, "Parton0Y[multi]/F");
  tree->Branch("Parton1BaryonNum", cell.fParton1BaryonNumber, "Parton1BaryonNumber[multi]/F");
  tree->Branch("Parton1X", cell.fParton1X, "Parton1X[multi]/F");
  tree->Branch("Parton1Y", cell.fParton1Y, "Parton1Y[multi]/F");
  tree->Branch("Parton2BaryonNum", cell.fParton2BaryonNumber, "Parton2BaryonNumber[multi]/F");
  tree->Branch("Parton2X", cell.fParton2X, "Parton2X[multi]/F");
  tree->Branch("Parton2Y", cell.fParton2Y, "Parton2Y[multi]/F");

  unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
  TRandom3* rndm = new TRandom3(seed);

  for (size_t iEvent = 0; iEvent < nEvents; iEvent++) {
    if (iEvent % 1000 == 0) std::cout << "Processing event # " << iEvent << endl;
    // Seed Partons
    vector<Parton> vecPartons;
    vector<Parton> vecAntiPartonsTemp;
    unsigned int n = 0;
    for (size_t i = 0; i < nPartons; i++) {
      Parton parton;
      parton.SetBaryonNumber(1. / 3.);
      float r = rndm->Uniform(0, 1);
      float phi = rndm->Uniform(0, TMath::TwoPi());
      parton.SetPositionVector(sqrt(r) * cos(phi), sqrt(r) * sin(phi));
      parton.SetSerialNumber(n);
      vecPartons.push_back(parton);
      if (rndm->Uniform(0,1) < frac_qq)    vecPartons.push_back(GetLBCFriendForThisParton(parton, sigmaLBC,  1./3.));
      if (rndm->Uniform(0,1) < frac_qqbar) vecPartons.push_back(GetLBCFriendForThisParton(parton, sigmaLBC, -1./3.));
      n++;
    }
    for (size_t i = 0; i < nAntiPartons; i++) {
      Parton antiparton;
      antiparton.SetBaryonNumber(-1. / 3.);
      float r = rndm->Uniform(0, 1);
      float phi = rndm->Uniform(0, TMath::TwoPi());
      antiparton.SetPositionVector(sqrt(r) * cos(phi), sqrt(r) * sin(phi));
      antiparton.SetSerialNumber(n);
      vecAntiPartonsTemp.push_back(antiparton);
      if (rndm->Uniform(0,1) < frac_qq)    vecAntiPartonsTemp.push_back(GetLBCFriendForThisParton(antiparton, sigmaLBC, -1./3.));
      if (rndm->Uniform(0,1) < frac_qqbar) vecAntiPartonsTemp.push_back(GetLBCFriendForThisParton(antiparton, sigmaLBC,  1./3.));
      n++;
    }
    shuffle(vecPartons.begin(),vecPartons.end(),std::default_random_engine(seed+1));
    shuffle(vecAntiPartonsTemp.begin(),vecAntiPartonsTemp.end(),std::default_random_engine(seed+2));
    vecPartons.resize(nPartons);
    vecAntiPartonsTemp.resize(nAntiPartons);
    vecPartons.insert(vecPartons.end(), vecAntiPartonsTemp.begin(), vecAntiPartonsTemp.end());

    vector<HardonCandidate> vecHardonCandidates;
    // Get Meson Candidates
    unsigned int serial_number = 0;
    for (size_t i = 0; i < vecPartons.size(); i++) {
      for (size_t j = i + 1; j < vecPartons.size(); j++) {
        bool ismeson = fabs(vecPartons[i].GetBaryonNumber() + vecPartons[j].GetBaryonNumber()) < 1.e-6;
        if (ismeson) {
          float distance = rM * GetDistence(vecPartons[i].GetPositionVector(), vecPartons[j].GetPositionVector());
          float x = 0.5 * (vecPartons[i].GetPositionVector() + vecPartons[j].GetPositionVector()).X();
          float y = 0.5 * (vecPartons[i].GetPositionVector() + vecPartons[j].GetPositionVector()).Y();

          vector<unsigned int> vec_parton_sn;
          vec_parton_sn.push_back(vecPartons[i].GetSerialNumber());
          vec_parton_sn.push_back(vecPartons[j].GetSerialNumber());
          vector<TVector2> vec_parton_pv;
          vec_parton_pv.push_back(vecPartons[i].GetPositionVector());
          vec_parton_pv.push_back(vecPartons[j].GetPositionVector());
          vector<float> vec_parton_bn;
          vec_parton_bn.push_back(vecPartons[i].GetBaryonNumber());
          vec_parton_bn.push_back(vecPartons[j].GetBaryonNumber());

          HardonCandidate meson_candidate;
          meson_candidate.SetBaryonNumber(0);
          meson_candidate.SetSerialNumber(serial_number);
          meson_candidate.SetVecPartonsSerialNumber(vec_parton_sn);
          meson_candidate.SetPositionVector(x, y);
          meson_candidate.SetVecPartonsPosition(vec_parton_pv);
          meson_candidate.SetMeanDistance(distance);
          meson_candidate.SetVecPartonsBaryonNumber(vec_parton_bn);
          vecHardonCandidates.push_back(meson_candidate);
          serial_number++;
        }
      }
    }

    // Get Baryon Candidates
    for (size_t i = 0; i < (size_t)vecPartons.size(); i++) {
      for (size_t j = i + 1; j < (size_t)vecPartons.size(); j++) {
        for (size_t k = j + 1; k < (size_t)vecPartons.size(); k++) {
          bool isbayron = fabs(vecPartons[i].GetBaryonNumber() + vecPartons[j].GetBaryonNumber() + vecPartons[k].GetBaryonNumber() - 1.) < 1.e-6;
          bool isantibaryon = fabs(vecPartons[i].GetBaryonNumber() + vecPartons[j].GetBaryonNumber() + vecPartons[k].GetBaryonNumber() + 1.) < 1.e-6;
          if (isbayron || isantibaryon) {
            float distance = GetDistence(vecPartons[i].GetPositionVector(), vecPartons[j].GetPositionVector(), vecPartons[k].GetPositionVector());
            float x = 1. / 3. * (vecPartons[i].GetPositionVector() + vecPartons[j].GetPositionVector() + vecPartons[k].GetPositionVector()).X();
            float y = 1. / 3. * (vecPartons[i].GetPositionVector() + vecPartons[j].GetPositionVector() + vecPartons[k].GetPositionVector()).Y();
            vector<unsigned int> vec_parton_sn;
            vec_parton_sn.push_back(vecPartons[i].GetSerialNumber());
            vec_parton_sn.push_back(vecPartons[j].GetSerialNumber());
            vec_parton_sn.push_back(vecPartons[k].GetSerialNumber());
            vector<TVector2> vec_parton_pv;
            vec_parton_pv.push_back(vecPartons[i].GetPositionVector());
            vec_parton_pv.push_back(vecPartons[j].GetPositionVector());
            vec_parton_pv.push_back(vecPartons[k].GetPositionVector());
            vector<float> vec_parton_bn;
            vec_parton_bn.push_back(vecPartons[i].GetBaryonNumber());
            vec_parton_bn.push_back(vecPartons[j].GetBaryonNumber());
            vec_parton_bn.push_back(vecPartons[k].GetBaryonNumber());

            HardonCandidate baryon_candidate;
            if (isbayron) baryon_candidate.SetBaryonNumber(1);
            if (isantibaryon) baryon_candidate.SetBaryonNumber(-1);
            baryon_candidate.SetSerialNumber(serial_number);
            baryon_candidate.SetVecPartonsSerialNumber(vec_parton_sn);
            baryon_candidate.SetPositionVector(x, y);
            baryon_candidate.SetVecPartonsPosition(vec_parton_pv);
            baryon_candidate.SetVecPartonsBaryonNumber(vec_parton_bn);
            baryon_candidate.SetMeanDistance(distance);
            vecHardonCandidates.push_back(baryon_candidate);
            serial_number++;
          }
        }
      }
    }

    if(!shuffle_mode) {
    // sort the hardon candidates by the distance
    sort(vecHardonCandidates.begin(), vecHardonCandidates.end());
    } else {
      // shuffle the hardon candidates
      shuffle(vecHardonCandidates.begin(),vecHardonCandidates.end(),std::default_random_engine(seed));
    }

    // to record the partons have been used
    vector<unsigned int> vecNodedBuffer;
    // to record the Hardons
    vector<Hardon> vecHardons;

    // for every hardon candidate, if it's partons have not been used, it should be freezed out
    for (auto hardon_candidates : vecHardonCandidates) {
      // assume no partons for this hardon candidate has been used
      bool is_every_parton_free = true;
      // loop it's partons. If a parton is found in the vecNodeBuffer, no hardon will form
      for (auto sm : hardon_candidates.GetVecPartonsSerialNumber()) {
        vector<unsigned int>::iterator it = find(vecNodedBuffer.begin(), vecNodedBuffer.end(), sm);
        if (it != vecNodedBuffer.end()) {
          is_every_parton_free = false;
          break;
        }
      }
      // If no parton is found in the vecNodeBuffer, a hardon will be formed
      if (is_every_parton_free) {
        for (auto sm : hardon_candidates.GetVecPartonsSerialNumber()) {
          vector<unsigned int>::iterator it = find(vecNodedBuffer.begin(), vecNodedBuffer.end(), sm);
          // record the partons have been used in this coalesence
          if (it == vecNodedBuffer.end()) vecNodedBuffer.push_back(sm);
        }
        // package hardons
        Hardon hardon;
        hardon = hardon_candidates;
        vecHardons.push_back(hardon);
      }
    }  // coalescence end, we have get all the hardons

    cell.nevent = iEvent;
    cell.multi = vecHardons.size();
    unsigned int i = 0;
    for (auto hardon : vecHardons) {
      cell.fBaryonNumber[i] = hardon.GetBaryonNumber();
      cell.fHadronX[i] = hardon.GetPositionVector().X();
      cell.fHadronY[i] = hardon.GetPositionVector().Y();
      cell.fDistance[i] = hardon.GetDistance();
      cell.fParton0BaryonNumber[i] = hardon.GetVecPartonsBaryonNumber()[0];
      cell.fParton0X[i] = hardon.GetVecPartonsPosition()[0].X();
      cell.fParton0Y[i] = hardon.GetVecPartonsPosition()[0].Y();
      cell.fParton1BaryonNumber[i] = hardon.GetVecPartonsBaryonNumber()[1];
      cell.fParton1X[i] = hardon.GetVecPartonsPosition()[1].X();
      cell.fParton1Y[i] = hardon.GetVecPartonsPosition()[1].Y();
      if (hardon.GetBaryonNumber() != 0) {
        cell.fParton2BaryonNumber[i] = hardon.GetVecPartonsBaryonNumber()[2];
        cell.fParton2X[i] = hardon.GetVecPartonsPosition()[2].X();
        cell.fParton2Y[i] = hardon.GetVecPartonsPosition()[2].Y();
        hardon.SampleRawMomentum(mass_proton);
        hardon.MomentumBoost(rho0,rho2);
        cell.fHadronPx[i] = hardon.GetP().Px();
        cell.fHadronPy[i] = hardon.GetP().Py();
      } else {
        cell.fParton2BaryonNumber[i] = -999;
        cell.fParton2X[i] = -999;
        cell.fParton2Y[i] = -999;
        hardon.SampleRawMomentum(mass_pion);
        hardon.MomentumBoost(rho0,rho2);
        cell.fHadronPx[i] = hardon.GetP().Px();
        cell.fHadronPy[i] = hardon.GetP().Py();
      }
      i++;
    }
    tree->Fill();
  }  // generate events end;
  
  delete rndm;
  rndm = nullptr;

  TFile* file = new TFile(FileOutput, "RECREATE");
  file->cd();
  tree->Write();
  file->Close();
}  // main end
