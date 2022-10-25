#include <algorithm>
#include <chrono>  // std::chrono::system_clock
#include <cmath>
#include <iostream>
#include <random>  // std::default_random_engine
#include <vector>

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
const int nPartons = 60;
const int nAntiPartons = 60;
const float rM = 1;

// class Parton {
//  private:
//   unsigned int serial_number;
//   float baryon_number;
//   TVector2 position_vector;

//  public:
//   Parton() {
//     serial_number = 9999;
//     baryon_number = -999;
//     position_vector.SetX(-999);
//     position_vector.SetY(-999);
//   };
//   void SetSerialNumber(unsigned int n) { serial_number = n; }
//   void SetBaryonNumber(float bn) { baryon_number = bn; }
//   void SetPositionVector(float x, float y) { position_vector.Set(x, y); }
//   void SetPositionVector(TVector2 pv) { position_vector = pv; }

//   unsigned int GetSerialNumber() { return serial_number; }
//   float GetBaryonNumber() { return baryon_number; }
//   TVector2 GetPositionVector() { return position_vector; }
// };

// class Hardon {
//  protected:
//   unsigned int serial_number;
//   int baryon_number;
//   TVector2 position_vector;
//   float mean_distance;
//   vector<unsigned int> vec_partons_serial_number;
//   vector<float> vec_partons_baryon_number;
//   vector<TVector2> vec_partons_position;

//  public:
//   Hardon() {
//     serial_number = -999;
//     baryon_number = 9999;
//     position_vector.SetX(-999);
//     position_vector.SetY(-999);
//     mean_distance = -999;
//     vector<unsigned int>().swap(vec_partons_serial_number);
//     vector<float>().swap(vec_partons_baryon_number);
//     vector<TVector2>().swap(vec_partons_position);
//   }
//   void SetSerialNumber(unsigned int n) { serial_number = n; }
//   void SetBaryonNumber(int bn) { baryon_number = bn; }
//   void SetPositionVector(float x, float y) { position_vector.Set(x, y); }
//   void SetPositionVector(TVector2 pv) { position_vector = pv; }
//   void SetMeanDistance(float dis) { mean_distance = dis; }
//   void SetVecPartonsSerialNumber(vector<unsigned int> sm) { vec_partons_serial_number.assign(sm.begin(), sm.end()); }
//   void SetVecPartonsBaryonNumber(vector<float> bn) { vec_partons_baryon_number.assign(bn.begin(), bn.end()); }
//   void SetVecPartonsPosition(vector<TVector2> vecpv) { vec_partons_position.assign(vecpv.begin(), vecpv.end()); }

//   unsigned int GetSerialNumber() { return serial_number; }
//   int GetBaryonNumber() { return baryon_number; }
//   TVector2 GetPositionVector() { return position_vector; }
//   float GetDistance() { return mean_distance; }
//   vector<unsigned int> GetVecPartonsSerialNumber() { return vec_partons_serial_number; }
//   vector<float> GetVecPartonsBaryonNumber() { return vec_partons_baryon_number; }
//   vector<TVector2> GetVecPartonsPosition() { return vec_partons_position; }
// };

// class HardonCandidate : public Hardon {
//  public:
//   friend bool operator<(const HardonCandidate& a, const HardonCandidate& b) { return a.mean_distance < b.mean_distance; }
// };

float GetDistence(float x0, float y0, float x1, float y1) { return rM * sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1)); }
float GetDistence(TVector2 p0, TVector2 p1) {
  float x0 = p0.X();
  float y0 = p0.Y();
  float x1 = p1.X();
  float y1 = p1.Y();
  return sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));
}
float GetDistence(TVector2 p0, TVector2 p1, TVector2 p2) {
  float x0 = p0.X();
  float y0 = p0.Y();
  float x1 = p1.X();
  float y1 = p1.Y();
  float x2 = p1.X();
  float y2 = p1.Y();
  return 1. / 3 *
         (sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1)) + sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)) +
          sqrt((x2 - x0) * (x2 - x0) + (y2 - y0) * (y2 - y0)));
}

const int MaxNumTracks = 100000;
typedef struct {
  unsigned int nevent;
  unsigned int multi;
  int fBaryonNumber[MaxNumTracks];
  float fHadronX[MaxNumTracks];
  float fHadronY[MaxNumTracks];
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

int main() {
  // Generate 10,000 events
  TTree* tree = new TTree("CoalData", "CoalData DST Tree");
  Cell_t cell;
  tree->Branch("Event", &cell.nevent, "nevent/I:multi/I");
  tree->Branch("BaryonNum", cell.fBaryonNumber, "BaryonNum[multi]/I");
  tree->Branch("HadronX", cell.fHadronX, "HadronX[multi]/F");
  tree->Branch("HadronY", cell.fHadronY, "HadronY[multi]/F");
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
  for (size_t iEvent = 0; iEvent < 20000; iEvent++) {
    if (iEvent % 100 == 0) std::cout << "Processing event # " << iEvent << endl;
    // cell.nevent = iEvent;
    // Seed Partons
    vector<Parton> vecPartons;
    unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
    TRandom3* rndm = new TRandom3(seed);
    unsigned int n = 0;
    for (int i = 0; i < nPartons; i++) {
      Parton parton;
      parton.SetBaryonNumber(1. / 3.);
      float r = rndm->Uniform(0, 1);
      float phi = rndm->Uniform(0, TMath::TwoPi());
      parton.SetPositionVector(sqrt(r) * cos(phi), sqrt(r) * sin(phi));
      parton.SetSerialNumber(n);
      vecPartons.push_back(parton);
      n++;
    }
    for (int i = 0; i < nAntiPartons; i++) {
      Parton antiparton;
      antiparton.SetBaryonNumber(-1. / 3.);
      float r = rndm->Uniform(0, 1);
      float phi = rndm->Uniform(0, TMath::TwoPi());
      antiparton.SetPositionVector(sqrt(r) * cos(phi), sqrt(r) * sin(phi));
      antiparton.SetSerialNumber(n);
      vecPartons.push_back(antiparton);
      n++;
    }
    delete rndm;
    rndm = nullptr;

    vector<HardonCandidate> vecHardonCandidates;
    // Get Meson Candidates
    unsigned int serial_number = 0;
    for (size_t i = 0; i < vecPartons.size(); i++) {
      for (size_t j = i + 1; j < vecPartons.size(); j++) {
        bool ismeson = fabs(vecPartons[i].GetBaryonNumber() + vecPartons[j].GetBaryonNumber()) < 1.e-6;
        if (ismeson) {
          float distance = GetDistence(vecPartons[i].GetPositionVector(), vecPartons[j].GetPositionVector());
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

    // sort the hardon candidates by the distance
    sort(vecHardonCandidates.begin(), vecHardonCandidates.end());
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
      } else {
        cell.fParton2BaryonNumber[i] = -999;
        cell.fParton2X[i] = -999;
        cell.fParton2Y[i] = -999;
      }
      i++;
    }
    tree->Fill();
  }  // generate events end;
  TFile* file = new TFile("CoaleEvent.root", "RECREATE");
  file->cd();
  tree->Write();
  file->Close();
}  // main end
