#include <algorithm>
#include <chrono>  // std::chrono::system_clock
#include <cmath>
#include <iostream>
#include <random>  // std::default_random_engine
#include <vector>

#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TVector2.h"

#include "Particle.h"

using namespace std;
const int nPartons = 24;
const int nAntiPartons = 24;
const float rM = 1.;

// TODO
// TVector2 GetCentre(TVector2 p0, TVector2 p1, TVector2 p2) {
//   float x0 = p0.X();
//   float y0 = p0.Y();
//   float x1 = p1.X();
//   float y1 = p1.Y();
//   float x2 = p1.X();
//   float y2 = p1.Y();

//   float d1 = (x1*x1 + y1*y1) - (x0*x0 + y0*y0);
//   float d2 = (x2*x2 + y2*y2) - (x1*x1 + y1*y1);
//   float fm = 2.*((y2-y1)*(x1-x0) - (y1-y0)*(x2-x1));

//   float x = ((y2-y1) * d1 - (y1-y0) * d2)/fm;
//   float y = ((y2-y1) * d1 - (y1-y0) * d2)/fm;

//   TVector2 centre(x,y);
//   return centre;
// }

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

int main() {
  TFile* outputfile = new TFile("ACoaleEvent.root", "RECREATE");
  TGraph* g_partons[2];
  g_partons[0] = new TGraph(nPartons);
  g_partons[1] = new TGraph(nAntiPartons);
  g_partons[1] -> SetName("partons");
  g_partons[0] -> SetName("antipartons");

  vector<Parton> vecPartons;
  // Seed Partons
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
    g_partons[0]->SetPoint(i, parton.GetPositionVector().X(), parton.GetPositionVector().Y());
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
    g_partons[1]->SetPoint(i, antiparton.GetPositionVector().X(), antiparton.GetPositionVector().Y());
  }
  delete rndm;
  rndm = nullptr;

  for (size_t i = 0; i < 2; i++) { 
    g_partons[i]->Write();
    delete g_partons[i];
    g_partons[i] = nullptr;
  }

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

        TGraph* gCandi = new TGraph(2);
        gCandi->SetPoint(0, vecPartons[i].GetPositionVector().X(), vecPartons[i].GetPositionVector().Y());
        gCandi->SetPoint(1, vecPartons[j].GetPositionVector().X(), vecPartons[j].GetPositionVector().Y());
        gCandi->SetName(Form("meson_candi_%u", serial_number));
        gCandi->Write();
        delete gCandi;
        gCandi = nullptr;

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

          TGraph* gCandi = new TGraph(3);
          gCandi->SetPoint(0, vecPartons[i].GetPositionVector().X(), vecPartons[i].GetPositionVector().Y());
          gCandi->SetPoint(1, vecPartons[j].GetPositionVector().X(), vecPartons[j].GetPositionVector().Y());
          gCandi->SetPoint(2, vecPartons[k].GetPositionVector().X(), vecPartons[k].GetPositionVector().Y());
          if (isbayron) gCandi->SetName(Form("baryon_candi_%u", serial_number));
          if (isantibaryon) gCandi->SetName(Form("antibaryon_candi_%u", serial_number));
          gCandi->Write();
          delete gCandi;
          gCandi = nullptr;

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

  unsigned int m = 0, b = 0, ab = 0;
  for (auto hardon : vecHardons) {
    if (hardon.GetBaryonNumber() == 0) {
      TGraph* gMeson = new TGraph(1);
      gMeson->SetName(Form("meson_%i", m));
      gMeson->SetPoint(0,hardon.GetPositionVector().X(),hardon.GetPositionVector().Y());
      gMeson->Write();
      delete gMeson;
      gMeson = nullptr;
      TGraph* gMesonPartons = new TGraph(2);
      gMesonPartons->SetName(Form("meson_partons_%i", m));
      gMesonPartons->SetPoint(0,hardon.GetVecPartonsPosition()[0].X(),hardon.GetVecPartonsPosition()[0].Y());
      gMesonPartons->SetPoint(1,hardon.GetVecPartonsPosition()[1].X(),hardon.GetVecPartonsPosition()[1].Y());
      gMesonPartons->Write();
      delete gMesonPartons;
      gMesonPartons = nullptr;
      m++;
    }
    if (hardon.GetBaryonNumber() == 1) {
      TGraph* gBaryon = new TGraph(1);
      gBaryon->SetName(Form("baryon_%i", b));
      gBaryon->SetPoint(0,hardon.GetPositionVector().X(),hardon.GetPositionVector().Y());
      gBaryon->Write();
      delete gBaryon;
      gBaryon = nullptr;
      TGraph* gBaryonPartons = new TGraph(4);
      gBaryonPartons->SetName(Form("baryon_partons_%i", b));
      gBaryonPartons->SetPoint(0,hardon.GetVecPartonsPosition()[0].X(),hardon.GetVecPartonsPosition()[0].Y());
      gBaryonPartons->SetPoint(1,hardon.GetVecPartonsPosition()[1].X(),hardon.GetVecPartonsPosition()[1].Y());
      gBaryonPartons->SetPoint(2,hardon.GetVecPartonsPosition()[2].X(),hardon.GetVecPartonsPosition()[2].Y());
      gBaryonPartons->SetPoint(3,hardon.GetVecPartonsPosition()[0].X(),hardon.GetVecPartonsPosition()[0].Y());
      gBaryonPartons->Write();
      delete gBaryonPartons;
      gBaryonPartons = nullptr;
      b++;
    }
    if (hardon.GetBaryonNumber() == -1) {
      TGraph* gAntiBaryon = new TGraph(1);
      gAntiBaryon->SetName(Form("antibaryon_%i", ab));
      gAntiBaryon->SetPoint(0,hardon.GetPositionVector().X(),hardon.GetPositionVector().Y());
      gAntiBaryon->Write();
      delete gAntiBaryon;
      gAntiBaryon = nullptr;
      TGraph* gAntiBaryonPartons = new TGraph(4);
      gAntiBaryonPartons->SetName(Form("antibaryon_partons_%i", ab));
      gAntiBaryonPartons->SetPoint(0,hardon.GetVecPartonsPosition()[0].X(),hardon.GetVecPartonsPosition()[0].Y());
      gAntiBaryonPartons->SetPoint(1,hardon.GetVecPartonsPosition()[1].X(),hardon.GetVecPartonsPosition()[1].Y());
      gAntiBaryonPartons->SetPoint(2,hardon.GetVecPartonsPosition()[2].X(),hardon.GetVecPartonsPosition()[2].Y());
      gAntiBaryonPartons->SetPoint(3,hardon.GetVecPartonsPosition()[0].X(),hardon.GetVecPartonsPosition()[0].Y());
      gAntiBaryonPartons->Write();
      delete gAntiBaryonPartons;
      gAntiBaryonPartons = nullptr;
      ab++;
    }
  }
  outputfile->Close();
}  // main end
