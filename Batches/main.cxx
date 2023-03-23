#include <algorithm>
#include <chrono>  // std::chrono::system_clock
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>  // std::default_random_engine
#include <vector>

#include "Particle.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPolyLine.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include "TVector2.h"

using namespace std;
const float mass_pion = 0.13957039;
const float mass_proton = 0.93827208816;
const bool shuffle_mode = false;
const bool global_sort_mode = false;  // very slow for a large nPartons
const bool debug = true;

void GetRestPartons(vector<Parton>& vecPartonsRest, vector<Parton>& vecPartons, vector<unsigned int>& vecNodedBuffer);
void GetHadrons(vector<Hardon>& vecHardons, vector<Parton>& vecPartons, vector<unsigned int>& vecNodedBuffer, float rM);

float GetDistence(TVector2 p0, TVector2 p1) {
  float dis = (p0 - p1).Mod();
  return dis;
}

float GetDistence(TVector2 p0, TVector2 p1, TVector2 p2) {
  float dis = 1. / 3. * ((p0 - p1).Mod() + (p1 - p2).Mod() + (p2 - p0).Mod());
  return dis;
}

int GetMeshCoordinate(double x) {
  if (fabs(x) > 1.) return nan("");
  if (x < -1. / 3.) return 0;
  if (x > +1. / 3.) return 2;
  return 1;
}

Parton GetLBCFriendForThisParton(Parton parton, float sigmaLBC, float bn_friend, unsigned int sm_offset) {
  Parton partonLBC;
  partonLBC.SetBaryonNumber(bn_friend);
  float r = gRandom->Uniform(0, 1);
  float phi = parton.GetPositionVector().Phi();
  float dphi = gRandom->Gaus(0, sigmaLBC);
  phi += dphi;
  partonLBC.SetPositionVector(sqrt(r) * cos(phi), sqrt(r) * sin(phi));
  partonLBC.SetSerialNumber(parton.GetSerialNumber() + sm_offset);
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
  if (!filein) return 0;

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
  cout << "frac_qq      = " << frac_qq << endl;
  cout << "frac_qqbar   = " << frac_qqbar << endl;
  cout << "sigmaLBC     = " << sigmaLBC << endl;
  sigmaLBC *= TMath::Pi();

  if (debug) {
    nEvents = 1;
    std::cout << "Debug Mode nEvent is fixed to " << nEvents << std::endl;
  }

  // For Debug
  TGraph* gParton = new TGraph();
  TGraph* gAntiParton = new TGraph();
  TGraph* gMeson = new TGraph();
  TGraph* gBaryon = new TGraph();
  TGraph* gAntiBaryon = new TGraph();
  gParton->SetName("gParton");
  gAntiParton->SetName("gAntiParton");
  gMeson->SetName("gMeson");
  gBaryon->SetName("gBaryon");
  gAntiBaryon->SetName("gAntiBaryon");
  vector<TGraph*> vecGraphPartonsInMeson;
  vector<TGraph*> vecGraphPartonsInBaryon;
  vector<TGraph*> vecGraphPartonsInAntiBaryon;

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

    // to record the Partons
    vector<Parton> vecPartons;
    // to record the Hardons
    vector<Hardon> vecHardons;

    // StopWatch for test efficiency
    TStopwatch timer;
    if (debug) timer.Start();

    // Seed Partons
    vector<Parton> vecPartonsTemp;
    vector<Parton> vecAntiPartonsTemp;
    unsigned int n = 1;  // serial number
    for (size_t i = 0; i < nPartons; i++) {
      Parton parton;
      parton.SetBaryonNumber(1. / 3.);
      float r = rndm->Uniform(0, 1);
      float phi = rndm->Uniform(0, TMath::TwoPi());
      parton.SetPositionVector(sqrt(r) * cos(phi), sqrt(r) * sin(phi));
      parton.SetSerialNumber(n);
      vecPartonsTemp.push_back(parton);
      if (rndm->Uniform(0, 1) < frac_qq) vecPartonsTemp.push_back(GetLBCFriendForThisParton(parton, sigmaLBC, 1. / 3., 1000));
      if (rndm->Uniform(0, 1) < frac_qqbar) vecAntiPartonsTemp.push_back(GetLBCFriendForThisParton(parton, sigmaLBC, -1. / 3., 2000));
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
      if (rndm->Uniform(0, 1) < frac_qq) vecAntiPartonsTemp.push_back(GetLBCFriendForThisParton(antiparton, sigmaLBC, -1. / 3., 1000));
      if (rndm->Uniform(0, 1) < frac_qqbar) vecPartonsTemp.push_back(GetLBCFriendForThisParton(antiparton, sigmaLBC, 1. / 3., 2000));
      n++;
    }

    shuffle(vecPartonsTemp.begin(), vecPartonsTemp.end(), std::default_random_engine(seed + 1));
    shuffle(vecAntiPartonsTemp.begin(), vecAntiPartonsTemp.end(), std::default_random_engine(seed + 2));
    vecPartonsTemp.resize(nPartons);
    vecAntiPartonsTemp.resize(nAntiPartons);
    vecPartons.insert(vecPartons.end(), vecPartonsTemp.begin(), vecPartonsTemp.end());
    vecPartons.insert(vecPartons.end(), vecAntiPartonsTemp.begin(), vecAntiPartonsTemp.end());

    if (debug) {
      timer.Stop();
      std::cout << "Time for seeding partons: " << timer.RealTime() << endl;
      timer.Start();
    }

    if (debug) {
      int iPartons = 0, iAntiPartons = 0;
      for (auto parton : vecPartons) {
        if (parton.GetBaryonNumber() > 0)
          gParton->SetPoint(iPartons++, parton.GetPositionVector().X(), parton.GetPositionVector().Y());
        else
          gAntiParton->SetPoint(iAntiPartons++, parton.GetPositionVector().X(), parton.GetPositionVector().Y());
      }
    }

    // Mesh Partons in 3x3 Grids
    vector<Parton> vecPartonsMeshed[3][3];
    if (!global_sort_mode) {
      for (auto parton : vecPartons) {
        double x = parton.GetPositionVector().X();
        double y = parton.GetPositionVector().Y();
        int iGrid = GetMeshCoordinate(x);
        int jGrid = GetMeshCoordinate(y);
        if (isnan(iGrid) || isnan(jGrid)) {
          std::cout << "Error: Parton's position out of range!" << endl;
          return 1;
        }
        vecPartonsMeshed[iGrid][jGrid].push_back(parton);
      }
    }

    if (debug) {
      timer.Stop();
      std::cout << "Time for meshing partons: " << timer.RealTime() << endl;
      timer.Start();
    }

    // Get Hadrons
    // vector to record the serial numbers of the used partons
    vector<unsigned int> vecUsedPartons;
    // buffer to record the serial numbers of the used partons in this loop
    vector<unsigned int> vecNodedBuffer;

    if (global_sort_mode) {
      GetHadrons(vecHardons, vecPartons, vecUsedPartons, rM);
      if (debug) cout << "Number of partons not been formed to hadrons: " << vecPartons.size() - vecUsedPartons.size() << endl;
      vector<unsigned int>().swap(vecUsedPartons);
    } else {
      vector<Parton> vecPartonsRest;
      for (size_t iGrid = 0; iGrid < 3; iGrid++) {
        for (size_t jGrid = 0; jGrid < 3; jGrid++) {
          vector<unsigned int>().swap(vecNodedBuffer);
          GetHadrons(vecHardons, vecPartonsMeshed[iGrid][jGrid], vecNodedBuffer, rM);
          GetRestPartons(vecPartonsRest, vecPartonsMeshed[iGrid][jGrid], vecNodedBuffer);
          vecUsedPartons.insert(vecUsedPartons.end(), vecNodedBuffer.begin(), vecNodedBuffer.end());
        }
      }
      vector<unsigned int>().swap(vecNodedBuffer);
      if (debug) cout << "Number of partons not been used after coalescence in grid: " << vecPartons.size() - vecUsedPartons.size() << endl;
      GetHadrons(vecHardons, vecPartonsRest, vecNodedBuffer, rM);
      vecUsedPartons.insert(vecUsedPartons.end(), vecNodedBuffer.begin(), vecNodedBuffer.end());
      if (debug) cout << "Number of partons not been used finaly: " << vecPartons.size() - vecUsedPartons.size() << endl;
      vector<unsigned int>().swap(vecNodedBuffer);
    }

    if (debug) {
      timer.Stop();
      std::cout << "Time for getting hadrons: " << timer.RealTime() << endl;
      timer.Start();
    }

    // Fill the Tree and Sample Momentum
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
        hardon.MomentumBoost(rho0, rho2);
        cell.fHadronPx[i] = hardon.GetP().Px();
        cell.fHadronPy[i] = hardon.GetP().Py();
      } else {
        cell.fParton2BaryonNumber[i] = -999;
        cell.fParton2X[i] = -999;
        cell.fParton2Y[i] = -999;
        hardon.SampleRawMomentum(mass_pion);
        hardon.MomentumBoost(rho0, rho2);
        cell.fHadronPx[i] = hardon.GetP().Px();
        cell.fHadronPy[i] = hardon.GetP().Py();
      }
      i++;
    }
    tree->Fill();

    if (debug) {
      timer.Stop();
      std::cout << "Time for filling tree: " << timer.RealTime() << endl;
    }

    if (debug) {
      int iMeson = 0, iBaryon = 0, iAntiBaryon = 0;
      for (auto hadron : vecHardons) {
        if (hadron.GetBaryonNumber() == 0) {
          gMeson->SetPoint(iMeson++, hadron.GetPositionVector().X(), hadron.GetPositionVector().Y());
          TGraph* partonsInMeson = new TGraph();
          for (size_t i = 0; i < 2; i++) partonsInMeson->SetPoint(i, hadron.GetVecPartonsPosition()[i].X(), hadron.GetVecPartonsPosition()[i].Y());
          vecGraphPartonsInMeson.push_back(partonsInMeson);

        } else if (hadron.GetBaryonNumber() == 1) {
          gBaryon->SetPoint(iBaryon++, hadron.GetPositionVector().X(), hadron.GetPositionVector().Y());
          TGraph* partonsInBaryon = new TGraph();
          for (size_t i = 0; i < 3; i++) partonsInBaryon->SetPoint(i, hadron.GetVecPartonsPosition()[i].X(), hadron.GetVecPartonsPosition()[i].Y());
          partonsInBaryon->SetPoint(3, hadron.GetVecPartonsPosition()[0].X(), hadron.GetVecPartonsPosition()[0].Y());
          vecGraphPartonsInBaryon.push_back(partonsInBaryon);

        } else if (hadron.GetBaryonNumber() == -1) {
          gAntiBaryon->SetPoint(iAntiBaryon++, hadron.GetPositionVector().X(), hadron.GetPositionVector().Y());
          TGraph* partonsInAntiBaryon = new TGraph();
          for (size_t i = 0; i < 3; i++) partonsInAntiBaryon->SetPoint(i, hadron.GetVecPartonsPosition()[i].X(), hadron.GetVecPartonsPosition()[i].Y());
          partonsInAntiBaryon->SetPoint(3, hadron.GetVecPartonsPosition()[0].X(), hadron.GetVecPartonsPosition()[0].Y());
          vecGraphPartonsInAntiBaryon.push_back(partonsInAntiBaryon);

        } else {
          std::cout << "Error: Baryon Number is not 0, 1 or -1!" << endl;
          return 1;
        }
      }
    }
  }  // generate events end;

  delete rndm;
  rndm = nullptr;

  // Write the Tree to File
  TFile* file = new TFile(FileOutput, "RECREATE");
  file->cd();
  tree->Write();
  file->Close();

  if (debug) {
    // Draw the Particles' Distribution
    gStyle->SetOptStat(0);
    // Load Color Palette
    int ci[6];
    TColor* color[6];
    ci[0] = TColor::GetFreeColorIndex();
    color[0] = new TColor(ci[0], 240 / 255., 102 / 255., 70 / 255.);  // red
    ci[1] = TColor::GetFreeColorIndex();
    color[1] = new TColor(ci[1], 79 / 255., 194 / 255., 216 / 255.);  // blue
    ci[2] = TColor::GetFreeColorIndex();
    color[2] = new TColor(ci[2], 254 / 255., 198 / 255., 101 / 255.);  // yellow
    ci[3] = TColor::GetFreeColorIndex();
    color[3] = new TColor(ci[3], 146 / 255., 100 / 255., 140 / 255.);  // purple
    ci[4] = TColor::GetFreeColorIndex();
    color[4] = new TColor(ci[4], 125 / 255., 200 / 255., 165 / 255.);  // green
    ci[5] = TColor::GetFreeColorIndex();
    color[5] = new TColor(ci[5], 64 / 255., 64 / 255., 64 / 255.);  // black

    TLegend* leg = new TLegend(0.77, 0.77, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(gParton, "Parton", "p");
    leg->AddEntry(gAntiParton, "AntiParton", "p");
    leg->AddEntry(gMeson, "Meson", "LP");
    leg->AddEntry(gBaryon, "Baryon", "LP");
    leg->AddEntry(gAntiBaryon, "AntiBaryon", "LP");

    TCanvas* c1 = new TCanvas("c1", "c1", 1000, 1000);
    c1->cd();
    TH2D* dummy = new TH2D("dummy", "", 1, -1.5, 1.5, 1, -1.5, 1.5);
    dummy->Draw("SAME");
    TBox* square = new TBox(-1., -1., 1., 1.);
    square->SetFillColor(0);
    square->SetLineWidth(1);
    square->SetLineColor(kBlack);
    square->Draw("L SAME");

    if (!global_sort_mode) {
      TLine* line[4];
      line[0] = new TLine(-1. / 3., -1., -1. / 3., 1.);
      line[1] = new TLine(1. / 3., -1., 1. / 3., 1.);
      line[2] = new TLine(-1., -1. / 3., 1., -1. / 3.);
      line[3] = new TLine(-1., 1. / 3., 1., 1. / 3.);
      for (int i = 0; i < 4; i++) {
        line[i]->SetLineColor(kBlack);
        line[i]->SetLineWidth(1);
        line[i]->SetLineStyle(2);
        line[i]->Draw("L SAME");
      }
    }

    gParton->SetMarkerColor(kBlack);
    gParton->SetMarkerStyle(kFullCircle);
    gParton->SetMarkerSize(0.5);
    gParton->Draw("P SAME");

    gAntiParton->SetMarkerColor(kBlack);
    gAntiParton->SetMarkerStyle(kOpenCircle);
    gAntiParton->SetMarkerSize(0.5);
    gAntiParton->Draw("P SAME");

    gMeson->SetMarkerColor(ci[2]);
    gMeson->SetMarkerStyle(kFullCircle);
    gMeson->SetMarkerSize(1);
    gMeson->Draw("P SAME");

    gBaryon->SetMarkerColor(ci[0]);
    gBaryon->SetMarkerStyle(kFullCircle);
    gBaryon->SetMarkerSize(1);
    gBaryon->Draw("P SAME");

    gAntiBaryon->SetMarkerColor(ci[1]);
    gAntiBaryon->SetMarkerStyle(kFullCircle);
    gAntiBaryon->SetMarkerSize(1);
    gAntiBaryon->Draw("P SAME");

    for (auto graph : vecGraphPartonsInMeson) {
      graph->SetLineColor(ci[2]);
      graph->SetLineWidth(1);
      graph->Draw("L SAME");
    }
    for (auto graph : vecGraphPartonsInBaryon) {
      graph->SetLineColor(ci[0]);
      graph->SetLineWidth(1);
      graph->Draw("L SAME");
    }
    for (auto graph : vecGraphPartonsInAntiBaryon) {
      graph->SetLineColor(ci[1]);
      graph->SetLineWidth(1);
      graph->Draw("L SAME");
    }
    leg->Draw("SAME");
    c1->SaveAs("ParticleDistribution.pdf");
  }
}  // main end

void GetHadrons(vector<Hardon>& vecHardons, vector<Parton>& vecPartons, vector<unsigned int>& vecNodedBuffer, float rM) {
  // vector to record all the hardon candidates
  vector<HardonCandidate> vecHardonCandidates;

  // Get Meson Candidates
  unsigned int serial_number = 1;
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
  if (shuffle_mode) {
    // shuffle the hardon candidates
    unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
    shuffle(vecHardonCandidates.begin(), vecHardonCandidates.end(), std::default_random_engine(seed));
  } else {
    // sort the hardon candidates by the distance
    sort(vecHardonCandidates.begin(), vecHardonCandidates.end());
  }

  // for every hardon candidate, if it's partons have not been used, it should be freezed out
  for (auto hardon_candidates : vecHardonCandidates) {
    // assume no partons for this hardon candidate has been used
    bool is_every_parton_free = true;
    // loop it's partons. If a parton is found in the vecNodedBuffer, no hardon will form
    for (auto sm : hardon_candidates.GetVecPartonsSerialNumber()) {
      auto it = find(vecNodedBuffer.begin(), vecNodedBuffer.end(), sm);
      if (it != vecNodedBuffer.end()) {  // find a parton in the buffer
        is_every_parton_free = false;
        break;
      }
    }
    // If no parton is found in the vecNodedBuffer, a hardon will be formed
    if (is_every_parton_free) {
      for (auto sm : hardon_candidates.GetVecPartonsSerialNumber()) {
        auto it = find(vecNodedBuffer.begin(), vecNodedBuffer.end(), sm);
        // record the partons have been used in this coalesence
        if (it == vecNodedBuffer.end()) vecNodedBuffer.push_back(sm);
      }
      // package hardons
      Hardon hardon;
      hardon = hardon_candidates;
      vecHardons.push_back(hardon);
    }
  }
}

void GetRestPartons(vector<Parton>& vecPartonsRest, vector<Parton>& vecPartons, vector<unsigned int>& vecNodedBuffer) {
  for (auto parton : vecPartons) {
    if (none_of(vecNodedBuffer.begin(), vecNodedBuffer.end(), [&](unsigned int sn) { return sn == parton.GetSerialNumber(); })) {
      vecPartonsRest.push_back(parton);
    }
  }
}
