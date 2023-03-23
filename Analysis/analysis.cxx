//--------------------------------
// load header
//--------------------------------

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ class PlotFile;
#endif
#ifndef __CINT__
#include <TFile.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <algorithm> // std::shuffle
#include <chrono>    // std::chrono::system_clock
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random> // std::default_random_engine
#include <vector>

#include "TBits.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TComplex.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TObject.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TUnixSystem.h"
#include "TVector.h"
#include "TVector2.h"
#include "TVector3.h"
#include "math.h"
#include "string.h"
#endif
#include "CoalData.h"
using namespace std;

//--------------------------------
// global variables
//--------------------------------
double RangePhi(double angle)
{
  while (angle >= 1.5 * TMath::Pi())
    angle -= TMath::TwoPi();
  while (angle < -0.5 * TMath::Pi())
    angle += TMath::TwoPi();
  return angle;
}
//--------------------------------
// main function starts here
//--------------------------------
int main(int argc, char **argv)
{
  if (argc != 3)
    return 1;

  TString inputFile;
  TString outputFile;
  if (argc == 3)
  {
    inputFile = argv[1];
    outputFile = argv[2];
  }
  //----------------------------------
  // open files and add to chain
  //----------------------------------
  int fileNumber = 0;
  char fileList[512];
  TChain *chain = new TChain("CoalData");
  if (inputFile.Contains(".list"))
  {
    ifstream *inputStream = new ifstream;
    inputStream->open(inputFile);
    if (!(inputStream))
    {
      std::cout << "can not open file list" << endl;
      return 1;
    }
    while (inputStream->good())
    {
      inputStream->getline(fileList, 512);
      if (inputStream->eof())
        break;
      TFile *fTmp = new TFile(fileList);
      if (!fTmp || !(fTmp->IsOpen()) || !(fTmp->GetNkeys()))
      {
        std::cout << "open file list error" << endl;
        return 1;
      }
      else
      {
        std::cout << "reading file " << fileList << endl;
        chain->Add(fileList);
        fileNumber++;
      }
      delete fTmp;
    }
    std::cout << fileNumber << " files read in" << endl;
  }
  else if (inputFile.Contains(".root"))
  {
    chain->Add(inputFile.Data());
  }
  //--------------------------------
  // define global histograms
  //--------------------------------
  TH1D *h_bn = new TH1D("h_bn", "h_bn", 3, -1, 2);
  TH1D *h_parton_phi[2];
  for (int i = 0; i < 2; i++)
  {
    h_parton_phi[i] = new TH1D(Form("h_parton_phi_%i", i), Form("h_parton_phi_%i", i), 100, 0, TMath::TwoPi());
  }
  TH1F *h_x[3];
  TH1F *h_y[3];
  TH1F *h_r[3];
  TH1F *h_phix[3];
  TH1F *h_px[3];
  TH1F *h_py[3];
  TH1F *h_pT[3];
  TH1F *h_phip[3];
  TH1F *h_distance[3];
  for (int i = 0; i < 3; i++)
  {
    h_x[i] = new TH1F(Form("h_x_%i", i), Form("h_x_%i", i), 100, -1, 1);
    h_y[i] = new TH1F(Form("h_y_%i", i), Form("h_y_%i", i), 100, -1, 1);
    h_r[i] = new TH1F(Form("h_r_%i", i), Form("h_r_%i", i), 100, 0., 1.0);
    h_phix[i] = new TH1F(Form("h_phix_%i", i), Form("h_phix_%i", i), 100, 0., TMath::TwoPi());
    h_px[i] = new TH1F(Form("h_px_%i", i), Form("h_px_%i", i), 100, -3, 3);
    h_py[i] = new TH1F(Form("h_py_%i", i), Form("h_py_%i", i), 100, -3, 3);
    h_pT[i] = new TH1F(Form("h_pT_%i", i), Form("h_pT_%i", i), 100, 0., 8.);
    h_phip[i] = new TH1F(Form("h_phip_%i", i), Form("h_phip_%i", i), 100, 0., TMath::TwoPi());
    h_distance[i] = new TH1F(Form("h_distance_%i", i), Form("h_distance_%i", i), 100, 0, 3);
  }
  TH1F *h_del_phi_coaled_partons_meson = new TH1F("h_del_phi_coaled_partons_meson", "h_del_phi_coaled_partons_meson", 50, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH1F *h_del_phi_coaled_partons_baryon = new TH1F("h_del_phi_coaled_partons_baryon", "h_del_phi_coaled_partons_baryon", 50, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH1F *h_del_phi_coaled_partons_antibaryon = new TH1F("h_del_phi_coaled_partons_antibaryon", "h_del_phi_coaled_partons_antibaryon", 50, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH1F *h_CDelPhiX_meson_particle[3]; //[m-m][m-b][m-bbar]
  TH1F *h_CDelPhiP_meson_particle[3];
  for (int i = 0; i < 3; i++)
  {
    h_CDelPhiX_meson_particle[i] = new TH1F(Form("h_CDelPhiX_meson_particle_%i", i), Form("h_CDelPhiX_meson_particle_%i", i), 50, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
    h_CDelPhiP_meson_particle[i] = new TH1F(Form("h_CDelPhiP_meson_particle_%i", i), Form("h_CDelPhiP_meson_particle_%i", i), 50, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  }
  TH1F *h_CDelPhiX_baryon_baryon[3]; //[b-b][b-bbar][bbar-bbar]
  TH1F *h_CDelPhiP_baryon_baryon[3];
  for (int i = 0; i < 3; i++)
  {
    h_CDelPhiX_baryon_baryon[i] = new TH1F(Form("h_CDelPhiX_baryon_baryon_%i", i), Form("h_CDelPhiX_baryon_baryon_%i", i), 50, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
    h_CDelPhiP_baryon_baryon[i] = new TH1F(Form("h_CDelPhiP_baryon_baryon_%i", i), Form("h_CDelPhiP_baryon_baryon_%i", i), 50, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  }

  TH1::SetDefaultSumw2(kTRUE);

  //--------------------------------
  // loop events
  //--------------------------------
  int nEvents = (int)chain->GetEntries();
  //nEvents = 300000;
  std::cout << "Total events : " << nEvents << endl;
  CoalData *coal = new CoalData(chain);
  for (int iEvent = 0; iEvent < nEvents; ++iEvent)
  {
    vector<float> vecPhiProton;
    vector<float> vecPhiAntiProton;
    coal->GetEntry(iEvent);
    if (iEvent % 10000 == 0)
      std::cout << "Processing event # " << iEvent << endl;
    int nTracks = coal->Event_multi;
    for (int iTrack = 0; iTrack < nTracks; ++iTrack)
    {
      float x_0 = coal->HadronX[iTrack];
      float y_0 = coal->HadronY[iTrack];
      TVector2 X_0(x_0, y_0);
      float phix_0 = X_0.Phi();
      float r_0 = X_0.Mod();

      float px_0 = coal->HadronPx[iTrack];
      float py_0 = coal->HadronPy[iTrack];
      TVector2 P_0(px_0, py_0);
      float phip_0 = P_0.Phi();
      float pT_0 = P_0.Mod();

      float distance = coal->Distance[iTrack];
      // partons 0
      float bn_0_0 = coal->Parton0BaryonNum[iTrack];
      float x_0_0 = coal->Parton0X[iTrack];
      float y_0_0 = coal->Parton0Y[iTrack];
      TVector2 X_0_0(x_0_0, y_0_0);
      float phi_0_0 = X_0_0.Phi();
      if (bn_0_0 > 0)
        h_parton_phi[0]->Fill(phi_0_0);
      if (bn_0_0 < 0)
        h_parton_phi[1]->Fill(phi_0_0);
      // partons 1
      float bn_0_1 = coal->Parton1BaryonNum[iTrack];
      float x_0_1 = coal->Parton1X[iTrack];
      float y_0_1 = coal->Parton1Y[iTrack];
      TVector2 X_0_1(x_0_1, y_0_1);
      float phi_0_1 = X_0_1.Phi();
      if (bn_0_1 > 0)
        h_parton_phi[0]->Fill(phi_0_1);
      if (bn_0_1 < 0)
        h_parton_phi[1]->Fill(phi_0_1);
      // partons 2
      float bn_0_2 = coal->Parton2BaryonNum[iTrack];
      if (bn_0_2 > -1)
      {
        float x_0_2 = coal->Parton2X[iTrack];
        float y_0_2 = coal->Parton2Y[iTrack];
        TVector2 X_0_2(x_0_2, y_0_2);
        float phi_0_2 = X_0_2.Phi();
        if (bn_0_0 > 0)
          h_parton_phi[0]->Fill(phi_0_2);
        if (bn_0_0 < 0)
          h_parton_phi[1]->Fill(phi_0_2);
      }

      int baryon_number_0 = coal->BaryonNum[iTrack];
      h_bn->Fill(baryon_number_0);
      if (baryon_number_0 == 0)
      {
        h_x[0]->Fill(x_0);
        h_y[0]->Fill(y_0);
        h_r[0]->Fill(r_0);
        h_phix[0]->Fill(phix_0);
        h_px[0]->Fill(px_0);
        h_py[0]->Fill(py_0);
        h_pT[0]->Fill(pT_0);
        h_phip[0]->Fill(phip_0);
        h_distance[0]->Fill(distance);
        h_del_phi_coaled_partons_meson->Fill(RangePhi(phi_0_0 - phi_0_1));
      }
      if (baryon_number_0 == 1)
      {
        h_x[1]->Fill(x_0);
        h_y[1]->Fill(y_0);
        h_r[1]->Fill(r_0);
        h_phix[1]->Fill(phix_0);
        h_px[1]->Fill(px_0);
        h_py[1]->Fill(py_0);
        h_pT[1]->Fill(pT_0);
        h_phip[1]->Fill(phip_0);
        h_distance[1]->Fill(distance);
        h_del_phi_coaled_partons_baryon->Fill(RangePhi(phi_0_0 - phi_0_1));
      }
      if (baryon_number_0 == -1)
      {
        h_x[2]->Fill(x_0);
        h_y[2]->Fill(y_0);
        h_r[2]->Fill(r_0);
        h_phix[2]->Fill(phix_0);
        h_px[2]->Fill(px_0);
        h_py[2]->Fill(py_0);
        h_pT[2]->Fill(pT_0);
        h_phip[2]->Fill(phip_0);
        h_distance[2]->Fill(distance);
        h_del_phi_coaled_partons_antibaryon->Fill(RangePhi(phi_0_0 - phi_0_1));
      }

      for (int jTrack = 0; jTrack < nTracks; jTrack++)
      {
        if (iTrack == jTrack)
          continue;
        float x_1 = coal->HadronX[jTrack];
        float y_1 = coal->HadronY[jTrack];
        TVector2 X_1(x_1, y_1);
        float phix_1 = X_1.Phi();

        float px_1 = coal->HadronPx[jTrack];
        float py_1 = coal->HadronPy[jTrack];
        TVector2 P_1(px_1, py_1);
        float phip_1 = P_1.Phi();

        double del_phix = RangePhi(phix_0 - phix_1);
        double del_phip = RangePhi(phip_0 - phip_1);

        int baryon_number_1 = coal->BaryonNum[jTrack];
        if ((baryon_number_0 == 0) && (baryon_number_1 == 0))
        {
          h_CDelPhiX_meson_particle[0]->Fill(del_phix);
          h_CDelPhiP_meson_particle[0]->Fill(del_phip);
        }
        if ((baryon_number_0 + baryon_number_1) == 1)
        {
          h_CDelPhiX_meson_particle[1]->Fill(del_phix);
          h_CDelPhiP_meson_particle[1]->Fill(del_phip);
        }
        if ((baryon_number_0 + baryon_number_1) == -1)
        {
          h_CDelPhiX_meson_particle[2]->Fill(del_phix);
          h_CDelPhiP_meson_particle[2]->Fill(del_phip);
        }
        if ((baryon_number_0 == 1) && (baryon_number_1 == 1))
        {
          h_CDelPhiX_baryon_baryon[0]->Fill(del_phix);
          h_CDelPhiP_baryon_baryon[0]->Fill(del_phip);
        }
        if ((baryon_number_0 * baryon_number_1) == -1)
        {
          h_CDelPhiX_baryon_baryon[1]->Fill(del_phix);
          h_CDelPhiP_baryon_baryon[1]->Fill(del_phip);
        }
        if ((baryon_number_0 == -1) && (baryon_number_1 == -1))
        {
          h_CDelPhiX_baryon_baryon[2]->Fill(del_phix);
          h_CDelPhiP_baryon_baryon[2]->Fill(del_phip);
        }
      }
    }
  } // 本事件结束

  TFile *f = new TFile(outputFile, "RECREATE");
  f->cd();
  h_bn->Write();
  for (int i = 0; i < 2; i++)
  {
    h_parton_phi[i]->Write();
  }

  for (int i = 0; i < 3; i++)
  {
    h_x[i]->Write();
    h_y[i]->Write();
    h_r[i]->Write();
    h_phix[i]->Write();
    h_px[i]->Write();
    h_py[i]->Write();
    h_pT[i]->Write();
    h_phip[i]->Write();
    h_distance[i]->Write();
  }
  h_del_phi_coaled_partons_meson->Write();
  h_del_phi_coaled_partons_baryon->Write();
  h_del_phi_coaled_partons_antibaryon->Write();
  for (int i = 0; i < 3; i++)
  {
    h_CDelPhiX_meson_particle[i]->Scale(h_CDelPhiX_meson_particle[i]->GetNbinsX() / h_CDelPhiX_meson_particle[i]->GetEntries());
    h_CDelPhiP_meson_particle[i]->Scale(h_CDelPhiP_meson_particle[i]->GetNbinsX() / h_CDelPhiP_meson_particle[i]->GetEntries());
    h_CDelPhiX_meson_particle[i]->Write();
    h_CDelPhiP_meson_particle[i]->Write();
  }
  for (int i = 0; i < 3; i++)
  {
    h_CDelPhiX_baryon_baryon[i]->Scale(h_CDelPhiX_baryon_baryon[i]->GetNbinsX() / h_CDelPhiX_baryon_baryon[i]->GetEntries());
    h_CDelPhiP_baryon_baryon[i]->Scale(h_CDelPhiP_baryon_baryon[i]->GetNbinsX() / h_CDelPhiP_baryon_baryon[i]->GetEntries());
    h_CDelPhiX_baryon_baryon[i]->Write();
    h_CDelPhiP_baryon_baryon[i]->Write();
  }
  f->Close();
  std::cout << ".root has been created" << endl;
}
