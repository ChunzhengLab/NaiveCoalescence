#include "TCanvas.h"
#include "TEllipse.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TStyle.h"

using namespace std;

template <class TH>
void SetStyle(TH& hist, unsigned int color, unsigned int markerStyle, double markerSize = 1, double lineWidth = 1) {
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(markerStyle);
  hist->SetMarkerSize(markerSize);
  hist->SetLineWidth(lineWidth);
}

void Draw() {
  gStyle->SetOptStat(0);
  int ci[6];
  TColor* color[6];
  ci[0] = TColor::GetFreeColorIndex();
  color[0] = new TColor(ci[0], 240 / 255., 102 / 255., 70 / 255.);  //红
  ci[1] = TColor::GetFreeColorIndex();
  color[1] = new TColor(ci[1], 79 / 255., 194 / 255., 216 / 255.);  //蓝
  ci[2] = TColor::GetFreeColorIndex();
  color[2] = new TColor(ci[2], 254 / 255., 198 / 255., 101 / 255.);  //黄
  ci[3] = TColor::GetFreeColorIndex();
  color[3] = new TColor(ci[3], 146 / 255., 100 / 255., 140 / 255.);  //紫
  ci[4] = TColor::GetFreeColorIndex();
  color[4] = new TColor(ci[4], 125 / 255., 200 / 255., 165 / 255.);  //绿
  ci[5] = TColor::GetFreeColorIndex();
  color[5] = new TColor(ci[5], 64 / 255., 64 / 255., 64 / 255.);  //黑

  TFile* f[6];
  f[0] = TFile::Open("o.root");
  f[1] = TFile::Open("output20.root");
  f[2] = TFile::Open("output30.root");
  f[3] = TFile::Open("output60.root");
  f[4] = TFile::Open("output100.root");
  f[5] = TFile::Open("output300.root");

  
  TH1D* h_bn = new TH1D("h_bn", "h_bn", 3, -1, 2);
  TH1D* h_parton_phi[2];
  for (int i = 0; i < 2; i++) {
    h_parton_phi[i] = new TH1D(Form("h_parton_phi_%i",i),Form("h_parton_phi_%i",i),100,0,TMath::TwoPi());
  }
  TH1F* h_x[3];
  TH1F* h_y[3];
  TH1F* h_r[3];
  TH1F* h_phix[3];
  TH1F* h_px[3];
  TH1F* h_py[3];
  TH1F* h_pT[3];
  TH1F* h_phip[3];
  TH1F* h_distance[3];
  for (int i = 0; i < 3; i++) {
   h_x[i] = new TH1F(Form("h_x_%i",i), Form("h_x_%i",i), 100, -1, 1);
   h_y[i] = new TH1F(Form("h_y_%i",i), Form("h_y_%i",i), 100, -1, 1);
   h_r[i] = new TH1F(Form("h_r_%i",i), Form("h_r_%i",i), 100, 0., 1.0);
   h_phix[i] = new TH1F(Form("h_phix_%i",i), Form("h_phix_%i",i), 100, 0., TMath::TwoPi());
   h_px[i] = new TH1F(Form("h_px_%i",i), Form("h_px_%i",i), 100, -3, 3);
   h_py[i] = new TH1F(Form("h_py_%i",i), Form("h_py_%i",i), 100, -3, 3);
   h_pT[i] = new TH1F(Form("h_pT_%i",i), Form("h_pT_%i",i), 100, 0., 8.);
   h_phip[i] = new TH1F(Form("h_phip_%i",i), Form("h_phip_%i",i), 100, 0., TMath::TwoPi());
   h_distance[i] = new TH1F(Form("h_distance_%i",i), Form("h_distance_%i",i), 100, 0, 3);
  }
  TH1F* h_del_phi_coaled_partons_meson  = new TH1F("h_del_phi_coaled_partons_meson", "h_del_phi_coaled_partons_meson", 50, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH1F* h_del_phi_coaled_partons_baryon = new TH1F("h_del_phi_coaled_partons_baryon", "h_del_phi_coaled_partons_baryon", 50, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH1F* h_del_phi_coaled_partons_antibaryon = new TH1F("h_del_phi_coaled_partons_antibaryon", "h_del_phi_coaled_partons_antibaryon", 50, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  TH1F* h_CDelPhiX_meson_particle[3]; //[m-m][m-b][m-bbar]
  TH1F* h_CDelPhiP_meson_particle[3];
  for (int i = 0; i < 3; i++) {
    h_CDelPhiX_meson_particle[i] = new TH1F(Form("h_CDelPhiX_meson_particle_%i",i), Form("h_CDelPhiX_meson_particle_%i",i), 50, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
    h_CDelPhiP_meson_particle[i] = new TH1F(Form("h_CDelPhiP_meson_particle_%i",i), Form("h_CDelPhiP_meson_particle_%i",i), 50, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  }
  TH1F* h_CDelPhiX_baryon_baryon[3]; //[b-b][b-bbar][bbar-bbar]
  TH1F* h_CDelPhiP_baryon_baryon[3];
  for (int i = 0; i < 3; i++) {
    h_CDelPhiX_baryon_baryon[i] = new TH1F(Form("h_CDelPhiX_baryon_baryon_%i",i), Form("h_CDelPhiX_baryon_baryon_%i",i), 50, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
    h_CDelPhiP_baryon_baryon[i] = new TH1F(Form("h_CDelPhiP_baryon_baryon_%i",i), Form("h_CDelPhiP_baryon_baryon_%i",i), 50, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
  }

  TH2D* dummy = new TH2D("", "", 1, -0.5 * TMath::Pi(), 1.5 * TMath::Pi(), 1, 0.5, 1.5);
  dummy->SetTitle("#phi for x");
  TH2D* dummy2 = new TH2D("", "", 1, -0.5 * TMath::Pi(), 1.5 * TMath::Pi(), 1, 0.95, 1.05);
  dummy2->SetTitle("#phi for p");
  TH2D* dummy3 = new TH2D("", "", 1, 0, 3, 1, 0, 1);
  for (size_t i = 0; i < 6; i++) SetStyle(h_del_phi_coaled_partons[i], ci[i], kFullCircle, 0.7, 1);
  for (size_t i = 0; i < 6; i++) SetStyle(h_CDelPhiX_meson_meson[i], ci[i], kFullCircle, 0.5, 1);
  for (size_t i = 0; i < 6; i++) SetStyle(h_CDelPhiP_meson_meson[i], ci[i], kFullCircle, 0.5, 1);
  for (size_t i = 0; i < 6; i++) SetStyle(h_distance[i], ci[i], kFullCircle, 0.5, 1);

  TLegend* legend_0 = new TLegend(0.6, 0.15, 0.85, 0.45);
  legend_0->SetName("num of partons");
  legend_0->AddEntry(h_CDelPhiX_meson_meson[0], "n = 10", "LEP");
  legend_0->AddEntry(h_CDelPhiP_meson_meson[1], "n = 20", "LEP");
  legend_0->AddEntry(h_CDelPhiX_meson_meson[2], "n = 30", "LEP");
  legend_0->AddEntry(h_CDelPhiP_meson_meson[3], "n = 60", "LEP");
  legend_0->AddEntry(h_CDelPhiP_meson_meson[4], "n = 100", "LEP");
  legend_0->AddEntry(h_CDelPhiP_meson_meson[5], "n = 300", "LEP");

  TCanvas* c_meson_x = new TCanvas("c_meson_x", "c_meson_x", 800, 400);
  c_meson_x->Divide(2);
  c_meson_x->cd(1);
  dummy->Draw("SAME");
  for (size_t i = 0; i < 6; i++) {
    h_CDelPhiX_meson_meson[i]->Draw("SAME");
    h_CDelPhiX_meson_meson[i]->Draw("SAME HIST LF2");
  }
  legend_0->Draw();
  TPad* pad = (TPad*)c_meson_x->GetPad(2);
  //pad->SetLogy(1);
  pad->cd();
  dummy2->Draw("SAME");
  for (size_t i = 0; i < 6; i++) {
    h_CDelPhiP_meson_meson[i]->Draw("SAME");
    h_CDelPhiP_meson_meson[i]->Draw("SAME HIST LF2");
  }

  TCanvas* c_dis = new TCanvas("c_dis", "c_dis", 400, 400);
  c_dis->cd();
  c_dis->SetLogy();
  //for (size_t i = 0; i < 1; i++) h_distance[i]->SetMaximum(50);
  for (size_t i = 0; i < 6; i++) h_distance[i]->Draw("SAME");
  for (size_t i = 0; i < 6; i++) h_distance[i]->Draw("SAME HIST LF2");
  legend_0->Draw();
}