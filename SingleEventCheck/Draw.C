#include "TFile.h"
#include "TLine.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TEllipse.h"

void Draw() {
  int ci[4];
  TColor* color[4];
  ci[0] = TColor::GetFreeColorIndex();
  color[0] = new TColor(ci[0],   0/255.,  24/255., 113/255.);//dark blue
  ci[1] = TColor::GetFreeColorIndex();
  color[1] = new TColor(ci[1], 255/255.,  88/255.,  93/255.);//red
  ci[2] = TColor::GetFreeColorIndex();
  color[2] = new TColor(ci[2], 255/255., 181/255.,  73/255.);//yellow
  ci[3] = TColor::GetFreeColorIndex();
  color[3] = new TColor(ci[3], 65/255.,  182/255., 230/255.);//light blue
  gStyle->SetOptStat(false);
  TFile* f = TFile::Open("ACoaleEvent.root");

  TGraph* g_partons[2];
  g_partons[0] = (TGraph*)f->Get("partons");
  g_partons[1] = (TGraph*)f->Get("antipartons");
  g_partons[0] -> SetLineWidth(0);
  g_partons[1] -> SetLineWidth(0);
  g_partons[0] -> SetMarkerStyle(kFullCircle);
  g_partons[1] -> SetMarkerStyle(kOpenCircle);
  g_partons[0] -> SetMarkerSize(0.7);
  g_partons[1] -> SetMarkerSize(0.7);

  TH2D* dummy = new TH2D("","",1,-1.2,1.2,1,-1.2,1.2);
  TEllipse* circle = new TEllipse(0,0,1,1);
  TCanvas* cPartons = new TCanvas("cPartons","cPartons",400,400);
  dummy->Draw("SAME");
  circle->Draw("SAME");
  g_partons[0]->Draw("SAME LEP");
  g_partons[1]->Draw("SAME LEP");


  TCanvas* c_candi = new TCanvas("c_candi","c_candi",1200,400);
  c_candi->Divide(3);
  c_candi->cd(1);
  dummy->Draw("SAME");
  circle->Draw("SAME");
  g_partons[0]->Draw("SAME LEP");
  g_partons[1]->Draw("SAME LEP");
  TGraph* g_meson_candi[100000];
  TGraph* g_baryon_candi[100000];
  TGraph* g_antibaryon_candi[100000];
  for (size_t i = 0; i < 100000; i++) {
    g_meson_candi[i] = (TGraph*)f->Get(Form("meson_candi_%zu",i));
    if(g_meson_candi[i] != nullptr) g_meson_candi[i] ->SetLineColor(ci[2]);
    if(g_meson_candi[i] != nullptr) g_meson_candi[i] ->Draw("SAME LEP");
  }
  c_candi->cd(2);
  dummy->Draw("SAME");
  circle->Draw("SAME");
  for (size_t i = 0; i < 100000; i++){
    g_baryon_candi[i] = (TGraph*)f->Get(Form("baryon_candi_%zu",i));
    if(g_baryon_candi[i] != nullptr) g_baryon_candi[i] ->SetLineColor(ci[0]);
    if(g_baryon_candi[i] != nullptr) g_baryon_candi[i] ->Draw("SAME LEP");
  }
  c_candi->cd(3);
  dummy->Draw("SAME");
  circle->Draw("SAME");
  for (size_t i = 0; i < 100000; i++){
    g_antibaryon_candi[i] = (TGraph*)f->Get(Form("antibaryon_candi_%zu",i));
    if(g_antibaryon_candi[i] != nullptr) g_antibaryon_candi[i] ->SetLineColor(ci[1]);
    if(g_antibaryon_candi[i] != nullptr) g_antibaryon_candi[i] ->Draw("SAME LEP");
  }

  TCanvas* c_hardon = new TCanvas("c_hardon","c_hardon",400,400);
  c_hardon->cd();
  dummy->Draw("SAME");
  circle->Draw("SAME");
  g_partons[0]->Draw("SAME LEP");
  g_partons[1]->Draw("SAME LEP");
  TGraph* g_meson[100000];
  TGraph* g_baryon[100000];
  TGraph* g_antibaryon[100000];
  TGraph* g_meson_partons[100000];
  TGraph* g_baryon_partons[100000];
  TGraph* g_antibaryon_partons[100000];
  for (size_t i = 0; i < 100000; i++) {
    g_meson[i] = (TGraph*)f->Get(Form("meson_%zu",i));
    if(g_meson[i] != nullptr) g_meson[i] ->SetMarkerColor(ci[2]);
    if(g_meson[i] != nullptr) g_meson[i] ->SetMarkerStyle(kFullCircle);
    if(g_meson[i] != nullptr) g_meson[i] ->Draw("SAME LEP");

    g_meson_partons[i] = (TGraph*)f->Get(Form("meson_partons_%zu",i));
    if(g_meson_partons[i] != nullptr) g_meson_partons[i] ->SetLineColor(ci[2]);
    if(g_meson_partons[i] != nullptr) g_meson_partons[i] ->Draw("SAME LEP");
  }
  for (size_t i = 0; i < 100000; i++){
    g_baryon[i] = (TGraph*)f->Get(Form("baryon_%zu",i));
    if(g_baryon[i] != nullptr) g_baryon[i] ->SetMarkerColor(ci[0]);
    if(g_baryon[i] != nullptr) g_baryon[i] ->SetMarkerStyle(kFullCircle);
    if(g_baryon[i] != nullptr) g_baryon[i] ->Draw("SAME LEP");

    g_baryon_partons[i] = (TGraph*)f->Get(Form("baryon_partons_%zu",i));
    if(g_baryon_partons[i] != nullptr) g_baryon_partons[i] ->SetLineColor(ci[0]);
    if(g_baryon_partons[i] != nullptr) g_baryon_partons[i] ->Draw("SAME LEP");
  }
  for (size_t i = 0; i < 100000; i++){
    g_antibaryon[i] = (TGraph*)f->Get(Form("antibaryon_%zu",i));
    if(g_antibaryon[i] != nullptr) g_antibaryon[i] ->SetMarkerColor(ci[1]);
    if(g_antibaryon[i] != nullptr) g_antibaryon[i] ->SetMarkerStyle(kFullCircle);
    if(g_antibaryon[i] != nullptr) g_antibaryon[i] ->Draw("SAME LEP");

    g_antibaryon_partons[i] = (TGraph*)f->Get(Form("antibaryon_partons_%zu",i));
    if(g_antibaryon_partons[i] != nullptr) g_antibaryon_partons[i] ->SetLineColor(ci[1]);
    if(g_antibaryon_partons[i] != nullptr) g_antibaryon_partons[i] ->Draw("SAME LEP");
  }
}