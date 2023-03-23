#include <iostream>
#include <fstream>
#include "TH2.h"
#include "TBox.h"
#include "TLine.h"
#include "TCanvas.h"

void TestSquare() {
  TCanvas* c1 = new TCanvas("c1", "c1", 1000, 1000);
  c1->cd();
  TH2D* dummy = new TH2D("dummy", "dummy", 1, -1.5, 1.5, 1, -1.5, 1.5);
  dummy->Draw("SAME");
  TBox* square = new TBox(-1., -1., 1., 1.);
  square->SetFillColor(0);
  square->SetLineWidth(1);
  square->SetLineColor(kBlack);
  square->Draw("L SAME");
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