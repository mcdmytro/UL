#include "config.h"

void UP(){

  gSystem->Load("libRooFit");                             
  using namespace RooFit;  

  gROOT->ProcessLine(".x datafile.C"); 
  gROOT->ProcessLine(".x cutdata.C");
  gROOT->ProcessLine(".x X_states.C");

  gStyle->SetOptStat(1);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(0);
  // White BG
  gStyle->SetCanvasColor(10);
  // Format for axes
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetLabelSize(0.06,"xyz");
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetNdivisions(510,"xyz");
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleColor(1,"xyz");
  gStyle->SetTitleSize(0.07,"xyz");
  gStyle->SetTitleOffset(1.1,"xyz");
  
  // No pad borders
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadBorderSize(0);
  // White BG
  gStyle->SetPadColor(10);
  // Margins for labels etc.
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetPadBottomMargin(0.17);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.10);
  // No error bars in x direction
  gStyle->SetErrorX(0);
  // Format legend
  gStyle->SetLegendBorderSize(0);

  TH1F *h1=new TH1F("h1","", 30, 4, 11);
  TH1F *h2=new TH1F("h2","", 30, 4, 11);
  if(decay==2317){
    data.Draw("(m_dsp-m_dsi-m_ds17+1.9683+2.3178)>>h1",cut_fin);
    data.Draw("(m_dsp-m_dsi-m_ds17+1.9683+2.3178)>>h2",cut_fin);
  }
  if(decay==2460){
    data.Draw("(m_mp-m_dsw-m_ds60+1.9683+2.4595)>>h1",cut_fin);
    data.Draw("(m_mp-m_dsw-m_ds60+1.9683+2.4595)>>h2",cut_fin);
  }
  h2->Scale(0.5);
  h2->SetFillColor(kViolet);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(0.5);
  h1->SetMarkerColor(kBlack);
  h1->SetLineColor(kBlack);
  h1->Draw("E1");
  if(decay==2317){h1->GetXaxis()->SetTitle("M(D_{s}D_{s0}^{*}(2317))(GeV/c^{2})");}
  if(decay==2460){h1->GetXaxis()->SetTitle("M(D_{s}D_{s1}(2460))(GeV/c^{2})");}
  h1->GetYaxis()->SetTitle("Events/233 MeV/c^{2}");
  h1->GetYaxis()->CenterTitle();
  h1->GetXaxis()->CenterTitle();

  const double x_min= 4;
  const double x_max= 11;
  const double x_nbins= 30;

  RooRealVar x("x", "x", x_min, x_max);

  RooRealVar mean("mean", "Mean of gaussian", X_mass);
  RooRealVar sigma("sigma", "Sigma of gaussian", X_width);
  RooBreitWigner bw("bw", "Breit-Wigner for sig", x, mean, sigma);

  RooRealVar a0("a0", "a0", 0.5, 0., 1.);
  RooRealVar a1("a1", "a1", 0.5, 0., 1.);
  RooRealVar a2("a2", "a2", 0.5, 0., 1.);
  RooRealVar a3("a3", "a3", 0.5, 0., 1.);
  RooRealVar a4("a4", "a4", 0.5, 0., 1.);
  RooRealVar a5("a5", "a5", 0.5, 0., 1.);
  RooRealVar a6("a6", "a6", 0.5, 0., 1.);
  RooBernstein bernstein("bernstein", "Bernstein polynomial for bkg", x, RooArgList(a0, a1, a2, a3, a4, a5, a6));

  RooRealVar nsig("nsig", "number of chic1 signal events", 0, 0, 10);
  RooRealVar nbkg("nbkg", "number of background events", 100, 0, 200);

  RooAddPdf sum("sum","bkg+bw", RooArgList(bw,bernstein), RooArgList(nsig,nbkg));

  RooDataHist hdata("hdata","hdata",x,h1);
  // RooDataHist hdata("hdata","hdata",x,h1);
 
  Double_t _nll0;
  nsig.setVal(0.0);
  nsig.setConstant(kTRUE);
  rr = (RooFitResult*) sum.fitTo(hdata, Extended(kTRUE),Save());
  _nll0=rr->minNll();

  ofstream fout1("b_up_lik.dat");
  ofstream fout2("b_up_lik_origin.dat");
  
  Double_t n(0.0);
  Int_t nscan(11);
  Double_t _nll[1000], _br[1000];
  for (int i=0; n<nscan; )
  {

    nsig.setVal(n);
    nsig.setConstant(kTRUE);
    rr = (RooFitResult*) sum.fitTo(hdata,Extended(kTRUE),Save());
    _br[i] = nsig.getVal();
    _nll[i]=rr->minNll();
    char buffer1[256];
    char buffer2[256];
    sprintf(buffer1,"%6.4f",_br[i]);
    sprintf(buffer2,"%10.4f",n);
    fout1 << buffer1 <<"  " << exp(-1.0*(_nll[i]-_nll0)) << endl;
    fout2 << buffer1 <<"  " << _nll[i] << endl;
    n += 1;
  }
  nsig.setConstant(kFALSE);

  RooFitResult* rr = sum.fitTo(hdata, Extended(kTRUE),Save());
  sum.fitTo(hdata,Minos(kTRUE),Save(kTRUE), Hesse(kFALSE),Strategy(2), PrintLevel(1), Extended());
  RooPlot* xframe = x.frame();

  hdata.plotOn(xframe);
  sum.plotOn(xframe);

  sum.plotOn(xframe,Components(bernstein),LineStyle(kDashed));

  xframe->Draw("");
  rr->Print("v") ;
} 
