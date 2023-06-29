{

  gSystem->Load("libRooFit");                             
  using namespace RooFit;  

  gROOT->ProcessLine(".x datafile.C"); 
  gROOT->ProcessLine(".x cutdata.C");

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

  TH1F *h1=new TH1F("h1","", 300, 3.6, 4.2);
  data.Draw("mgfit01>>h1",cut06&&cut07&&cut08&&cut09&&cut19&&cut21);

  TH1F *h2=new TH1F("h2","", 300, 3.6, 4.2);
  data.Draw("mgfit01>>h2",cut06&&cut07&&cut08&&cut12&&cut19&&cut21);
  h2->Scale(0.5);
  h2->SetFillColor(kViolet);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(0.5);
  h1->SetMarkerColor(kBlack);
  h1->SetLineColor(kBlack);
  h1->Draw("E1");
  h1->GetXaxis()->SetTitle("MM(Z_{c}(3900))(GeV/c^{2})");
  h1->GetYaxis()->SetTitle("Events/1 MeV/c^{2}");
  h1->GetYaxis()->CenterTitle();
  h1->GetXaxis()->CenterTitle();

  const double x_min= 3.6;
  const double x_max= 4.2;
  const double x_nbins= 300;

  RooRealVar x("x", "x", x_min, x_max);
  RooRealVar sigmean("sigmean","mass", 3.88845e+00);//,3.7866,3.9866) ;
  RooRealVar sigwidth("sigwidth","width", 3.80486e-03);//,0.00,0.05) ;
  RooGaussian gauss("gauss","gaussian PDF",x,sigmean,sigwidth) ;
  RooRealVar width("width1","B-W Width", 2.91996e-02);//,0.00,0.1);
  RooVoigtian voigtian("voitian1","voitian1",x,sigmean,width,sigwidth);

  RooRealVar c0("c0","coefficient #0", 0,-10,10);
  RooRealVar c1("c1","coefficient #1", 0,-10,10);
  RooChebychev bkg("bkg", "background p.d.f.", x, RooArgList(c0));

  RooRealVar nsig("nsig","number of chic1 signal events",0,-100,100) ;
  RooRealVar nbkg("nbkg","number of background events",10,0,100) ;

  RooAddPdf sum("sum","bkg+bw",RooArgList(voigtian,bkg),RooArgList(nsig,nbkg)) ;

  RooDataHist hdata("hdata","hdata",x,h1);
 
  Double_t _nll0;
  nsig.setVal(0.0);
  nsig.setConstant(kTRUE);
  rr = (RooFitResult*) sum.fitTo(hdata, Extended(kTRUE),Save());
  _nll0=rr->minNll();

  ofstream fout1("b_up_lik.dat");
  ofstream fout2("b_up_lik_origin.dat");
  
  Double_t n(0.0);
  Int_t nscan(100);
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

  sum.plotOn(xframe,Components(bkg),LineStyle(kDashed));


  xframe->Draw("");
  rr->Print("v") ;
} 
