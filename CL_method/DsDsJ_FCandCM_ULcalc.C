#include "StandardHypoTestInvDemo.C"
#include <tuple>
#include "Riostream.h"

using namespace RooFit;

// Performs fit and returns nobs, s, b, b_err
tuple <double, double, double, double> perform_fit(int name, double mass_mean, double width_mean, RooRealVar *mass, RooDataSet *data){
    printf("Name: X(%i), Mass: %.3f GeV, Width: %.3f Gev", name, mass_mean, width_mean);

    RooRealVar *mean  = new RooRealVar("mean", "Mean of gaussian", mass_mean);
    RooRealVar *sigma = new RooRealVar("sigma", "Sigma of gaussian",width_mean);
    RooBreitWigner *bw = new RooBreitWigner("bw","BW PDF",*mass, *mean, *sigma);

    RooRealVar *b0 = new RooRealVar("b00","b00",0.,1.);
    RooRealVar *b1 = new RooRealVar("b1","b1",0.,1.);
    RooRealVar *b2 = new RooRealVar("b2","b2",0.,1.);
    RooRealVar *b3 = new RooRealVar("b3","b3",0.,1.);
    RooRealVar *b4 = new RooRealVar("b4","b4",0.,1.);
    RooRealVar *b5 = new RooRealVar("b5","b5",0.,1.);
    RooRealVar *b6 = new RooRealVar("b6","b6",0.,1.);

    RooBernstein *bernstein = new RooBernstein("bernstein", "Bernstein polynomial for background",*mass,RooArgList(*b0,*b1,*b2,*b3,*b4,*b5,*b6));

    RooRealVar *sig  = new RooRealVar("N_sig",  "N true signal", 3, 0, 10);
    RooRealVar *bkg  = new RooRealVar("N_bkg",  "N background", 100, 10, 300);

    RooAddPdf *pdf = new RooAddPdf("pdf", "Gaussian + Pol",RooArgList(*bw, *bernstein), RooArgList(*sig, *bkg));

    RooFitResult *fitresult  = pdf->fitTo(*data, Save(kTRUE), Range("signal"), PrintEvalErrors(0));
    fitresult->Print("v");

    // Integral calculation
    Double_t
    CL90 = 1.64486,
    CL95 = 1.95996;

    Double_t
    I0 = mean->getVal()-CL90*(sigma->getVal()),
    I1 = mean->getVal()+CL90*(sigma->getVal());

    mass->setRange("int", I0, I1);

    RooAbsReal* int_sig = bw->createIntegral(*mass, Range("int"), NormSet(*mass));
    RooAbsReal* int_bkg = bernstein->createIntegral(*mass, Range("int"), NormSet(*mass));
    RooAbsReal* int_all = pdf->createIntegral(*mass, Range("int"), NormSet(*mass));

    Double_t
    int_val_scaled_sig = int_sig->getVal()*sig->getVal(),
    int_err_scaled_sig = int_sig->getVal()*sig->getError(),
    int_val_scaled_bkg = int_bkg->getVal()*bkg->getVal(),
    int_err_scaled_bkg = int_bkg->getVal()*bkg->getError(),
    int_val_scaled_all = int_all->getVal()*(sig->getVal()+bkg->getVal()),
    int_err_scaled_all = int_all->getVal()*(sig->getError()+bkg->getError());

    printf("\n--> Integral range:  [%.3f, %.3f]\n", I0, I1);
    printf("\tSgn integral value:\t %.3f +- %.3f \n", int_val_scaled_sig, int_err_scaled_sig);
    printf("\tBkg integral value:\t %.3f +- %.3f \n", int_val_scaled_bkg, int_err_scaled_bkg);
    printf("\tAll integral value:\t %.3f +- %.3f \n", int_val_scaled_all, int_err_scaled_all);

    // return int_val_scaled_all;
    return make_tuple(int_val_scaled_all, int_val_scaled_sig, int_val_scaled_bkg, int_err_scaled_bkg);
}

// Pepforms scan
double CountUL(int name, double nobs, double s, double b, double sigmab){
  RooWorkspace w("w");

  // make Poisson model * Gaussian constraint
  w.factory("sum:nexp(s[0,0,15],b[1,0,10])");
  w.factory("Poisson:pdf(nobs[0,50],nexp)");
  w.factory("Gaussian:constraint(b0[0,10],b,sigmab[5])");
  w.factory("PROD:model(pdf,constraint)");

  w.var("b")->setVal(b);
  w.var("nobs")->setVal(nobs);
  w.var("s")->setVal(s);

  w.var("b0")->setVal(b);
  w.var("b0")->setConstant(true);
  w.var("sigmab")->setVal(sigmab*b);

  ModelConfig mc("ModelConfig",&w);
  mc.SetPdf(*w.pdf("model"));
  mc.SetParametersOfInterest(*w.var("s"));
  mc.SetObservables(*w.var("nobs"));
  mc.SetNuisanceParameters(*w.var("b"));

  // these are needed for the hypothesis tests
  mc.SetSnapshot(*w.var("s"));
  mc.SetGlobalObservables(*w.var("b0"));

//   mc.Print();
  w.import(mc);

  // make data set with the number of observed events
  RooDataSet dtst("dtst","", *w.var("nobs"));

  dtst.add(*w.var("nobs") );
  w.import(dtst);

//   w.Print();
  TString fileName = Form("CountingModel_%i.root", name);
  w.writeToFile(fileName, true);

  optHTInv.noSystematics = true;

  StandardHypoTestInvDemo(
                          fileName,               // file name
                          "w",                    // workspace name
                          "ModelConfig",          // S+B modelconfig name
                          "",                     // B model name
                          "dtst",                 // data set name
                          0,                      // calculator type (0 = Freq)
                          2,                      // test statistic type (2 = two-sided Profile Likelihood)
                          true,                  // use CLs
                          25,                     // number of points
                          0,                      // xmin
                          5,                      // xmax
                          1e4                     // number of toys
                          );

    TFile f(Form("Freq_Cls+b_grid_ts2_CountingModel_%i.root",name));
    RooStats::HypoTestInverterResult* r;
    f.GetObject("result_s", r);
    double ul = r->UpperLimit();

    return ul;
}

// Performs FCs method and returns
tuple <double, double> FeldmanCousinUL(double ntot, double nbkg){
  // Feldman-Cousin trust intervas calculation
  TFeldmanCousins f;

  Double_t ul = f.CalculateUpperLimit(ntot, nbkg);
  Double_t ll = f.GetLowerLimit();

  return make_tuple(ul, ll);
}

// Set X_states parameters
tuple <double, double, double> SetXstatesParams(int X_name, int mass_err_fact, int width_err_fact){
    Double_t X_mass_mean, X_mass_err, X_mass, X_width_mean, X_width_err, X_width, LHCb_err;

    // Considered X states: X(4274), X(4685), X(4630), X(4500) and X(4700)
    if(X_name==4274){
        // X(4274)
        X_mass_mean = 4.294;
        X_mass_err  = 0.004;
        X_mass = X_mass_mean+X_mass_err*mass_err_fact;
        X_width_mean = 0.053;
        X_width_err  = 0.005;
        X_width = X_width_mean+X_width_err*width_err_fact;
        LHCb_err = 0.02658;
    }
    if(X_name==4685){
        // X(4685)
        X_mass_mean = 4.684;
        X_mass_err  = 0.007;
        X_mass = X_mass_mean+X_mass_err*mass_err_fact;
        X_width_mean = 0.126;
        X_width_err  = 0.015;
        X_width = X_width_mean+X_width_err*width_err_fact;
        LHCb_err = 0.05536;
    }
    if(X_name==4630){
        // X(4630)
        X_mass_mean = 4.626;
        X_mass_err  = 0.016;
        X_mass = X_mass_mean+X_mass_err*mass_err_fact;
        X_width_mean = 0.174;
        X_width_err  = 0.027;
        X_width = X_width_mean+X_width_err*width_err_fact;
        LHCb_err = 0.03847;
    }
    if(X_name==4500){
        // X(4500)
        X_mass_mean = 4.474;
        X_mass_err  = 0.003;
        X_mass = X_mass_mean+X_mass_err*mass_err_fact;
        X_width_mean = 0.077;
        X_width_err  = 0.006;
        X_width = X_width_mean+X_width_err*width_err_fact;
        LHCb_err = 0.02408;
    }
    if(X_name==4700){
        // X(4700)
        X_mass_mean = 4.694;
        X_mass_err  = 0.004;
        X_mass = X_mass_mean+X_mass_err*mass_err_fact;
        X_width_mean = 0.087;
        X_width_err  = 0.008;
        X_width = X_width_mean+X_width_err*width_err_fact;
        LHCb_err = 0.05400;
    }

    return make_tuple(X_mass, X_width, LHCb_err);
}

// Process parameters assignment, fit, run counting model, run FCs and ofstream results
double run_both_methods(int name, int mass_err_fact, int width_err_fact, RooRealVar *mass, RooDataSet *data, bool syst_scan){
    // Set X parameters
    double X_mass, X_width, LHCb_err;
    tie(X_mass, X_width, LHCb_err) = SetXstatesParams(name, mass_err_fact, width_err_fact);
    // Perform fit
    double nobs, s, b, sigmab;
    tie(nobs, s, b, sigmab)  = perform_fit(name, X_mass, X_width, mass, data);
    sigmab /= b;

    double tot_syst_err=0, glob_syst=0, effBrBr=0;
    if(name==4274 || name==4685){
        glob_syst = 0.0950;
        effBrBr = 2.890e-5;
        tot_syst_err = sqrt(pow(glob_syst,2)+pow(sigmab,2)+pow(LHCb_err, 2));}
    else if(name==4630 || name==4500 || name==4700){
        glob_syst = 0.1122;
        effBrBr = 1.30e-5;
        tot_syst_err = sqrt(pow(glob_syst,2)+pow(sigmab,2)+pow(LHCb_err, 2));}

    // Run counting model
    double CM_ul=0, CM_ul_ps=0, CM_ul_ns=0;

    CM_ul    = CountUL(name, int(nobs), s, b, sigmab);
    if(syst_scan){
        CM_ul_ps = CountUL(name, int(nobs), s, b*(1+tot_syst_err), sigmab);
        CM_ul_ns = CountUL(name, int(nobs), s, b*(1-tot_syst_err), sigmab);
    }

    printf("\n\n\n\t >>>>> X(%i) <<<<<", name);

    if(syst_scan){
        printf("\n\t\t===== Systematics =====\n");
        printf("Glob. sys. err.:  %.2f%%\n", glob_syst*100);
        printf("Stat. err.:\t  %.2f%%\n", sigmab*100);
        printf("LHCb. err.:\t  %.2f%%\n", LHCb_err*100);
        printf("Tot. syst. err.:  %.2f%%\n", tot_syst_err*100);
    }

    printf("\n--> Counting model result (at 90%% CL):");
    printf("\n\t Upper Limit (no syst.) = %.3f \n", CM_ul);
    if(syst_scan){
        printf("\t Upper Limit (+ syst)   = %.3f \n", CM_ul_ps);
        printf("\t Upper Limit (- syst)   = %.3f \n", CM_ul_ns);
    }
    // Run FCs methos
    double FC_ul=0, FC_ll=0, FC_ul_ps=0, FC_ll_ps=0, FC_ul_ns=0, FC_ll_ns=0;
    tie(FC_ul, FC_ll) = FeldmanCousinUL(nobs, b);
    if(syst_scan){
        tie(FC_ul_ps, FC_ll_ps) = FeldmanCousinUL(nobs, b*(1+tot_syst_err));
        tie(FC_ul_ns, FC_ll_ns) = FeldmanCousinUL(nobs, b*(1-tot_syst_err));
    }
    // cout << "--> FC confidence intervals (at the 90\% CL):\n";
    // cout << "1. w/o syst." << endl;
    // cout << "\tUpper Limit = " <<  FC_ul << endl;
    // cout << "\tLower Limit = " <<  FC_ll << endl;
    // if(syst_scan){
    //     cout << "2. + syst." << endl;
    //     cout << "\tUpper Limit = " <<  FC_ul_ps << endl;
    //     cout << "\tLower Limit = " <<  FC_ll_ps << endl;
    //     cout << "3. - syst." << endl;
    //     cout << "\tUpper Limit = " <<  FC_ul_ns << endl;
    //     cout << "\tLower Limit = " <<  FC_ll_ns << endl;
    // }
    // Cross-sections calculation
    Double_t lum = 980.;
    printf("--> Cross-sections:\n");
    // printf("\tFeldman-Cousins method:\t%.2f fb\n", max(FC_ul, max(FC_ul_ps, FC_ul_ns))/effBrBr/lum);
    printf("\tCounting model method:\t%.2f fb\n", max(CM_ul, max(CM_ul_ps, CM_ul_ns))/effBrBr/lum);

    ofstream res("results.txt", ios::app);
    if(!mass_err_fact && !width_err_fact){
        res << "\n\t >>>>> X(" << name << ") <<<<<\n";
        if(syst_scan){
            res << "===== Systematics =====\n";
            res << "Glob. sys. err.: " << glob_syst*100 << "\%\n";
            res << "Stat. err.: \t " << sigmab*100 << "\%\n";
            res << "LHCb. err.: \t " << LHCb_err*100 << "\%\n";
            res << "Tot. syst. err.: " << tot_syst_err*100 << "\%\n";
            res << "=======================\n";
        }
    }
    else res << "!!! mass_err_fact = " << mass_err_fact << ", width_err_fact = " << width_err_fact << " !!! \n";
    res << "--> Counting model result (at 90\% CL): \n";
    res << "\t\tUpper Limit (no syst.) = " << CM_ul << endl;
    if(syst_scan){
        res << "\t\tUpper Limit (+ syst.)  = " << CM_ul_ps << endl;
        res << "\t\tUpper Limit (- syst.)  = " << CM_ul_ns << endl;
    }
    // res << "--> FC confidence intervals (at the 90\% CL):\n";
    // res << "1. w/o syst." << endl;
    // res << "\tUpper Limit = " <<  FC_ul << endl;
    // res << "\tLower Limit = " <<  FC_ll << endl;
    // if(syst_scan){
    //     res << "2. + syst." << endl;
    //     res << "\tUpper Limit = " <<  FC_ul_ps << endl;
    //     res << "\tLower Limit = " <<  FC_ll_ps << endl;
    //     res << "3. - syst." << endl;
    //     res << "\tUpper Limit = " <<  FC_ul_ns << endl;
    //     res << "\tLower Limit = " <<  FC_ll_ns << endl;
    // }
    res << "--> Cross-sections:\n";
    // res << "\t\tFeldman-Cousins method:\t" << max(FC_ul, max(FC_ul_ps, FC_ul_ns))/effBrBr/lum << " fb\n";
    res << "\t\tCounting model method:\t"  << max(CM_ul, max(CM_ul_ps, CM_ul_ns))/effBrBr/lum << " fb\n";
    res.close();

    return CM_ul/effBrBr/lum;
}

// Main function
void DsDsJ_FCandCM_ULcalc(){

    bool LHCb_syst_scan = false;
    bool syst_scan = false;

    remove("results.txt");

    gROOT->Reset();
    gSystem->Load("libRooFit");

    TChain chain_2317("out");
    TChain chain_2460("out");
    chain_2317.Add("../../TMVA/Ds2317/4.Classification_NEW/root_files/Data_outTMVA.root");
    chain_2460.Add("../../TMVA/DsDs2460/root_files/Data_outTMVA.root");

    RooRealVar *mass_2317 = new RooRealVar("mass_2317", "", 4, 11);
    RooRealVar *mass_2460 = new RooRealVar("mass_2460", "", 4, 11);

    Float_t m_dsp_2317, m_dsi_2317, m_ds17_2317, mvaValue_2317;
    chain_2317.SetBranchAddress("m_dsp",   &m_dsp_2317);
    chain_2317.SetBranchAddress("m_ds17",  &m_ds17_2317);
    chain_2317.SetBranchAddress("m_dsi",   &m_dsi_2317);
    chain_2317.SetBranchAddress("mvaValue",&mvaValue_2317);

    Float_t m_mp_2460, m_ds60_2460, m_dsw_2460, mvaValue_2460;
    chain_2460.SetBranchAddress("m_mp",    &m_mp_2460);
    chain_2460.SetBranchAddress("m_ds60",  &m_ds60_2460);
    chain_2460.SetBranchAddress("m_dsw",   &m_dsw_2460);
    chain_2460.SetBranchAddress("mvaValue",&mvaValue_2460);

    RooDataSet *data_2317 = new RooDataSet("data_2317","", RooArgSet(*mass_2317));
    RooDataSet *data_2460 = new RooDataSet("data_2460","", RooArgSet(*mass_2460));

    for(int i=0; i<chain_2317.GetEntries(); i++){
        chain_2317.GetEntry(i);
        float m_var_2317 = m_dsp_2317-m_dsi_2317-m_ds17_2317+1.9683+2.3178;
        if(m_var_2317>4 && m_var_2317<11 && m_ds17_2317>2.305 && m_ds17_2317<2.329 && m_dsi_2317>1.955 && m_dsi_2317<1.979 && mvaValue_2317>0.2){
            mass_2317->setVal(m_var_2317);
            data_2317->add(RooArgSet(*mass_2317));
    }}
    for(int i=0; i<chain_2460.GetEntries(); i++){
        chain_2460.GetEntry(i);
        float m_var_2460  = m_mp_2460-m_dsw_2460-m_ds60_2460+1.9683+2.4595;
        if(m_var_2460>4 && m_var_2460<11 && m_ds60_2460>2.448 && m_ds60_2460<2.472 && m_dsw_2460>1.955 && m_dsw_2460<1.979 && mvaValue_2460>0.2){
            mass_2460->setVal(m_var_2460);
            data_2460->add(RooArgSet(*mass_2460));
        }
    }

    int X_names[5] = {4274, 4685, 4630, 4500, 4700};
    double cross_sect_err[5] = {0};

    for(int X_i=0; X_i<1; X_i++){
        double cross_sect_ul_0 = run_both_methods(X_names[X_i], 0, 0, mass_2317, data_2317, syst_scan);
        if(LHCb_syst_scan){
            double cross_sect_ul_i = run_both_methods(X_names[X_i], 1, 0, mass_2317, data_2317, syst_scan);
            cross_sect_err[X_i] += pow(cross_sect_ul_i-cross_sect_ul_0, 2);
            cross_sect_ul_i = run_both_methods(X_names[X_i], -1, 0, mass_2317, data_2317, syst_scan);
            cross_sect_err[X_i] += pow(cross_sect_ul_i-cross_sect_ul_0, 2);
            cross_sect_ul_i = run_both_methods(X_names[X_i], 0, 1, mass_2317, data_2317, syst_scan);
            cross_sect_err[X_i] += pow(cross_sect_ul_i-cross_sect_ul_0, 2);
            cross_sect_ul_i = run_both_methods(X_names[X_i], 0, -1, mass_2317, data_2317, syst_scan);
            cross_sect_err[X_i] += pow(cross_sect_ul_i-cross_sect_ul_0, 2);
            cross_sect_err[X_i] = sqrt(cross_sect_err[X_i])/cross_sect_ul_0;
            printf("----- Cross-section error for X(%i) = %.2f %% -----", X_names[X_i], cross_sect_err[X_i]*100);
            ofstream res("results.txt", ios::app);
            res << "----- Cross-section error for X(" << X_names[X_i] << "):\t" << cross_sect_err[X_i]*100 << "\% -----\n";
            res.close();
        }
    }
    for(int X_i=5; X_i<5; X_i++){
        double cross_sect_ul_0 = run_both_methods(X_names[X_i], 0, 0, mass_2460, data_2460, syst_scan);
        if(LHCb_syst_scan){
            double cross_sect_ul_i = run_both_methods(X_names[X_i], 1, 0, mass_2460, data_2460, syst_scan);
            cross_sect_err[X_i] += pow(cross_sect_ul_i-cross_sect_ul_0, 2);
            cross_sect_ul_i = run_both_methods(X_names[X_i], -1, 0, mass_2460, data_2460, syst_scan);
            cross_sect_err[X_i] += pow(cross_sect_ul_i-cross_sect_ul_0, 2);
            cross_sect_ul_i = run_both_methods(X_names[X_i], 0, 1, mass_2460, data_2460, syst_scan);
            cross_sect_err[X_i] += pow(cross_sect_ul_i-cross_sect_ul_0, 2);
            cross_sect_ul_i = run_both_methods(X_names[X_i], 0, -1, mass_2460, data_2460, syst_scan);
            cross_sect_err[X_i] += pow(cross_sect_ul_i-cross_sect_ul_0, 2);
            cross_sect_err[X_i] = sqrt(cross_sect_err[X_i])/cross_sect_ul_0;
            printf("----- Cross-section error for X(%i) = %.2f %% -----", X_names[X_i], cross_sect_err[X_i]*100);
            ofstream res("results.txt", ios::app);
            res << "----- Cross-section error for X(" << X_names[X_i] << "):\t" << cross_sect_err[X_i]*100 << "\% -----\n";
            res.close();
        }
    }
}