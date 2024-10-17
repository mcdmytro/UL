#include "StandardHypoTestInvDemo.C"

void CountUL(double nobs, double s, double b, double sigmab){

  RooWorkspace w("w");
    
  // make Poisson model * Gaussian constraint
  w.factory("sum:nexp(s[0,0,5],b[1,0,5])");
  w.factory("Poisson:pdf(nobs[0,50],nexp)");
  w.factory("Gaussian:constraint(b0[0,10],b,sigmab[5])");
  w.factory("PROD:model(pdf,constraint)");

  w.var("b")->setVal(b);
  w.var("nobs")->setVal(nobs); 
  w.var("s")->setVal(s);


  w.var("b0")->setVal(0);
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

  mc.Print();
  w.import(mc);

  // make data set with the number of observed events
  RooDataSet dtst("dtst","", *w.var("nobs")); 

  dtst.add(*w.var("nobs") );
  w.import(dtst);

  w.Print();

  w.writeToFile("CountingModel.root", true);

  // optHTInv.noSystematics = true;

  StandardHypoTestInvDemo(
                          "CountingModel.root",   // file name
                          "w",                    // workspace name
                          "ModelConfig",          // S+B modelconfig name
                          "",                     // B model name
                          "dtst",                 // data set name
                          0,                      // calculator type (0 = Freq)
                          2,                      // test statistic type (2 = two-sided Profile Likelihood)
                          false,                  // use CLs
                          25,                     // number of points
                          0,                      // xmin
                          3,                      // xmax
                          1e4                     // number of toys
                          ); 
}

void CountingModel(){
  // CountUL(0, 0, 0.024, 0.088);
  CountUL(0, 0, 0.602, 0.088);
}