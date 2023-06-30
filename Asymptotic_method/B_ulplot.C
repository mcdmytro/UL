#include "Riostream.h"
#include "TH1.h"
#include "config.h"

#define NN 1000			//max number of the lppt
#define CONVOLVE_FLAG 1		//0-no covolve / 1-covolved
#define E_SYS syst_err        	//if covolved, the systematic error
// #define E_SYS 0.04087        	//if covolved, the systematic error
#define N_CUT 5			//choose how many values will be cutted from n=0. *if covolved, the value of the first serval n is incorrect, we could cut them and then fit them relative points with other n values.
void B_ulplot() {
//   example of macro to read data from an ascii file and
//   create a root file with an histogram and an ntuple.
// the file basic.dat has 3 columns of float data
   ifstream in;
   in.open("b_up_lik.dat");	//input file
   ifstream ini;
   ini.open("b_up_num.dat");	//another in put file, do not use in thes code.

   Int_t sum = 0 ;
   Int_t n=6;
   Int_t n_event;
   Int_t nlines = 0;
   Int_t ii = 0 ;
   
   Int_t iiflag = 0;
   Double_t xflag = 0;
   Double_t yflag = 0;

   Double_t tx;
   Double_t ty;
   Double_t teyl;
   Double_t teyh;
   Double_t SUM;
   Double_t partsum;
   Double_t ratio;
   Double_t step = 0.5;


   Double_t x[NN] = {0};
   Double_t y[NN] = {0};
   Double_t yy[NN] = 0;
   Double_t exl[NN] = {0};
   Double_t exh[NN] = {0};
   Double_t eyl[NN] = {0};
   Double_t eyh[NN] = {0};

   

//   Float_t c,d;
//   Char_t e[4];
//   TTree *tree2 = new TTree("test","each tree");
//   TH1F *h2 = new TH1F("distribution","distribution of NO",33,0.5,33.5);

//Canvas
   TCanvas *c1 = new TCanvas("Canvas 1 ","A Simple Graph with error bars",10,10,800,640);
/*
   //Set the Branches of tree
      tree->Branch("npart",&a,"npart/I");
      tree->Branch("ntrack",&b,"ntrack/I");
      tree->Branch("px",&c,"px/F");
      tree->Branch("py",&d,"py/F");
      tree->Branch("ch",&e,"ch/C");
*/
//   tree1->Branch("sum",&sum,"sum/I");
//   tree2->Branch("each",&k,"each/I");
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//input
   SUM=0;
   partsum=0;
   ratio=0;
   for(ii=0;1;ii++)
   {
      in >>x[ii]>>y[ii]  ;
//      ini >>x[ii];     
 
      if(!in.good())break;
      exl[ii]=0;
      exh[ii]=0;
      eyl[ii]=0;
      eyh[ii]=0;
      
      if (nlines < 1) printf("                    x       y       eyl       eyh\n");
      if (nlines < 5) printf("line %d  :  x=%6.2f, y=%6.2f, eyl=%6.2f, eyh=%6.2f\n",ii+1,x[ii],y[ii],eyl[ii],eyh[ii]);

      nlines++;
      n_event++;
    }
   printf("total lines : %d \n",n_event);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~convolve the delta L~~~~~~~~~~~~~~~
   if(CONVOLVE_FLAG)// 1:convolve the result ; 0:no convolved result
   {
	  Double_t E_sys = E_SYS; // error;
	  for(int i=0; i<n_event; i++)
	  {
	    for(int j=1; j<n_event; j++)
		    {
		      Double_t Cof  = (x[j]-x[j-1])/(sqrt(2*3.141593)*E_sys*x[j]);
		      Double_t Expo = ((x[j]-x[i])**2)/(2*((E_sys*x[j])**2));
//		      if(j>10)
   //			{
            yy[i] += Cof*y[j]*TMath::Exp(-Expo);
   //			}
   //			if((i == 25)||(i == 24))
   //			{
   //				printf("y2[%d]%d = %6.2f\n",i,j,yy[i]);
   //				printf("cof = %6.2f ,Expo = %6.2f y2[%d] = %6.2f\n",Cof,Expo,i,yy[i]);
   //			}
		    }
	  printf("cof = %6.2f ,Expo = %6.2f y2[%d] = %6.2f\n",Cof,Expo,i,yy[i]);
	  }
	  for(int i=0; i<n_event; i++)
	  {
		      y[i] = yy[i];
//		      if(i<100)
//			{
		      printf("y[%d] = %6.2f\n",i,y[i]);
//			}
	  }	
    

//~~~~~~~~~~~~smooth~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for(i = N_CUT;i>0;i--)
    {
    	y[i-1] = y[i] - ((y[i+2]-y[i])/(x[i+2]-x[i]))*(x[i]-x[i-1]);
    }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for(ii=1;ii<n_event;ii++)
      {
      SUM+=(y[ii]+y[ii-1])/2*(x[ii]-x[ii-1]);
      }


    for(ii=0;ii<n_event;ii++)
   {
      if (ii>0) partsum+=(y[ii]+y[ii-1])/2*(x[ii]-x[ii-1]);
      ratio=partsum/SUM*100;
      if((ratio>90)&&(iiflag == 0))
      {
      iiflag = 1;
      xflag = x[ii-1]+(90 - yflag)*(x[ii] - x[ii-1])/(ratio - yflag);
      }
      yflag = ratio;
      printf("GE(%i) = %4.2f  ,  ratio = %6.4f \% \n",X_name, x[ii], ratio);

    }

      printf("GE(%i) = %4.3f  ,  ratio =  90.00\% \n", X_name, xflag);
      if(CONVOLVE_FLAG)
      {
      printf("convoloved result, systematic uncertainty = %6.2f\%, smooth number = %d  \n",E_sys*100,N_CUT);
      }
      else
      {
      printf("no-convoloved result  \n",);
      }

//   Double_t x[ii] = {0};
//   Double_t y[ii] = {0};
//   Double_t exl[ii] = {0};
//   Double_t exh[ii] = {0};
//   Double_t eyl[ii] = {0};
//   Double_t eyh[ii] = {0};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   in.close();
   ini.close();
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//draw
   gr = new TGraphAsymmErrors(n_event,x,y,exl,exh,eyl,eyh);
   gPad->SetLeftMargin(0.13);
   gPad->SetRightMargin(0.05);
   gPad->SetBottomMargin(0.14);
   gPad->SetTopMargin(0.05);
   gr->SetTitle("");
   gr->GetXaxis()->SetNdivisions(507);
   gr->GetXaxis()->SetLabelFont(132);
   gr->GetXaxis()->SetLabelSize(0.05);
   gr->GetXaxis()->SetTitleSize(0.07);
   gr->GetXaxis()->SetTitleFont(132);
   gr->GetXaxis()->SetTitleOffset(0.85);
   gr->GetYaxis()->SetNdivisions(507);
   gr->GetYaxis()->SetLabelFont(132);
   gr->GetYaxis()->SetLabelSize(0.05);
   gr->GetYaxis()->SetTitleSize(0.07);
   gr->GetYaxis()->SetTitleOffset(0.7);
   gr->GetYaxis()->SetTitleFont(132);
   // gr->SetMarkerColor(4);
   // gr->SetFillColor(kWhite);
   gr->SetLineColor(kBlue);
   gr->SetMarkerStyle(20);
   gr->SetMarkerSize(0.8);
   gr->SetLineWidth(4);
   gr->GetXaxis()->SetTitle("N_{sig}");
   gr->GetYaxis()->SetTitle("#delta L");
//   gr->GetYaxis()->SetRangeUser(0,1.5);
   gr->GetXaxis()->SetRangeUser(0.,50);
   gr->Draw("AC");

   c1->Update();

   Double_t eff = 6.3e-3;
   Double_t lum = 982.1;

   // Save results in a txt file 
   ofstream res("results.txt", ios::app);

   if(mass_err_fact==0 && width_err_fact==0){
      printf("\n--> X(%i) cross-section UL:\n\t>>>>> %.3f fb <<<<<\n", X_name, xflag/eff/lum);
      res << "\nDecay: Ds(" << decay << ")\t\tState: X(" << X_name << ")\n";
      res << "\t N^UL:\t" << xflag << "\n";
      res << "\t Cross-section:\t" << xflag/eff/lum << " fb\n";
   }
   if(mass_err_fact!=0 && width_err_fact==0){
      printf("\t Mass variation %i sigma:  UL = %.3f fb\n", mass_err_fact, xflag/eff/lum);
      res << "\t Mass  variation " << mass_err_fact << "sigma:\t UL = " << xflag/eff/lum << "\n";
   }
   if(mass_err_fact==0 && width_err_fact!=0){
      printf("\t Width variation %i sigma:  UL = %.3f fb\n", width_err_fact, xflag/eff/lum);
      res << "\t Width variation " << width_err_fact << "sigma:\t UL = " << xflag/eff/lum << "\n";
   }
   if(mass_err_fact!=0 && width_err_fact!=0)
      printf("\t ERR: You can't vary mass and width simultaneously! \n");

   res.close();

   ofstream res_short("results_short.txt", ios::app);
   res_short << xflag/eff/lum << "\n";
   res_short.close();
}

