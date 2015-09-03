#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TChain.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <iostream>
#include <fstream>
#include <TMath.h>
#include <TF1.h>
#include <TMinuit.h>
#include <TCanvas.h>
#include <TSystem.h>

using namespace std;

Double_t GausOnBase(Double_t *, Double_t *);
void parFromHuman(Double_t beamMeV = 1500.0, Double_t edgeMeV = 750.0, Double_t spreadMeV = 20.0, 
		  Double_t colliDist_m = 2.5, Double_t colliRad_mm = 3.0, Int_t nVec = 2, Double_t *par=NULL);
void enhFromParams(Double_t *par=NULL);
Double_t efit(const Double_t *);


TCanvas *genCanvas=NULL;
TF1 *gausFit;

//Some enumerators and names
enum {      
  THETA,  // [0] theta      main angle responsible for coherent edge cutoffs
  SIGMA,  // [1] sigma      smearing of theta
  THETAR, // [2] thetar     relative angle resonsible for colli cutoffs
  SIGMAR, // [3] sigmar     smearing of colli cutoff angle
  E0MEV,  // [4] beam energy
  NVEC,   // [5] nvec       no of vectors contributing
  IVEC};  // [6] ivec[]     array of intensities of vectors up to nvec.

//Approx Form factor is F(g^2) = (q^2 + b^(-2)) ^ -2
//Where b= 111 x Z^(-1/3) (x 925 to get into units of crystal lattice)

const Double_t k=26.5601;        //put in formula for k later (my own stonehenge paper)

Double_t beamMeV = 1500.0;
Double_t edgeGuess = 750.0;
Double_t colliDist_m = 2.5;
Double_t colliRad_mm = 3.0;
Int_t nVec=2;
Double_t scalingN = 0.006;

const Int_t VECTORS[]={2,4,6,8,10};    //list of the vectors to be included (022,044);

Int_t THETASTEPS = 201;          //no of steps in convoluting with gaussian
Double_t LOWFIT = 400.0;         //how far below the peak in MeV to start fitting

Double_t fitMinEnergy;
Int_t fitMinBin;
Double_t fitMaxEnergy;
Int_t fitMaxBin;
Int_t verbose=0;
Double_t bestPar[10];
Double_t bestChisq;

Int_t counter = 1;
 
//these are really just convenient arrays - they don't ever get plotted.
TH1F *weightHist  = NULL;
TH2F *thetaWeight = NULL;
TH2F *thetaPol    = NULL;
TH2F *thetaTtot   = NULL;
TH2F *thetaItot   = NULL;

//All histograms and related are globals
TH1F *histP = NULL;  //pol
TH1F *histE = NULL;  //enh from calculation

TH1F *histD = NULL;  //enh from data to be fitted


//fit of a gaussian on a baseline, for gausFit,
Double_t GausOnBase(Double_t *x, Double_t *par) {
  Double_t arg = 0;
  if (par[2] != 0) arg = (x[0] - par[1])/par[2];
   Double_t fitval = par[3] + par[0]*TMath::Exp(-0.5*arg*arg);
   return fitval;
}

void CheckPol3( TString inputFile = "/home/pauli/acqu/acqu_user/data/fitEnhancement.root", TString outFile = "/home/pauli/acqu/acqu_user/data/EnhancementFitted.root" ){

  //Setup Canvas for checking stuff
  if(genCanvas) delete genCanvas;
  genCanvas = new TCanvas("genCanvas","genCanvas",50,50,800,1000);
  genCanvas->Divide(1,2);   
  genCanvas->GetPad(1)->SetGridx(1);
  genCanvas->GetPad(1)->SetGridy(1);
  genCanvas->GetPad(2)->SetGridx(1);
  genCanvas->GetPad(2)->SetGridy(1);

  //Open input files
  TFile*  iFile   = new TFile(inputFile,"read"    );
  TFile*  oFile   = new TFile(outFile,  "recreate");
  
  TH1F*    hFileE = (TH1F*)iFile->Get("LinPol_GatedCurrEnh");

  histD    = (TH1F*)hFileE->Clone("EnhacementData");
  histD->SetDirectory(0);
  
  UInt_t nBins = hFileE->GetNbinsX();
  Double_t diff1, diff2, lowmean,fitedge,scalefac;
  Double_t par[10];
  ROOT::Math::Minimizer* min;
  Char_t name[30];
  
  gausFit=new TF1("gausFit",GausOnBase,0,100,4);

  // find a reasonable minumum spot to set to 1 for the baseline.
  // the lowest 5 channel mean between 0.2 and 0.95 of the range
  lowmean=1000000.0;
  for(int n=(int)(0.15*(float)nBins);n<=(int)(0.95*(float)nBins);n++){
    if((histD->Integral(n-2,n+2)<lowmean)){
      lowmean=histD->Integral(n-2,n+2);
      cout << "lowmean: " << lowmean << endl;
    }
  }
  if(lowmean<1) lowmean=500.0;
  //histD->Scale(1.06*5.0/(lowmean));
    histD->Scale(5.0/(lowmean));

  //energy dependent scaling to improve fit
  for(UInt_t n=1;n<=nBins-1;n++){
    Double_t binenergy = histD->GetBinCenter(n);
    Double_t bintemp = histD->GetBinContent(n);
    if(binenergy==0) continue;
    bintemp = bintemp*TMath::Power((1/binenergy),scalingN);
    histD->SetBinContent(n,bintemp);
  }
  
  //Get rid of zeros
  for(UInt_t n=1;n<=nBins-1;n++){
    if(histD->GetBinContent(n)<0.1)histD->SetBinContent(n,histD->GetBinContent(n+1));
  }
  //Get rid of zeros 2nd pass
  for(UInt_t n=1;n<=nBins-1;n++){
    if(histD->GetBinContent(n)<0.1)histD->SetBinContent(n,histD->GetBinContent(n+1));
  }
  //  Get rid of spikes up and down
  for(UInt_t n=2;n<=nBins-1;n++){
    diff1=(histD->GetBinContent(n)-histD->GetBinContent(n-1))/histD->GetBinContent(n-1);
    diff2=(histD->GetBinContent(n)-histD->GetBinContent(n+1))/histD->GetBinContent(n+1);

    if (((fabs(diff1)>0.03)&&(fabs(diff2)>0.03))&&(fabs(diff1-diff2)<0.1)){
      histD->SetBinContent(n,0.5*(histD->GetBinContent(n-1)+histD->GetBinContent(n+1)));
    }
  }

  histD->SetMaximum(1.2*histD->GetMaximum());
  histD->SetMinimum(0.0);

  histE = (TH1F*)histD->Clone("histE");
  histE->SetDirectory(0);
  histE->Reset();
  histE->SetMaximum(1.2*histD->GetMaximum());
  histE->SetMinimum(0.0);
  histP = (TH1F*)histD->Clone("histP");
  histP->SetDirectory(0);
  histP->Reset();
  histP->SetMaximum(1);

  genCanvas->cd(1);
  histD->Draw("P");
  //Fit
  //Now try to make some guesses at the initial parameters
  histD->GetXaxis()->SetRange(100,300);
  gausFit->SetRange(histD->GetBinCenter(histD->GetMaximumBin()),histD->GetBinCenter(histD->GetMaximumBin())+100.0);
  gausFit->SetParameter(1,histD->GetBinCenter(histD->GetMaximumBin()));
  gausFit->SetParameter(2,10.0);
  gausFit->SetParameter(3,0.84);
  histD->Fit(gausFit,"rN");

  Int_t maxBin = histD->GetMaximumBin();

  histD->GetXaxis()->SetRange();

  //Get the edge from the fit
  fitedge=(gausFit->GetParameter(1)) + abs(gausFit->GetParameter(2));


  cout << beamMeV << endl;
  cout << fitedge << endl;
  cout << edgeGuess << endl;
  cout << gausFit->GetParameter(2) << endl;

  //Now we have enough information to set the basic parameters
  parFromHuman(beamMeV,fitedge,gausFit->GetParameter(2),colliDist_m,colliRad_mm,nVec,par);

  //set the intensities
  for(int v=0;v<par[NVEC];v++){                                               //give the vectors intensities
    par[IVEC+v] = histD->GetMaximum()*2.0/((Double_t)VECTORS[v]*(Double_t)VECTORS[v]);      //tailing off as 1/VECTORS[v]^2
  }
    

  enhFromParams(par);

  //Redo the intensities according to a the calc / data ration
  scalefac=histD->GetMaximum()/histE->GetMaximum();
  for(int v=0;v<par[NVEC];v++){                                               //give the vectors intensities
    par[IVEC+v]*=scalefac;
  }
  enhFromParams(par);
  histE->SetLineColor(2);
  histE->Draw("same");    
  genCanvas->cd(2);
  histP->SetLineColor(2);
  histP->Draw();
  genCanvas->Update();
  gSystem->ProcessEvents();
    
  //Set the range of the fit to be some sensible amount below peak and just past the 2nd peak.
  histD->GetXaxis()->SetRange(100,300);
  fitMinBin=histE->FindBin(histD->GetBinCenter(histD->GetMaximumBin())-LOWFIT);
  fitMaxBin=40+histE->FindBin(par[E0MEV]/((((2.0/4.0)*((par[E0MEV]/histD->GetBinCenter(histD->GetMaximumBin()))-1.0))+1.0)));
  cout << "fitMinBin " << fitMinBin << endl;
  cout << "fitMaxBin " << fitMaxBin << endl;
  histD->GetXaxis()->SetRange();
  min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simple");
    
  // set tolerance , etc...
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
  min->SetMaxIterations(10000);  // for GSL 
  min->SetTolerance(0.001);
  min->SetPrintLevel(1);
    

  ROOT::Math::Functor ft(&efit,IVEC+nVec);     
  min->SetFunction(ft);
    
  
  //Now set the variables
  min->SetLimitedVariable(THETA,   "Theta",   par[THETA],      par[THETA]/100.0,  0.95*par[THETA], 1.05*par[THETA]);
  min->SetLimitedVariable(SIGMA,   "Sigma",   2.5*par[SIGMA],  par[SIGMA]/100.0,  par[SIGMA],  5.0*par[SIGMA]);
  min->SetLimitedVariable(THETAR,  "Thetar",  2.0*par[THETAR],     par[THETAR]/100.0, 0.2*par[THETAR], 5.0*par[THETAR]);
  min->SetLimitedVariable(SIGMAR,  "Sigmar",  0.5*par[SIGMAR], par[SIGMAR]/100.0, 0.1*par[SIGMAR], 20.0*par[SIGMAR]);
  min->SetFixedVariable(E0MEV,     "E0MeV",   par[E0MEV]);  //no of vectors
  min->SetFixedVariable(NVEC,      "Nvec",    par[NVEC]);  //no of vectors
  for(int v=0;v<nVec;v++){
    sprintf(name,"Vec0%d%d", VECTORS[v],VECTORS[v]);
    min->SetVariable(v+IVEC, name, par[v+IVEC], par[v+IVEC]/100.0);
 
  }

  verbose=1;             //make it show the fit as it's improving
  bestChisq=100000.00;   //set this high for starters
    
  min->Minimize(); 
  
  //enhFromParams(bestPar);
 
  genCanvas->cd(1);
  histD->Draw("P");
  histE->SetLineColor(2);
  histE->Draw("same");    
  genCanvas->cd(2);
  histP->SetLineColor(2);
  histP->Draw();
  gSystem->Sleep(500);

  histD->SetDirectory(oFile);
  histD->Write();
  histE->SetDirectory(oFile);
  histE->Write();
  histP->SetDirectory(oFile);
  histP->Write();

  histD->SetDirectory(false);
  histE->SetDirectory(false);
  histP->SetDirectory(false);
  
  oFile->Close();
  iFile->Close();
  cout << endl;
  cout << inputFile << endl;
  cout << outFile << endl << endl;

}
  
void parFromHuman(Double_t beamMeV, Double_t edgeMeV, Double_t spreadMeV, Double_t colliDist_m, Double_t colliRad_mm, Int_t nVec, Double_t *par){

  //takes some physical quantities and makes them into parameters, then calls the 
  //enhFromParams function.
  
  
  //  Double_t par[10];                                                           //array of parameters
  Int_t g = 2;                                                                //variables used in CLAS note
  Double_t E0 = beamMeV;
  Double_t Eg = edgeMeV;
  
  
  par[THETA]  = k/(g*E0*E0*((1/Eg)-(1/E0)));                                  //theta from edge and beam energy
  par[SIGMA]  = (par[THETA]-(k/(g*E0*E0*((1/(Eg-spreadMeV))-(1/E0)))))/3.0;   //spread in theta from spread in edge 
  par[THETAR] = E0*0.001*5.0*colliRad_mm/colliDist_m;                         //cut from collimator
  par[SIGMAR] = par[THETAR]*par[SIGMA]/par[THETA];                            //smear in above same fractional sigma as above
  par[E0MEV]  = E0;                                                           //beam energy
  par[NVEC]   = (Double_t)nVec;                                                         //no of harmonics

  for(int v=0;v<par[NVEC];v++){                                               //give the vectors intensities
    par[IVEC+v] = 2.0/(Double_t)VECTORS[v];                                   //tailing off as 1/VECTORS[v]
    cout << IVEC+v << "  v   " << par[IVEC+v] << endl; 
  }
  cout << "Theta: " << par[THETA] << endl;
  cout << "Sigma: " << par[SIGMA] << endl;
  cout << "Thetar: " << par[THETAR] << endl;
  cout << "Sigmar: " << par[SIGMAR] << endl;
  cout << "E0: " << par[E0MEV] << endl;
  cout << "Nvec: " << par[NVEC] << endl;
}  


//The main customized fitting function which gets called by MINUIT
Double_t efit(const Double_t *parms){
  
  Double_t chisq = 1.0;
  Double_t delta;
  Double_t b1,b2;
  Double_t err;
  Double_t *par = (Double_t*)parms;

  histE->Reset("ICE"); //reset the histogram
  //  cout << "par[0]= " << par[0] << endl;

  //call the function to make the enhancement and polarization
  enhFromParams(par);

  chisq = 1.0;
  //loop over all the required bins in the histogram to work out a chisq
  for(int n=fitMinBin;n<=fitMaxBin;n++){
    b1=histE->GetBinContent(n);
    b2=histD->GetBinContent(n);
    err=1.0;
    delta=(b1-b2)/err;
    chisq+=(delta*delta);
    //note - not a proper chisq because its an enhancement
  }
   
  fprintf(stderr,"Chisq: \t%6.2f\t\r",chisq);


  if(chisq<bestChisq){
    bestChisq=chisq;
    for(int n=0;n<10;n++){
      bestPar[n]=par[n];
    }
    if(verbose){
      if(10%(counter++)){
	//if verbose, draw this on the canvas for every iteration to see how it's going
	genCanvas->cd(1);
	histD->SetLineColor(4);
	histD->Draw("P");
	histD->SetMinimum(0.5);
	histD->SetMaximum(2.8);
	genCanvas->cd(1);
	histE->Draw("same");
	
	cout << "HI" << endl;

	genCanvas->cd(2);
	histP->Draw();
	histP->SetMinimum(0);
	histP->SetMaximum(1);
	
	genCanvas->Draw();   
	
	genCanvas->Update();
	gSystem->ProcessEvents();
	counter=1;
      }
    }
  }
  return chisq;
  
}

void enhFromParams(Double_t *par){
  //make an enhancement and corresponding polarization from some the parameters as defined in the CLAS note.
  //this function is can be called stand alone, but will also ba called many times from the fitting function

  Double_t xd[10];
  Double_t xc[10];
  Double_t Q[10];
  Double_t cohContrib;
  Double_t cohTotal;
  Double_t phiTotal;
  Double_t etotal;
  Double_t ptotal;
  Double_t x=0.0;
  Int_t    g=0;
  Double_t weight=0.0;
  Double_t weightSum=0.0;
  Double_t polSum=0.0;
  Double_t phi,chi,cd;
  Double_t amo;
  Int_t jbin=0;
 
   
  if(!thetaPol){
    weightHist   = new TH1F("weightHist",  "weightHist", THETASTEPS+1, 0, THETASTEPS+1 );
    thetaWeight  = new TH2F("thetaWeight", "thetaWeight",histE->GetNbinsX(), histE->GetXaxis()->GetXbins()->GetArray(), THETASTEPS+1,0, THETASTEPS+1);
    thetaPol     = new TH2F("thetaPol",    "thetaPol",   histE->GetNbinsX(), histE->GetXaxis()->GetXbins()->GetArray(), THETASTEPS+1,0, THETASTEPS+1);
    thetaItot    = new TH2F("thetaItot",   "thetaItot",  histE->GetNbinsX(), histE->GetXaxis()->GetXbins()->GetArray(), THETASTEPS+1,0, THETASTEPS+1);
  }

  //reset them all for fresh filling
  histE->Reset("ICE");
  histP->Reset("ICE");
  thetaPol->Reset("ICE");
  thetaItot->Reset("ICE");
  weightHist->Reset("ICE");
  thetaWeight->Reset("ICE");

  for(Double_t j=par[THETA]-3.0*par[SIGMA];j<=par[THETA]+3.001*par[SIGMA];j+=(6.0*par[SIGMA])/THETASTEPS){
    
    weight=TMath::Gaus(j,par[THETA],par[SIGMA]);   //get the weight from the gaussian
    weightSum+=weight;                             //add to sum      
    
    //find the discontinuity for each vector
    for(int v=0;v<par[NVEC];v++){
      g=VECTORS[v];
      xd[v]=1.0/((k/(g*par[E0MEV]*j))+1.0);
      Q[v]=(1.0-xd[v])/xd[v];
      xc[v]=xd[v]/(1+((par[THETAR]*par[THETAR])*(1-xd[v])));
    }

    //loop over all bins in the histogram
    for(int bin=1;bin<=histE->GetNbinsX();bin++){
      x=histE->GetBinCenter(bin)/par[E0MEV];            //find the value of the bin
      amo=1/x;                                    //assume amo = inc = 1/x over regio of interest
      //amo=TMath::Power(x,-0.5);                                    //assume amo = inc = 1/x over regio of interest
      
      cohTotal=0.0;
      phiTotal=0.0;
      
      //loop over all the vectors
      for(int v=0;v<par[NVEC];v++){
	if(x>xd[v]) continue;           //only do up to x_dg
	 
	//work out chi and phi
	phi=(2*Q[v]*Q[v]*x*x)/((1-x)*(1+((1-x)*(1-x))-((4*Q[v]*Q[v]*x*x/(1-x))*(((1-x)/(Q[v]*x))-1))));
	chi=((Q[v]*Q[v]*x)/(1-x))*(1+((1-x)*(1-x))-((4*Q[v]*Q[v]*x*x/(1-x))*(((1-x)/(Q[v]*x))-1)));
	//	cout << j  << "  " << chi << endl;
	cd=0.5*(1+TMath::Erf((x-xc[v])/(TMath::Sqrt(2)*par[SIGMAR])));

	//get coherent contrib for the vector
	cohContrib=cd*par[IVEC+v]*chi;

	//add to the total and update the phi total
	cohTotal+=cohContrib;
	phiTotal+=cohContrib*phi;

      }
      if(cohTotal>0.0) {
	phiTotal/=cohTotal;   //divide by the cohTotal to get the weighted dmean phi
	//cout << x << " " << phiTotal << " " << cohTotal << " " << weight << endl;	 
      }

      //enhancement = coherent total + inc (or amo).
      etotal=(amo+cohTotal)/amo;
      //and pol like this
      //      ptotal=phiTotal*cohTotal/(cohTotal + amo);
      ptotal=phiTotal*cohTotal;

      //add the weighted contribution to the enhancement
      histE->Fill(x*par[E0MEV],weight*etotal);

      //keep the pol for this x,theta coord
      thetaPol->Fill(x*par[E0MEV],jbin,ptotal);

      //keep the total intensity for this x,theta coord
      thetaItot->Fill(x*par[E0MEV],jbin,cohTotal+amo);
    }
    
    //save the weight for this theta point
    weightHist->Fill(jbin,weight);
    jbin++;

  }
  //normalize the sum of the weighted enhancements
  histE->Scale(1.0/weightSum);
  
  
  //loop over each x bin, adding the weighted contribs from each theta pos
  for(int bin=1; bin<=histP->GetNbinsX(); bin++){
    weightSum=0.0;
    polSum=0.0;
    
    for(int jb=1;jb<=weightHist->GetNbinsX();jb++){
      weight=weightHist->GetBinContent(jb);

      //      polSum+=thetaPol->GetBinContent(bin,jb)*thetaItot->GetBinContent(bin,jb)*weight;
      polSum+=thetaPol->GetBinContent(bin,jb)*weight;
      weightSum+=thetaItot->GetBinContent(bin,jb)*weight;
      //polSum+=thetaPol->GetBinContent(bin,jb)*weight;
      //weightSum+=weight;
    }
    polSum/=weightSum;
    histP->Fill(histP->GetBinCenter(bin),polSum);
  } 
}
