// Code to fit an enhancement spectrum
#include <stdio.h>
#include <iostream> 
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TH1.h>
#include <TH2.h>
#include <TLine.h>
#include <TSystem.h>
#include <TStyle.h>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TGraphErrors.h"
#include "TError.h"
#include "TVirtualFitter.h"
using namespace std;

Double_t GausOnBase(Double_t *, Double_t *);
Double_t efit(const Double_t *);
void fitData(Char_t *fname, Char_t *outfile = NULL, Double_t beamMeV = 1557.0, Double_t colliDist_m = 2.5,Double_t colliRad_mm = 1.0, Int_t nVec=2);
void drawAnother(Double_t beamMeV = 1557.0, Double_t edgeMeV = 350.0, Double_t spreadMeV = 10.0, 
		 Double_t colliDist_m = 2.5, Double_t colliRad_mm = 1.0, Int_t nVec = 3, Char_t *opt="");
void enhFromHuman(Double_t beamMeV = 1557.0, Double_t edgeMeV = 350.0, Double_t spreadMeV = 10.0, 
		  Double_t colliDist_m = 2.5, Double_t colliRad_mm = 1.0, Int_t nVec = 3);
void parFromHuman(Double_t beamMeV = 1557.0, Double_t edgeMeV = 350.0, Double_t spreadMeV = 10.0, 
		  Double_t colliDist_m = 2.5, Double_t colliRad_mm = 1.0, Int_t nVec = 3, Double_t *par=NULL);
void enhFromParams(Double_t *par=NULL);


//Some enumerators and names
enum {      
  THETA,  // [0] theta      main angle responsible for coherent edge cutoffs
  SIGMA,  // [1] sigma      smearing of theta
  THETAR, // [2] thetar     relative angle resonsible for colli cutoffs
  SIGMAR, // [3] sigmar     smearing of colli cutoff angle
  E0MEV,  // [4] beam energy
  NVEC,   // [5] nvec       no of vectors contributing
  IVEC};  // [6] ivec[]     array of intensities of vectors up to nvec.


// Some basic consts etc first
// Consts are all UPPER CASE


//Approx Form factor is F(g^2) = (q^2 + b^(-2)) ^ -2
//Where b= 111 x Z^(-1/3) (x 925 to get into units of crystal lattice)
const Double_t B = 0.247892436;  //where did I get that ? Timm ?
const Double_t A=0.03;           //made up for now, need to get the actual no for this later
const Double_t k=26.5601;        //put in formula for k later (my own stonehenge paper)

const Int_t VECTORS[]={2,4,6,8,10};    //list of the vectors to be included (022,044);

Int_t counter=1;

//THESE NEED TO BE CHANGED FOR EACH SETTING (ie comment in/out)

Int_t THETASTEPS = 201;          //no of steps in convoluting with gaussian
Double_t LOWFIT = 300.0;         //how far below the peak in MeV to start fitting

//All histograms and related are globals
TH1F *histP = NULL;  //pol
TH1F *histE = NULL;  //enh from calculation

TH1F *histD = NULL;  //enh from data to be fitted

 
//these are really just convenient arrays - they don't ever get plotted.
TH1F *weightHist  = NULL;
TH2F *thetaWeight = NULL;
TH2F *thetaPol    = NULL;
TH2F *thetaTtot   = NULL;
TH2F *thetaItot   = NULL;

TH1F *histEArray[100]; //to save a load of hists for later use
TH1F *histPArray[100];
TH1F *histDArray[100];
Double_t edges[100];
Double_t allparamsg[11][100]; 
Double_t allparams[100][11]; 
TGraph *allGraphs[10];

Int_t histCount=0;
 
Double_t energy[1000];
Double_t energyBins[1000];
Int_t nBins=0;

Double_t fitMinEnergy;
Int_t fitMinBin;
Double_t fitMaxEnergy;
Int_t fitMaxBin;
Int_t verbose=0;
Double_t bestPar[10];
Double_t bestChisq;

TF1 *gausFit;

TCanvas *genCanvas=NULL;
TCanvas *graphCanvas=NULL;

//general one off init things
Bool_t isInit=kFALSE;

void init(){
  gStyle->SetOptStat(kFALSE);

  if(genCanvas) delete genCanvas;
  genCanvas = new TCanvas("genCanvas","genCanvas",50,50,800,1000);
  genCanvas->Divide(1,2);   
  genCanvas->GetPad(1)->SetGridx(1);
  genCanvas->GetPad(1)->SetGridy(1);
  genCanvas->GetPad(2)->SetGridx(1);
  genCanvas->GetPad(2)->SetGridy(1);

  if(isInit) return;
  gausFit=new TF1("gausFit",GausOnBase,0,100,4);
  isInit=kTRUE;
}

//fit of a gaussian on a baseline, for gausFit,
Double_t GausOnBase(Double_t *x, Double_t *par) {
  Double_t arg = 0;
  if (par[2] != 0) arg = (x[0] - par[1])/par[2];
   Double_t fitval = par[3] + par[0]*TMath::Exp(-0.5*arg*arg);
   return fitval;
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
	histD->SetMinimum(0.0);
	histD->SetMaximum(6.0);
	genCanvas->cd(1);
	histE->Draw("same");
	
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

void fitData(Char_t *fname, Char_t *outfile, Double_t beamMeV, Double_t colliDist_m, Double_t colliRad_mm, Int_t nVec){
  //read in a file to the global data histogram then fit it
  //file should be of the form Egamma, Enhancement
  //if file begins -L if will be interpreted as a list of files eg -Lmyfilelist
  
  FILE *fp=NULL;
  Double_t enhancement[1000];
  Char_t line[250];
  Double_t diff1,diff2,fitedge,scalefac;
  Char_t name[30];
  Double_t lowmean=100000000.0;
  Double_t par[10];
  ROOT::Math::Minimizer* min;
  Char_t filename[100][200]; // allow up to 100 files in a list
  Int_t nFiles=0;
  Char_t tempname[200];
  Char_t plane[8];
  Double_t theta;
  Double_t lowPolLimit;
  Double_t highPolLimit;
  Double_t pol;
  
  //if(!isInit) init(); // init anything if needed
  init();
  //open all files and make enhancements, then fit them 
  if(strstr(fname,"-L")){   // check if it's a single or a list of files
    if((fp=fopen(fname+2,"r"))==NULL){
      cout << "Error, couldn't open " << fname+2 << endl;
      return;
    }
    while(fgets(line,200,fp)){
      if((line[0] == '*')||(line[0] == '#')) continue; //skip comments    
      sscanf(line,"%s",filename[nFiles++]);
    }
    fclose(fp);
  }
  else{                    // single file only
    strcpy(filename[0],fname);
    nFiles=1;
  }

  for(int f=0;f<nFiles;f++){                         //go over the all data files    
    if((fp=fopen(filename[f],"r"))==NULL){
      cout << "Error, couldn't open " << filename[f] << endl;
      continue;
    }
    cout << "Opened " << filename[f] << endl;

    nBins=0;
    while(fgets(line,200,fp)){                         //read all the data points
      if((line[0] == '*')||(line[0] == '#')) continue; //skip comments    
      sscanf(line,"%*d%lg%lg",&energy[nBins],&enhancement[nBins]);
      nBins++;
    }
    fclose(fp);
  
    if(f==0){       //only do this for the 1st file in the list

      //make some bin edges half way between the energy values.
      for(int b=0;b<nBins-1;b++){
	energyBins[b+1]=0.5*(energy[b]+energy[b+1]);
      }
      //and the top and bottom have width of the adjacent bin
      energyBins[0]     = energyBins[1]       - (energyBins[2]      - energyBins[1]);
      energyBins[nBins] = energyBins[nBins-1] + (energyBins[nBins-1]- energyBins[nBins-2]);

      //Now make any required histograms if neccessary
      histE        = new TH1F("Enhancement", "Enhancement;Energy (MeV)",nBins,energyBins);
      histP        = new TH1F("Polarization", "Polarization;Energy (MeV);Deg. of Pol.",nBins,energyBins);
      histD        = new TH1F("EnhancementData", "EnhancementData;Energy (MeV)",nBins,energyBins);
      
      histE->SetMinimum(0);
      histD->SetMinimum(0);
      histP->SetMinimum(0);
      histP->SetMaximum(1);
      histD->SetMarkerStyle(20);
      histD->SetMarkerSize(0.7);
      histD->GetXaxis()->SetNdivisions(20);
      histE->GetXaxis()->SetNdivisions(20);
      histP->GetXaxis()->SetNdivisions(20);
      histP->GetYaxis()->SetNdivisions(20);
    }
    
    //fill the histD with the enhancement
    for(int n=0;n<nBins;n++){
      histD->SetBinContent(n+1,enhancement[n]);
    }
    //Get rid of zeros
    for(int n=1;n<=nBins-1;n++){
      if(histD->GetBinContent(n)<0.1)histD->SetBinContent(n,histD->GetBinContent(n+1));
    }
    //Get rid of zeros 2nd pass
    for(int n=1;n<=nBins-1;n++){
      if(histD->GetBinContent(n)<0.1)histD->SetBinContent(n,histD->GetBinContent(n+1));
    }
    //  Get rid of spikes up and down
    for(int n=2;n<=nBins-1;n++){
      diff1=(histD->GetBinContent(n)-histD->GetBinContent(n-1))/histD->GetBinContent(n-1);
      diff2=(histD->GetBinContent(n)-histD->GetBinContent(n+1))/histD->GetBinContent(n+1);
      //    cout << histD->GetBinCenter(n) << " " << histD->GetBinContent(n) << " " << diff1 << " " << diff2 << endl;
      if (((fabs(diff1)>0.1)&&(fabs(diff2)>0.1))&&(fabs(diff1-diff2)<0.1)){
	//cout << "****" << endl;
	histD->SetBinContent(n,0.5*(histD->GetBinContent(n-1)+histD->GetBinContent(n+1)));
      }
    }
    
    // find a reasonable minumum spot to set to 1 for the baseline.
    // the lowest 5 channel mean between 0.2 and 0.95 of the range
    lowmean=1000000.0;
    for(int n=(int)(0.05*(float)nBins);n<=(int)(0.95*(float)nBins);n++){
      if((histD->Integral(n-2,n+2)<lowmean)){
	lowmean=histD->Integral(n-2,n+2);
	//cout << lowmean << endl;
      }
    }
    histD->Scale(5.0/(lowmean));
    genCanvas->cd(1);
    
    if(f==0){
      histD->SetMaximum(1.2*histD->GetMaximum());
      histD->SetMinimum(0.0);
    }
    
    histD->Draw("P");
    //Now try to make some guesses at the initial parameters
    gausFit->SetRange(histD->GetBinCenter(histD->GetMaximumBin()),histD->GetBinCenter(histD->GetMaximumBin())+100.0);
    gausFit->SetParameter(1,histD->GetBinCenter(histD->GetMaximumBin()));
    gausFit->SetParameter(2,10.0);
    gausFit->SetParameter(3,1.0);
    histD->Fit(gausFit,"rN");
    
    lowmean=0.0;
    //Get the edge from the derivative
    for(float d = histD->GetBinCenter(histD->GetMaximumBin());d < histD->GetBinCenter(histD->GetMaximumBin()+90.0);d+=0.1){ 
      if(gausFit->Derivative(d)<lowmean){
	lowmean=gausFit->Derivative(d);
	fitedge=d;
      }
    }
    //  cout << "edge = " << fitedge << " MeV" << endl;
    
    //Now we have enough information to set the basic parameters
    parFromHuman(beamMeV,fitedge,gausFit->GetParameter(2),colliDist_m,colliRad_mm,nVec,par);
    
    //set the intensities
    for(int v=0;v<par[NVEC];v++){                                               //give the vectors intensities
      par[IVEC+v] = histD->GetMaximum()*2.0/((Double_t)VECTORS[v]*(Double_t)VECTORS[v]);      //tailing off as 1/VECTORS[v]^2
      //cout << IVEC+v << "  v   " << par[IVEC+v] << endl; 
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
    gSystem->Sleep(500);

    //Set the range of the fit to be some sensible amount below peak and just past the 2nd peak.
    fitMinBin=histE->FindBin(histD->GetBinCenter(histD->GetMaximumBin())-LOWFIT);

    fitMaxBin    = histE->FindBin(par[E0MEV]/((((2.0/4.0)*((par[E0MEV]/histD->GetBinCenter(histD->GetMaximumBin()))-1.0))+1.0)));
    cout << "fitMaxBin " << fitMaxBin << endl;
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
    for(int n=0;n<nVec;n++){
      sprintf(name,"Vec0%d%d", VECTORS[n],VECTORS[n]);
      min->SetVariable(n+IVEC, name, par[n+IVEC], par[n+IVEC]/100.0); 
    }
    
    verbose=1;             //make it show the fit as it's improving
    bestChisq=100000.00;   //set this high for starters
    
    min->Minimize(); 
    
    genCanvas->cd(1);
    enhFromParams(bestPar);
    histD->Draw("P");
    histE->SetLineColor(2);
    histE->Draw("same");
    
    genCanvas->cd(2);
    histP->SetLineColor(2);
    histP->Draw();
    //cout << "1bestChisq: " << bestChisq << endl;

    //now save everything for
    edges[f]=fitedge;
    for(int n=0;n<10;n++){
      allparamsg[n][f]=bestPar[n];
      allparams[f][n]=bestPar[n];
      cout << bestPar[n] << endl;
}
    allparamsg[10][f]=bestChisq;
    allparams[f][10]=bestChisq;

    cout << par[E0MEV] << endl;

    histEArray[f]=(TH1F*)histE->Clone();
    histDArray[f]=(TH1F*)histD->Clone();
    histPArray[f]=(TH1F*)histP->Clone();
  }

  //Now we have save the data, fitted enh, and correponding pol, fitted edge and parameters, and chisq
  //The world is our oyster


  // if(nFiles>=4){
  //   //make graphs of the essential parameters with polynomials vs edge;
  //   //.... which we may need to use later.
  //   graphCanvas = new TCanvas("ParamFits","ParamFits",200,50,800,1200);
  //   graphCanvas->Divide(2,5);
  //   graphCanvas->Draw();
  //   for(int g=0;g<IVEC+nVec;g++){
  //     allGraphs[g]=new TGraph(nFiles,edges,allparamsg[g]);
  //     graphCanvas->cd(g+1);
  //     allGraphs[g]->Draw("ALP");
  //     allGraphs[g]->Fit("pol2");
  //   }
  // }

  //Now do some stuff to record things
  //Save the canvas as .root and .gif, then write out the parameters and table

  //  return;



  for(int f=0;f<nFiles;f++){
    sprintf(tempname,"%s.C",filename[f]);
    genCanvas->SaveAs(tempname);
    sprintf(tempname,"%s.gif",filename[f]);
    genCanvas->SaveAs(tempname);
    
    
    //assume that the plane is in the file title, and grab it
    sprintf(plane,"");
    if(strcasestr(filename[f],"para")) sprintf(plane,"Para"); 
    if(strcasestr(filename[f],"perp")) sprintf(plane,"Perp"); 
    if(strlen(plane)<2){
      fprintf(stdout,"Warning: Can't tell from the name whether file %s is for Para or Perp. Skipping table\n",filename[f]);
      continue;
    }

    if(outfile==NULL){
      sprintf(tempname,"%s.table",filename[f]);
      fp=fopen(tempname,"w");
    }
      
    else{
      //open existing file in append more
      fp=fopen(outfile,"a");
    }
    

    if(outfile==NULL){
      
      //First write all the photon energies
      fprintf(fp,"#Photon energies (MeV)\n");
      for(int n=0;n<nBins;n++){
	fprintf(fp,"PhotonEnergy:\t%6.2f\n",energy[n]);
      }
    }
    
    //Now write comments,  parameters etc
    fprintf(fp,"#Pol table based on enhancement in %s\n",filename[f]);
    fprintf(fp,"#Edge = %6.2f MeV, Beam Energy = %6.2f\n",edges[f],allparams[f][E0MEV]);
    
    //Say what plane it is
    fprintf(fp,"Plane:\t%s\n",plane);


    //State parameters used
    fprintf(fp,"Parameters: ");
    for(int n=0;n<11;n++){
      fprintf(fp,"%f ",allparams[f][n]);
    }
    fprintf(fp,"\n");
    
    //Print out errors on parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fprintf(fp,"#Parameter errors: ");
    for(int n=0;n<11;n++){
      fprintf(fp,"%f ",allparams[f][n]);
    }
    fprintf(fp,"\n");
    
    //Now we take a range in theta which corresponds to the edge position +/- 10%  
    //
    enhFromParams(allparams[f]);                   //make the central one first
    //and work out the range to trust (LOWFIT to edge of 2nd harmonic)
    lowPolLimit=histE->GetBinCenter(histE->GetMaximumBin())-LOWFIT;
    highPolLimit=allparams[f][E0MEV]/((((2.0/4.0)*((allparams[f][E0MEV]/histD->GetBinCenter(histD->GetMaximumBin()))-1.0))+1.0));
    
    theta=allparams[THETA][f];
    for(Double_t t=0.9*theta;t<1.1*theta;t+=0.01*theta){   //in 0.01 steps
      cout << t << endl;
      allparams[f][THETA]=t;                //change only the theta param
      enhFromParams(allparams[f]);          //and make the enhancement / pol
      
      //work out the edege for each by fitting with the gaus
      gausFit->SetRange(histE->GetBinCenter(histE->GetMaximumBin()),histE->GetBinCenter(histE->GetMaximumBin())+100.0);
      gausFit->SetParameter(1,histE->GetBinCenter(histE->GetMaximumBin()));
      gausFit->SetParameter(2,10.0);
      gausFit->SetParameter(3,1.0);
      histE->Fit(gausFit,"rQ");
      
      genCanvas->cd(1);
      histD->SetLineColor(4);
      histD->Draw("P");
      histD->SetMinimum(0.0);
      histD->SetMaximum(6.0);
      genCanvas->cd(1);
      histE->Draw("same");
      
      genCanvas->cd(2);
      histP->Draw();
      histP->SetMinimum(0);
      histP->SetMaximum(1);
      
      genCanvas->Draw();   
      
      genCanvas->Update();
      gSystem->ProcessEvents();
      
      lowmean=0.0;
      //Get the edge from the derivative
      for(float d = histE->GetBinCenter(histE->GetMaximumBin());d < histE->GetBinCenter(histE->GetMaximumBin()+90.0);d+=0.1){ 
	if(gausFit->Derivative(d)<lowmean){
	  lowmean=gausFit->Derivative(d);
	  fitedge=d;
	}
      }
      fprintf(fp,"#Coherent Edge position (MeV)\n");
      fprintf(fp,"Edge:\t%6.2f\n",fitedge);
      fprintf(fp,"#Polarisation for each energy bin\n");
      for(int b=0;b<nBins;b++){                                             //go over all the bins
	if((energy[b]<lowPolLimit)||(energy[b]>highPolLimit)) pol=-1.0; //-1 if out of range
	else(pol=histP->GetBinContent(histP->FindBin(energy[b])));      //otherwise, pol
	fprintf(fp,"Pol:\t%5.2f\n",pol);
      }
    }
  }
}


void drawAnother(Double_t beamMeV, Double_t edgeMeV, Double_t spreadMeV, 
		 Double_t colliDist_m, Double_t colliRad_mm, Int_t nVec, Char_t *opt){

  if(!isInit) init(); // init anything if needed
  
  
  //call the calculation with the next histogram in the list
  enhFromHuman(beamMeV,edgeMeV,spreadMeV,colliDist_m,colliRad_mm,nVec);
  
  histEArray[histCount] = (TH1F*)histE->Clone();
  histPArray[histCount] = (TH1F*)histP->Clone();

  genCanvas->Draw();
  genCanvas->cd(1);
  histEArray[histCount]->SetMinimum(0.0);
  histEArray[histCount]->Draw(opt);
  
  genCanvas->cd(2);
  histPArray[histCount]->Draw(opt);

  histCount++;
}

void enhFromHuman(Double_t beamMeV, Double_t edgeMeV, Double_t spreadMeV, Double_t colliDist_m, Double_t colliRad_mm, Int_t nVec){
  Double_t par[10];
  parFromHuman(beamMeV,edgeMeV,spreadMeV,colliDist_m,colliRad_mm,nVec,par);
  enhFromParams(par);                                                         //call the function to make the enh and pol

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
    //cout << IVEC+v << "  v   " << par[IVEC+v] << endl; 
  }
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
 
  //loop over sigma
  // for(int p=0;p<10;p++){
  //  cout << p << ": " << par[p] << ", ";
  //}
  //cout << endl;

  // if needed, make some hists
  if(!histE){
    histE        = new TH1F("Enhancement", "Enhancement",  1000, 0, par[E0MEV]);
    histP        = new TH1F("Polarization", "Polarization",1000, 0, par[E0MEV]);
    histE->SetMinimum(0);
    histP->SetMinimum(0);
    histP->SetMaximum(1);
  }    
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
