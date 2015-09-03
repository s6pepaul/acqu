#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TChain.h"
#include <iostream>
#include <fstream>

using namespace std;

void CheckPol2( TString inputFile = "/scratch/simong/FileCheckPol.root", TString outFile = "/scratch/simong/FileCheckPol2.root" ){


  //Open input files
  TFile*  iFile   = new TFile(inputFile,"read"    );
  TFile*  oFile   = new TFile(outFile,  "recreate");

  TH2F*    hFileEG     = (TH2F*)iFile->Get("FileEG");

  TH1F*    hFileEGamor = (TH1F*)iFile->Get("FileEGamor");

  TH2F*    hFileEGenh  = (TH2F*)hFileEG->Clone("FileEGenh");

  hFileEGenh->Reset();

  Double_t ratio   = 0;
  Int_t    normBin = 215;

  //Get enhancement
  for(Int_t i=1;i<hFileEG->GetNbinsY();i++){

    if(hFileEG->GetBinContent(normBin,i)==0) continue;

    ratio  = hFileEGamor->GetBinContent(normBin)/hFileEG->GetBinContent(normBin,i);
    ratio += hFileEGamor->GetBinContent(normBin+1)/hFileEG->GetBinContent(normBin+1,i);
    ratio += hFileEGamor->GetBinContent(normBin+2)/hFileEG->GetBinContent(normBin+2,i);
    ratio += hFileEGamor->GetBinContent(normBin+3)/hFileEG->GetBinContent(normBin+3,i);
    ratio += hFileEGamor->GetBinContent(normBin+4)/hFileEG->GetBinContent(normBin+4,i);
    ratio += hFileEGamor->GetBinContent(normBin+5)/hFileEG->GetBinContent(normBin+5,i);

    ratio = ratio/6;

    for(Int_t j=1;j<hFileEG->GetNbinsX();j++){

      if(hFileEG->GetBinContent(j,i)==0) continue;
      if(hFileEGamor->GetBinContent(j)==0) continue;

      hFileEGenh->SetBinContent(j,i,(ratio*hFileEG->GetBinContent(j,i))/hFileEGamor->GetBinContent(j));

    }
  }


  hFileEG->Write();
  hFileEGamor->Write();
  hFileEGenh->Write();

  oFile->Close();
  iFile->Close();

  cout << endl;
  cout << inputFile << endl;
  cout << outFile << endl << endl;

}
