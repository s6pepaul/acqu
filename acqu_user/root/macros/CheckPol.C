#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TChainElement.h"
#include "TChain.h"
#include <iostream>
#include <fstream>

using namespace std;

void CheckPol( TString inputFile = "/w/work14/simong/EdamOutSimple/PARA193*.root", TString inputFile3 = "/w/work14/simong/EdamOutSimple/AM*.root", TString outFile = "/scratch/simong/FileCheckPol.root" ){

  //Open input files

  //  cout << inputFile  << endl;
  TChain* oTree   = new TChain("Recon","datachain");
  TChain* oTree2  = new TChain("Recon","datachain");
  oTree->Add(inputFile);
  oTree2->Add(inputFile3);
  //TFile*  iFile   = new TFile(inputFile,"read");
  TFile*  oFile   = new TFile(outFile,"recreate");
  //if(!(oTree = (TTree*)iFile->Get("Recon"))) return;
  TTree*  outTree = new TTree("React", "Reaction channel filtered final data");

  const Int_t fMaxParticles = 6;
  const Int_t fMaxTagged    = 100;

  //Input Tree
  Int_t          fNParticleTree;
  TClonesArray   *fParticles;
  Int_t          fPDG[fMaxParticles];
  Int_t          fRecoil;
  Double_t       fTargetM;
  TLorentzVector *fRecoilP4;
  TClonesArray   *fMMP4;
  TClonesArray   *fRealMMP4;
  Double_t       fEGamma[fMaxTagged];
  Double_t       fTGamma[fMaxTagged];
  Int_t          fTaggedN;
  Int_t          fPrompt[fMaxTagged];
  Int_t          fHeli[fMaxTagged];
  Int_t          fNPhotons;
  Int_t          fNUnknowns;
  Double_t       fEventQual;
  Double_t       fDeltaE[fMaxTagged];
  Double_t       fAngleDiff[fMaxTagged];
  
  //Calculation vectors
  TLorentzVector *fInput     = new TLorentzVector();
  TLorentzVector *fTargetP4  = new TLorentzVector();
  TLorentzVector *fTaggedP4  = new TLorentzVector();
  TLorentzVector *fTempPart  = new TLorentzVector();
  TLorentzVector *fTempPart2 = new TLorentzVector();
  TVector3       *fCMBoost   = new TVector3();
  TVector3       fTempVec1;
  TVector3       fTempVec2;
  TVector3       fMMVec   ;
  TVector3       fTagVec  ;
  
  //For Recoilangle Calculation
  TRotation      RotateFrame;
  TVector3       XAxis;
  TVector3       YAxis;
  TVector3       ZAxis; 
  TVector3       fScattered;
  TVector3       fReconstructed;
  TVector3       fDetected;

  Double_t fTempMM;
  Double_t fTempPhi;
  
  UInt_t bins      = 352;
  Double_t *fEBins = new Double_t[bins+1];
  //Beam Energy Bins
  ifstream f;
  Double_t tempE;
  f.open("/home/simong/A2Edam/data/EgBins.txt.FULL");
  for(UInt_t i = 0; i <= bins; i++){
    
    f >> tempE;
    fEBins[i] = tempE;
    //    fEBins[i] = ToCME(tempE);

  }
  f.close();

  fRecoilP4  = new TLorentzVector();
  fParticles = new TClonesArray("TLorentzVector",fMaxParticles);
  fMMP4      = new TClonesArray("TLorentzVector",fMaxTagged);
  fRealMMP4  = new TClonesArray("TLorentzVector",fMaxTagged);

  oTree->SetBranchAddress("FNParticle",&fNParticleTree);
  oTree->SetBranchAddress("Particles",&fParticles);
  oTree->SetBranchAddress("PDG",fPDG);
  oTree->SetBranchAddress("Recoil",&fRecoil);
  oTree->SetBranchAddress("RecoilP4",fRecoilP4);
  oTree->SetBranchAddress("TargetM",&fTargetM);
  oTree->SetBranchAddress("MMP4",&fMMP4);
  oTree->SetBranchAddress("RealMMP4",&fRealMMP4);
  oTree->SetBranchAddress("EGamma",fEGamma);
  oTree->SetBranchAddress("TGamma",fTGamma);
  oTree->SetBranchAddress("TaggedN",&fTaggedN);
  oTree->SetBranchAddress("Prompt",fPrompt);
  oTree->SetBranchAddress("Heli",fHeli);
  oTree->SetBranchAddress("NPhotons",&fNPhotons);
  oTree->SetBranchAddress("NUnknowns",&fNUnknowns);
  oTree->SetBranchAddress("EventQual",&fEventQual);
  oTree->SetBranchAddress("DeltaE",fDeltaE);
  oTree->SetBranchAddress("AngleDiff",fAngleDiff);

  oTree->SetBranchStatus("*",0);
  oTree->SetBranchStatus("FNParticle",1);
//   oTree->SetBranchStatus("Particles",1);
  oTree->SetBranchStatus("MMP4",1);
  oTree->SetBranchStatus("TaggedN",1);
  oTree->SetBranchStatus("EGamma",1);
  oTree->SetBranchStatus("TGamma",1);
  oTree->SetBranchStatus("Prompt",1);
//   oTree->SetBranchStatus("Heli",1);
  //   oTree->SetBranchStatus("NPhotons",1);
//   oTree->SetBranchStatus("PDG",1);
//   oTree->SetBranchStatus("TargetM",1);
//   oTree->SetBranchStatus("Recoil*",1);

  oTree2->SetBranchAddress("FNParticle",&fNParticleTree);
  oTree2->SetBranchAddress("Particles",&fParticles);
  oTree2->SetBranchAddress("PDG",fPDG);
  oTree2->SetBranchAddress("Recoil",&fRecoil);
  oTree2->SetBranchAddress("RecoilP4",fRecoilP4);
  oTree2->SetBranchAddress("TargetM",&fTargetM);
  oTree2->SetBranchAddress("MMP4",&fMMP4);
  oTree2->SetBranchAddress("RealMMP4",&fRealMMP4);
  oTree2->SetBranchAddress("EGamma",fEGamma);
  oTree2->SetBranchAddress("TGamma",fTGamma);
  oTree2->SetBranchAddress("TaggedN",&fTaggedN);
  oTree2->SetBranchAddress("Prompt",fPrompt);
  oTree2->SetBranchAddress("Heli",fHeli);
  oTree2->SetBranchAddress("NPhotons",&fNPhotons);
  oTree2->SetBranchAddress("NUnknowns",&fNUnknowns);
  oTree2->SetBranchAddress("EventQual",&fEventQual);
  oTree2->SetBranchAddress("DeltaE",fDeltaE);
  oTree2->SetBranchAddress("AngleDiff",fAngleDiff);

  oTree2->SetBranchStatus("*",0);
  oTree2->SetBranchStatus("FNParticle",1);
//   oTree2->SetBranchStatus("Particles",1);
  oTree2->SetBranchStatus("MMP4",1);
  oTree2->SetBranchStatus("TaggedN",1);
  oTree2->SetBranchStatus("EGamma",1);
  oTree2->SetBranchStatus("TGamma",1);
  oTree2->SetBranchStatus("Prompt",1);
//   oTree2->SetBranchStatus("Heli",1);
  //   oTree2->SetBranchStatus("NPhotons",1);
//   oTree2->SetBranchStatus("PDG",1);
//   oTree2->SetBranchStatus("TargetM",1);
//   oTree2->SetBranchStatus("Recoil*",1);



  Int_t    NPhotoCut  = 2;
  const  UInt_t    NFiles     = oTree->GetNtrees();
  
  Double_t *FileBins = new Double_t[NFiles+1];
  for(UInt_t i = 0; i <= NFiles; i++){
    
    FileBins[i] = i;

  }

  TObjArray *fileElements=oTree->GetListOfFiles();
  TIter next(fileElements);
  TChainElement *chEl=0;
  while (( chEl=(TChainElement*)next() )) {
    
    TString temp = chEl->GetTitle();
    size_t begin = temp.find_last_of("/");
    size_t end = temp.find_last_of(".");
    TString temp2 = temp.substr(begin, end - begin + 1);
    cout << temp << endl;
    cout << temp2 << endl;
  }

  TH2F*    hFileEG        = new TH2F("FileEG",
				     "FileEG",
				     bins, fEBins,
				     NFiles, FileBins);

  TH1F*    hFileEGamor    = new TH1F("FileEGamor",
				     "FileEGamor",
				     bins, fEBins);

  TH2F*    hFileEGenh     = new TH2F("FileEGenh",
				     "FileEGenh",
				     bins, fEBins,
				     NFiles, FileBins);


  Long64_t ientry = 0;
  Int_t outEntries = 0;

  cout << oTree->GetNtrees() << endl;
  cout << oTree->GetEntries() << endl;

//   TObjArray *fileElements = oTree->GetListOfFiles();
  
//   fileElements

//   for (Int_t i=0;i<NFiles;i++){

//     TIter next(fileElements);

//     cout << oTree->GetListOfFiles()[i] << endl;

//   }
  
//    for(Int_t i = 0; i<200000; i+=1){
// //  for(Int_t i = 0; i<oTree->GetEntries(); i+=1){
    
//     ientry = oTree->LoadTree(i);

//     if(ientry>200000) continue;

//     oTree->GetEntry(i);

//     if(i%5000 == 0){
//       cout << '\r' << "Reading File: " << oTree->GetTreeNumber() <<  " Tree Entry: " << ientry << " Events Processed: " << i << " Output Entries: " << outEntries << flush;
//     }

    
//     for(Int_t j=0;j<fTaggedN;j++){

//       if(fTGamma[j]<40 || fTGamma[j]>160) continue;
//       if(fPrompt[j]==1) hFileEG->Fill(fEGamma[j],oTree->GetTreeNumber(),1);
//       if(fPrompt[j]==2) hFileEG->Fill(fEGamma[j],oTree->GetTreeNumber(),-0.125); 
      
//     }
//     outEntries++;  
//   }

//   cout << oTree2->GetNtrees() << endl;
//   cout << oTree2->GetEntries() << endl;
  
//   for(UInt_t i = 0; i<oTree2->GetEntries(); i+=1){
    
//     ientry = oTree2->LoadTree(i);

//     if(ientry>200000) continue;

//     oTree2->GetEntry(i);

//     if(i%5000 == 0){
//       cout << '\r' << "Reading File: " << oTree2->GetTreeNumber() <<  " Tree Entry: " << ientry << " Events Processed: " << i << " Output Entries: " << outEntries << flush;
//     }
    

//     for(Int_t j=0;j<fTaggedN;j++){

//       if(fTGamma[j]<40 || fTGamma[j]>160) continue;
      
//       if(fPrompt[j]==1) hFileEGamor->Fill(fEGamma[j],1);
//       if(fPrompt[j]==2) hFileEGamor->Fill(fEGamma[j],-0.1);
      
//     }
//     outEntries++;  
//   }

  Double_t ratio   = 0;
  Int_t    normBin = 150;

  //Get enhancement
  for(UInt_t i=1;i<NFiles;i++){

    if(hFileEG->GetBinContent(normBin,i)==0) continue;

    ratio = hFileEGamor->GetBinContent(normBin)/hFileEG->GetBinContent(normBin,i);

    for(UInt_t j=1;j<bins;j++){

      if(hFileEG->GetBinContent(j,i)==0) continue;
      if(hFileEGamor->GetBinContent(j)==0) continue;
      hFileEGenh->SetBinContent(j,i,(ratio*hFileEG->GetBinContent(j,i))/hFileEGamor->GetBinContent(j));

    }
    
  }
  

  hFileEG->Write();
  hFileEGamor->Write();
  hFileEGenh->Write();

  oFile->Close();
  //iFile->Close();

  cout << endl;
  cout << inputFile << endl;
  cout << outFile << endl << endl;

}
