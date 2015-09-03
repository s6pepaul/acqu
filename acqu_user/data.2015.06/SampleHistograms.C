double GausOnBase(double *x, double *par) {
  double arg = 0;
  if (par[2] != 0) arg = (x[0] - par[1])/par[2];
  double fitval = par[3] + par[0]*TMath::Exp(-0.5*arg*arg);
  return fitval;
}

void SampleHistograms(){
  TFile *file = new TFile("sampleHisto.root", "READ");
  
  TH1F* HInc = (TH1F*) file->Get("LinPol_Incoherent;1");
  TH1F* HGatedInc = (TH1F*) file->Get("LinPol_GatedIncoherent;1");
  TH1F* HCoh = (TH1F*) file->Get("LinPol_Coherent;1");
  TH1F* HGatedCoh = (TH1F*) file->Get("LinPol_GatedCoherent;1");
  TH1F* HEnh = HInc->Clone("HEnh"); HEnh->Reset("ICES");
  TH1F* HGatedEnh = HInc->Clone("HGatedEnh"); HGatedEnh->Reset("ICES");
  TH1F* HSampleInc = HInc->Clone("HSampleInc"); HSampleInc->Reset("ICES");
  TH1F* HSampleEnh = HInc->Clone("HSampleEnh"); HSampleEnh->Reset("ICES");
  TH1F* HSampleEdge = new TH1F("HSampleEdge", "SampleEdge", 4000, 0, 4000);
  TH1F* HSampleGatedInc = HInc->Clone("HSampleGatedInc"); HSampleGatedInc->Reset("ICES");
  TH1F* HSampleGatedEnh = HInc->Clone("HSampleGatedEnh"); HSampleGatedInc->Reset("ICES");
  TH1F* HSampleGatedEdge = new TH1F("HSampleGatedEdge", "SampleGatedEdge", 4000, 0, 4000);
  
  double EdgeMin = 735.0;
  double EdgeMax = 800.0;
  double normenergy = 900.0;
  double normbin = 200; //default

  TCanvas* CInc = new TCanvas("CInc", "Incoherent", 1);
  HInc->SetLineColor(4);
  HInc->SetStats(false);
  HInc->Draw();
  CInc->Update();

  TCanvas* CGatedInc = new TCanvas("CGatedInc", "GatedIncoherent", 1);
  HGatedInc->SetLineColor(2);
  HGatedInc->SetStats(false);
  HGatedInc->Draw();
  CGatedInc->Update();

  TCanvas* CCoh = new TCanvas("CCoh", "Coherent", 1);
  HCoh->SetLineColor(4);
  HCoh->SetStats(false);
  HCoh->Draw();
  CCoh->Update();

  TCanvas* CGatedCoh = new TCanvas("CGatedCoh", "GatedCoherent", 1);
  HGatedCoh->SetLineColor(2);
  HGatedCoh->SetStats(false);
  HGatedCoh->Draw();
  CGatedCoh->Update();

  TCanvas* CEnh = new TCanvas("CEnh", "Enhancement", 1);
  TLegend* LEnh = new TLegend(0.7,0.75,0.88,0.88);
  HEnh->Divide(HCoh,HInc,100);
  //HEnh->Sumw2();
  HGatedEnh->Divide(HGatedCoh,HGatedInc,100);
  //HGatedEnh->Sumw2();
  if(HEnh->GetBinContent(1)<1.0){
    HEnh->SetBinContent(1,0);
  }
  if(HGatedEnh->GetBinContent(1)<1.0){
    HGatedEnh->SetBinContent(1,0);
  }
  for(int i=2;i<=352;i++){
    if(HEnh->GetBinContent(i)<1.0){
      HEnh->SetBinContent(i,(HEnh->GetBinContent(i-1)));
    }
    if(HGatedEnh->GetBinContent(i)<1.0){
      HGatedEnh->SetBinContent(i,(HGatedEnh->GetBinContent(i-1)));
    }
  }
  if(HEnh!=NULL){
    HEnh->SetStats(false);
    HEnh->SetMinimum(0.0);
    HEnh->SetMaximum(250);
    HEnh->SetLineColor(4);
    HEnh->Draw();
  }
  if(HGatedEnh!=NULL){
    HGatedEnh->SetStats(false);
    HGatedEnh->SetLineColor(2);
    HGatedEnh->Draw("SAME");
  }

  ////////////////////////////////////////////////////////////    find edge
  int binx=0;
  double ymax=0;
  for(int b=HEnh->FindBin(EdgeMin);b<HEnh->FindBin(EdgeMax);b++){
    if(HEnh->GetBinContent(b)>ymax){
      ymax=HEnh->GetBinContent(b);
      binx=b;
    }
  }
  double xmax=HEnh->GetBinCenter(binx);
  ymax=HEnh->GetBinContent(binx);
  
  TF1* edgeFit = new TF1("edgeFit",GausOnBase,0,100,4);
  edgeFit->SetRange(xmax,xmax+40.0);
  edgeFit->SetParameter(1,xmax);
  edgeFit->SetParameter(2,10.0);
  edgeFit->SetParameter(3,100.0);
  edgeFit->SetLineColor(4);
  edgeFit->SetLineWidth(2);
  HEnh->Fit(edgeFit,"QNR");
  edgeFit->Print();
  
  if(edgeFit!=NULL){
    edgeFit->Draw("SAME");
  }

  /////////////////////////////////////////////////////////    find gated edge
  binx=0;
  ymax=0;
  for(int b=HGatedEnh->FindBin(EdgeMin);b<HGatedEnh->FindBin(EdgeMax);b++){
    if(HGatedEnh->GetBinContent(b)>ymax){
      ymax=HGatedEnh->GetBinContent(b);
      binx=b;
    }
  }
  xmax=HGatedEnh->GetBinCenter(binx);
  ymax=HGatedEnh->GetBinContent(binx);
  
  TF1* edgeFitGated = new TF1("edgeFitGated",GausOnBase,0,100,4);
  edgeFitGated->SetRange(xmax,xmax+40.0);
  edgeFitGated->SetParameter(1,xmax);
  edgeFitGated->SetParameter(2,10.0);
  edgeFitGated->SetParameter(3,100.0);
  edgeFitGated->SetLineColor(2);
  edgeFitGated->SetLineWidth(2);
  HGatedEnh->Fit(edgeFitGated,"QNR");
  edgeFitGated->Print();

  if(edgeFitGated!=NULL){
    edgeFitGated->Draw("SAME");
  }
  LEnh->AddEntry(HEnh,"free running scalers", "l");
  LEnh->AddEntry(HGatedEnh, "gated scalers", "l");
  LEnh->Draw();
  CEnh->Update();

  ///////////////////////////////////////////////////////////////////  sample from histogram

  double edgepos;
  int samplesize = 50000; //number of events in one buffer
  int buffer = 10;
  int ncount = 0;
  double norm = 0;
  int ngcount = 0;
  double gnorm = 0;
  double incdummy, cohdummy;
  int filecount=1;

  TCanvas* CSample = new TCanvas("CSample", "Sample", 1);
  TLegend* LSample = new TLegend(0.7,0.75,0.88,0.88);
  TF1* edgeFitGatedSample = new TF1("edgeFitGatedSample",GausOnBase,0,100,4);
  TF1* edgeFitSample = new TF1("edgeFitSample",GausOnBase,0,100,4);

  for(int i=1; i<800;i++){
    HSampleInc->Reset("ICES");
    HSampleInc->FillRandom(HInc, buffer*samplesize);
    
    ncount = 0;               // calculate normalisation
    for(int n=0;n<352;n++){
      if(HSampleInc->GetBinCenter(n)>=normenergy){
	normbin=n;
	break;
      }
    }
    for(int v = normbin; v<(normbin+20); v++){
      incdummy = HSampleInc->GetBinContent(v);
      cohdummy = HCoh->GetBinContent(v);
      if((incdummy>1.0) && (cohdummy>1.0)){
	norm+=100.0*(incdummy/cohdummy);
	ncount++;
      }
    }
    norm/=ncount;
    
    HSampleEnh->Divide(HCoh,HSampleInc,norm);
    
    if(HSampleEnh->GetBinContent(1)<1.0){
      HSampleEnh->SetBinContent(1,0);
    }
    for(int n=2;n<=352;n++){
      if(HSampleEnh->GetBinContent(n)<1.0){
	HSampleEnh->SetBinContent(n,(HSampleEnh->GetBinContent(n-1)));
      }
    }
    for(int n=100;n<=352;n++){
      if(HSampleEnh->GetBinContent(n)>1.3*HSampleEnh->GetBinContent(n-1)){
	HSampleEnh->SetBinContent(n,(HSampleEnh->GetBinContent(n-1)));
      }
    }
    HSampleEnh->SetLineColor(4);
    HSampleEnh->SetMinimum(0);
    HSampleEnh->SetMaximum(300);
    HSampleEnh->SetStats(false);
    HSampleEnh->Draw();
    
    binx=0;// find edge
    ymax=0;
    for(int b=HSampleEnh->FindBin(EdgeMin);b<HSampleEnh->FindBin(EdgeMax);b++){
      if(HSampleEnh->GetBinContent(b)>ymax){
	ymax=HSampleEnh->GetBinContent(b);
	binx=b;
      }
    }
    xmax=HSampleEnh->GetBinCenter(binx);
    ymax=HSampleEnh->GetBinContent(binx);
    
    edgeFitSample->SetRange(xmax-5.0,xmax+50.0);
    edgeFitSample->SetParameter(1,xmax);
    edgeFitSample->SetParLimits(1,xmax-20,xmax+20);
    edgeFitSample->SetParameter(2,10.0);
    edgeFitSample->SetParLimits(2,0,50);
    edgeFitSample->SetParameter(3,100.0);
    edgeFitSample->SetLineColor(2);
    edgeFitSample->SetLineWidth(2);
    HSampleEnh->Fit(edgeFitSample,"QNR");
    edgeFitSample->Print();
    edgepos = edgeFitSample->GetParameter(1) + abs(edgeFitSample->GetParameter(2));

    
    if(edgeFitSample!=NULL){
      edgeFitSample->SetLineColor(4);
      edgeFitSample->Draw("SAME");
    }

    HSampleEdge->Fill(i,edgepos);
    
    
    HSampleGatedInc->Reset("ICES");
    HSampleGatedInc->FillRandom(HGatedInc, buffer*samplesize);
    ngcount = 0;               // calculate normalisation
    for(int n=0;n<352;n++){
      if(HSampleGatedInc->GetBinCenter(n)>=normenergy){
	normbin=n;
	break;
      }
    }
    for(int v = normbin; v<(normbin+20); v++){
      incdummy = HSampleGatedInc->GetBinContent(v);
      cohdummy = HGatedCoh->GetBinContent(v);
      if((incdummy>1.0) && (cohdummy>1.0)){
	gnorm+=100.0*(incdummy/cohdummy);
	ngcount++;
      }
    }
    gnorm/=ngcount;

    HSampleGatedEnh->Divide(HGatedCoh,HSampleGatedInc,gnorm);
    
    if(HSampleGatedEnh->GetBinContent(1)<1.0){
      HSampleGatedEnh->SetBinContent(1,0);
    }
    for(int n=2;n<=352;n++){
      if(HSampleGatedEnh->GetBinContent(n)<1.0){
	HSampleGatedEnh->SetBinContent(n,(HSampleGatedEnh->GetBinContent(n-1)));
      }
    }
    for(int n=100;n<=352;n++){
      if(HSampleGatedEnh->GetBinContent(n)>1.3*HSampleGatedEnh->GetBinContent(n-1)){
	HSampleGatedEnh->SetBinContent(n,(HSampleGatedEnh->GetBinContent(n-1)));
      }
    }
    HSampleGatedEnh->SetLineColor(2);
    HSampleGatedEnh->SetStats(false);
    HSampleGatedEnh->Draw("SAME");
    
    binx=0;// find edge
    ymax=0;
    for(int b=HSampleGatedEnh->FindBin(EdgeMin);b<HSampleGatedEnh->FindBin(EdgeMax);b++){
      if(HSampleGatedEnh->GetBinContent(b)>ymax){
	ymax=HSampleGatedEnh->GetBinContent(b);
	binx=b;
      }
    }
    xmax=HSampleGatedEnh->GetBinCenter(binx);
    ymax=HSampleGatedEnh->GetBinContent(binx);
    
    edgeFitGatedSample->SetRange(xmax-5.0,xmax+50.0);
    edgeFitGatedSample->SetParameter(1,xmax);
    edgeFitGatedSample->SetParLimits(1,xmax-20,xmax+20);
    edgeFitGatedSample->SetParameter(2,10.0);
    edgeFitGatedSample->SetParLimits(2,0,50);
    edgeFitGatedSample->SetParameter(3,100.0);
    edgeFitGatedSample->SetLineColor(2);
    edgeFitGatedSample->SetLineWidth(2);
    HSampleGatedEnh->Fit(edgeFitGatedSample,"QNR");
    edgeFitGatedSample->Print();
    edgepos = edgeFitGatedSample->GetParameter(1) + abs(edgeFitGatedSample->GetParameter(2));

    
    if(edgeFitGatedSample!=NULL){
      edgeFitGatedSample->SetLineColor(2);
      edgeFitGatedSample->Draw("SAME");
    }
    if(LSample->GetNRows()<2){
      LSample->AddEntry(HSampleEnh, "free running scalers", "l");
      LSample->AddEntry(HSampleGatedEnh, "gated scalers", "l");
    }
    LSample->Draw();

    if((i%100)==0){
      std::ostringstream filenameStream;
      filenameStream << "sample_" << buffer << "_" << filecount << ".pdf";
      std::string filename = filenameStream.str();
      CSample->SaveAs(filename.c_str());
      filecount++;
    }
    
    HSampleGatedEdge->Fill(i,edgepos);
  }
  HSampleEdge->SetLineColor(4);
  HSampleGatedEdge->SetLineColor(2);
  HSampleEdge->SetStats(false);
  HSampleGatedEdge->SetStats(false);
  HSampleGatedEdge->SetMinimum(740);
  HSampleGatedEdge->SetMaximum(760);
  HSampleGatedEdge->Draw();
  HSampleEdge->Draw("SAME");
  LSample->Draw();
  CSample->Update();
  
  // file->Close();
}
