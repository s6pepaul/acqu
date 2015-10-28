// SVN Info: $Id$

/*************************************************************************
 * Author: Irakli Keshelashvili, Dominik Werthmueller
   Modified: Farah Afzal, Karsten Spieker in October 2015
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCCalibEnergy                                                        //
//                                                                      //
// Base energy calibration module class.                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TCCalibEnergy.h"
#include "TLatex.h"

    Double_t lower_limit;		//rejection range while fitting background
    Double_t higher_limit;
//______________________________________________________________________________
TCCalibEnergy::TCCalibEnergy(const Char_t* name, const Char_t* title, const Char_t* data,
                             Int_t nElem)
    : TCCalib(name, title, data, nElem)
{
    // Constructor.
    fTAPSType=kTAPS_2009;
    // init members
    fPi0Pos = 0;
    fLine = 0;
    fDetectorView = 0;
}

//______________________________________________________________________________
TCCalibEnergy::~TCCalibEnergy()
{
    // Destructor. 
    
    if (fLine) delete fLine;
}

//______________________________________________________________________________
void TCCalibEnergy::Init()
{
    // Init the module.
    
    Char_t tmp[256];

    // init members

    fPi0Pos = 0;
    fLine = new TIndicatorLine();
    
    // get histogram name
    sprintf(tmp, "%s.Histo.Fit.Name", GetName());
    if (!TCReadConfig::GetReader()->GetConfig(tmp))
    {
        Error("Init", "Histogram name was not found in configuration!");
        return;
    }
    else fHistoName = *TCReadConfig::GetReader()->GetConfig(tmp);
    
    // read old parameters (only from first set)
    TCMySQLManager::GetManager()->ReadParameters(fData, fCalibration.Data(), fSet[0], fOldVal, fNelem);
    
    // copy to new parameters
    for (Int_t i = 0; i < fNelem; i++) fNewVal[i] = fOldVal[i];

    // sum up all files contained in this runset
    TCFileManager f(fData, fCalibration.Data(), fNset, fSet);
    
    // get the main calibration histogram
    fMainHisto = f.GetHistogram(fHistoName.Data());
    if (!fMainHisto)
    {
        Error("Init", "Main histogram does not exist!\n");
        return;
    }
    
    // create the overview histogram
    fOverviewHisto = new TH1F("Overview", ";Element;2#gamma inv. mass [MeV]", fNelem, 0, fNelem);
    fOverviewHisto->SetMarkerStyle(2);
    fOverviewHisto->SetMarkerColor(4);
    
    // draw main histogram
    fCanvasFit->Divide(1, 2, 0.001, 0.001);
    fCanvasFit->cd(1)->SetLogz();
    sprintf(tmp, "%s.Histo.Fit", GetName());
    TCUtils::FormatHistogram(fMainHisto, tmp);
    fMainHisto->Draw("colz");

    // draw the overview histogram
    fCanvasResult->cd();
    sprintf(tmp, "%s.Histo.Overview", GetName());
    TCUtils::FormatHistogram(fOverviewHisto, tmp);
    fOverviewHisto->Draw("P");
}


Double_t fit_poldeg3(Double_t *x, Double_t *par){
	Double_t xx = x[0];

  	return  par[0] + par[1]*xx + par[2]*TMath::Power(xx,2) + par[3]*TMath::Power(xx,3);
}

Double_t fit_poldeg3_bg(Double_t *x, Double_t *par){
	Double_t xx = x[0];

	if ((xx > lower_limit) && (xx < higher_limit)) {
         TF1::RejectPoint();
   	}

  	return  par[0] + par[1]*xx + par[2]*TMath::Power(xx,2) + par[3]*TMath::Power(xx,3);
}

Double_t fitgaus(Double_t *x, Double_t *par)
{
   Double_t arg = 0;
   if (par[2] != 0) arg = (x[0] - par[1])/par[2];

   Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
   return fitval;
}

//Novosibirsk function for fitting the pi0 peak
Double_t fit_Novo(Double_t *x, Double_t *par){

  Double_t Amp   = par[0];
  Double_t peak  = par[1];
  Double_t width = par[2];
  Double_t tail  = par[3];
  Double_t qa=0,qb=0,qc=0,qx=0,qy=0;


  	if(0 != width){
     if(TMath::Abs(par[3]) < 1.e-7){
       	qc = 0.5*TMath::Power(((x[0]-peak)/width),2);
       }
     else {
   	qa = tail*sqrt(log(4.));
   	qb = sinh(qa)/qa;
   	qx = (x[0]-peak)/width*qb;
   	qy = 1.+tail*qx;

       	if( qy > 1.E-7)
 	     	qc = 0.5*(TMath::Power((log(qy)/tail),2) + tail*tail);
       	else
            qc = 15.0;
     }
   return Amp*exp(-qc);
	}
	else return 0;

}

//fit function for the entire spectrum (background+fFitPeaknal)
Double_t myFitFct(Double_t *x, Double_t *par){
  return fit_poldeg3(x , par) + fit_Novo(x, &par[4]); 
}
        
TF1* bg1; 

//______________________________________________________________________________
void TCCalibEnergy::Fit(Int_t elem)
{
    // Perform the fit of the element 'elem'.
    
    Char_t tmp[256];
    
    // create histogram projection for this element
    sprintf(tmp, "ProjHisto_%i", elem);
    TH2* h2 = (TH2*) fMainHisto;
    if (fFitHisto) delete fFitHisto;
    fFitHisto = (TH1D*) h2->ProjectionX(tmp, elem+1, elem+1, "e");
    
    // draw histogram
    fFitHisto->SetFillColor(35);
    fCanvasFit->cd(2);
    sprintf(tmp, "%s.Histo.Fit", GetName());
    TCUtils::FormatHistogram(fFitHisto, tmp);
    fFitHisto->GetYaxis()->SetRangeUser(0.0, fFitHisto->GetMaximum()*1.05);
    fFitHisto->Draw("hist");
     
    // check for sufficient statistics
    if (fFitHisto->GetEntries() > 500)
    {
        // delete old function
        if (fFitPeak) delete fFitPeak;
        if (bg1) delete bg1;
        if (fFitBackGround) delete fFitBackGround;
        if (fFitFunc) delete fFitFunc;

        sprintf(tmp, "fEnergy_%i", elem);
        
        // estimate peak position
        fPi0Pos = fFitHisto->GetBinCenter(fFitHisto->GetMaximumBin());
        if (fPi0Pos < 100 || fPi0Pos > 160) fPi0Pos = 135;

        // configure fitting function
        initFitFunction(elem);

	fFitHisto->Fit(fFitFunc, "RBQ0");

	for (Int_t i=0; i<4; ++i) fFitBackGround->FixParameter(i,fFitFunc->GetParameter(i));
	for (Int_t i=0; i<4; ++i) fFitPeak->FixParameter(i,fFitFunc->GetParameter(i+4));

	// final results
       	fPi0Pos = fFitFunc->GetParameter(5); 
        
        // check if mass is in normal range
        if (fPi0Pos < 80 || fPi0Pos > 200) fPi0Pos = 135;
 
        // set indicator line
        fLine->SetY1(0);
        fLine->SetupY(0, fFitHisto->GetMaximum()+20);
        fLine->SetX1(fPi0Pos);
        fLine->SetX2(fPi0Pos);
   
        // draw fitting function
        if (fFitFunc) {

          fFitPeak->Draw("same");   
          fFitPeak->SetLineColor(2);
          fFitBackGround->Draw("same"); 
          fFitBackGround->SetLineColor(kOrange);
          fFitFunc->Draw("same");
          fFitFunc->SetLineColor(kGreen+1);

        }
    
        // draw indicator line
        fLine->Draw();

        double xmin,xmax;
        fFitFunc->GetRange(xmin,xmax);
        double maximum = fFitFunc->GetMaximumX();
        fFitOk=false;
        if(maximum < xmax && maximum > xmin) {
            puts("Fit OK!\n");
            fFitOk=true;
        } 
        else {
            puts("No Maximum found in fit!");
        }

    }



    if(fDetectorView) {
        fDetectorView->SetElement(elem, vCurrPos);
        fExtraCanvas->cd();
        fDetectorView->Draw("col");
        fExtraCanvas->Update();
    }

    // update canvas
    fCanvasFit->Update();

    // update overview
    if (elem % 20 == 0)
    {
        fCanvasResult->cd();
        fOverviewHisto->Draw("E1");
        fCanvasResult->Update();

    }   
}

//______________________________________________________________________________
void TCCalibEnergy::Calculate(Int_t elem)
{
    // Calculate the new value of the element 'elem'.
    
    Bool_t unchanged = kFALSE;

    // check if fit was performed
    if (fFitHisto->GetEntries() > 500)
    {
        // check if line position was modified by hand
        if (fLine->GetX1() != fPi0Pos) fPi0Pos = fLine->GetX1();
        
        // calculate the new gain
        // apply fConvergenceFactor only to the desired procentual change of fOldVal,
        // given by (TCConfig::kPi0Mass / fPi0Pos - 1)
        fNewVal[elem] = fOldVal[elem] + fOldVal[elem] * fConvergenceFactor * (TCConfig::kPi0Mass / fPi0Pos - 1);
        
        // if new value is negative take old
        if (fNewVal[elem] < 0) 
        {
            fNewVal[elem] = fOldVal[elem];
            unchanged = kTRUE;
        }

        // update overview histogram
        fOverviewHisto->SetBinContent(elem+1, fPi0Pos);
        fOverviewHisto->SetBinError(elem+1, 0.0000001);
    
        // update average calculation
        fAvr += fPi0Pos;
        fAvrDiff += TMath::Abs(fPi0Pos - TCConfig::kPi0Mass);
        fNcalc++;
    }
    else
    {   
        // do not change old value
        fNewVal[elem] = fOldVal[elem];
        unchanged = kTRUE;
    }

    // user information
    printf("Element: %03d    Pi0: %12.8f    "
           "old gain: %12.8f    new gain: %12.8f    diff: %6.2f %%",
           elem, fPi0Pos, fOldVal[elem], fNewVal[elem],
           TCUtils::GetDiffPercent(fOldVal[elem], fNewVal[elem]));
    if (unchanged) printf("    -> unchanged");
    if (this->InheritsFrom("TCCalibCBEnergy"))
    {
        if (TCUtils::IsCBHole(elem)) printf(" (hole)");
    }
    printf("\n");

    // show average
    if (elem == fNelem-1)
    {
        fAvr /= (Double_t)fNcalc;
        fAvrDiff /= (Double_t)fNcalc;
        printf("Average pi0 position           : %.3f MeV\n", fAvr);
        printf("Average difference to pi0 mass : %.3f MeV\n", fAvrDiff);
    }

    if(fDetectorView) {
        fDetectorView->SetElement(elem, fFitOk ? vFitOK : vFitFailed);
        fExtraCanvas->cd();
        fDetectorView->Draw("col");
        fExtraCanvas->Update();
    }
}   


void TCCalibCBEnergy::initFitFunction(Int_t elem)
{
   //limit for pi0 peak only, important for bg1 which excludes the area of peak from fit
   lower_limit = 100;
   higher_limit = 160;

   fFitPeak    = new TF1("peak",fit_Novo,55, 220,4);
   bg1    = new TF1("bg1",fit_poldeg3_bg, 55, 220,4); 
   fFitBackGround    = new TF1("fFitBackGround",fit_poldeg3, 55,220,4); 
   fFitFunc = new TF1("fFitFunc", myFitFct, 55,220,8);

   for (Int_t i=0; i<4; ++i) bg1->SetParameter(i,1.); 
   fFitHisto->Fit(bg1, "RBQ0");
   
   for (Int_t i=0; i<4; ++i) fFitFunc->FixParameter(i,bg1->GetParameter(i));

   fFitFunc->SetParLimits(4,1, 20000000);            // amplitude
   fFitFunc->SetParLimits(5,115, 165);   	//peak position
   fFitFunc ->SetParLimits(6, 5,20);   	//width
   fFitFunc ->SetParLimits(7, -0.3,0.3);   	//tail factor
}

void TCCalibTAPSEnergyLG::initFitFunction(Int_t elem)
{
   Int_t ring_nr = TCUtils::GetTAPSRing(elem, fTAPSType); //the fit range and rejection range for the peak is set separately for each ring see below

   Int_t limit_low[11]={118,118,112, 102, 102, 100, 100, 98, 98,98, 98};//range for the pi0 peak which is rejected
   Int_t limit_high[11]={155,155,155, 155, 155, 155, 160, 155, 155,158, 158};

   Int_t range_low[11]={90,70,60, 70, 50, 50, 40, 35, 35,40, 30}; //range for the fit range
   Int_t range_high[11]={250,250,250, 250, 250, 250, 250, 250, 250,250, 250};

   lower_limit = limit_low[ring_nr-1];//limit for peak only, important for bg1 which excludes the area of peak from fit
   higher_limit = limit_high[ring_nr-1];

   fFitPeak    = new TF1("peak",fit_Novo,range_low[ring_nr-1],range_high[ring_nr-1],4);
   bg1    = new TF1("bg1",fit_poldeg3_bg, range_low[ring_nr-1],range_high[ring_nr-1],4); 
   fFitBackGround    = new TF1("fFitBackGround",fit_poldeg3, range_low[ring_nr-1],range_high[ring_nr-1],4); 
   fFitFunc = new TF1("fFitFunc", myFitFct, range_low[ring_nr-1],range_high[ring_nr-1],8);

   for (Int_t i=0; i<4; ++i) bg1->SetParameter(i,1/pow(10.,4*i));
   fFitHisto->Fit(bg1, "RBQ0");
   
   for (Int_t i=0; i<4; ++i) fFitFunc->FixParameter(i,bg1->GetParameter(i));

   fFitFunc->SetParLimits(4,1, 200000); 
   fFitFunc->SetParLimits(5,115, 140);
   fFitFunc ->SetParLimits(6, 5,20);
   fFitFunc ->SetParLimits(7, -0.3,0.3);

   TLatex *m1 = new TLatex();
   m1->SetTextColor(kRed);
   m1->SetTextSize(0.07);
   m1->SetTextFont(62);
   m1->SetNDC();
   m1->DrawLatex(0.7,0.835,Form("Ring %i",ring_nr));
   m1->Draw("same");
}

ClassImp(TCCalibEnergy)
