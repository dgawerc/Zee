
#include <iostream>
#include <fstream> 
#include <sstream>
#include <map>
#include <utility>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TLatex.h"
#include <math.h>
#include <time.h>
#include <algorithm>
#include <functional>
#include "TRandom3.h"
#include "TError.h"
#include <cmath>


TF1* Fitter(TH1 *hist) {
  // Helper function for fitting DOUBLE Gaussian
  float xmin = hist->GetMean() - 3.0*hist->GetRMS();
  float xmax = hist->GetMean() + 3.0*hist->GetRMS();
  TF1 *DGfit = new TF1("DGfit", "gaus(0) + gaus(3)", xmin, xmax);
  DGfit->SetParameters(hist->GetMaximum()/2.0, 0, hist->GetRMS(),
                       hist->GetMaximum()/2.0, 0, hist->GetRMS());
  DGfit->SetParLimits(1, -1.0, 1.0 );
  DGfit->SetParLimits(4, -1.0, 1.0 );
  DGfit->SetParLimits(2, .001, 3.0 );
  DGfit->SetParLimits(5, .001, 3.0 );
  DGfit->SetParNames("A_{1}", "#mu_{1}", "#sigma_{1}", "A_{2}", "#mu_{2}", "#sigma_{2}");
  gPrintViaErrorHandler = kTRUE; // Deal with MINOS print statements in error handler (for next line to be useful)
  gErrorIgnoreLevel = kWarning; //Suppress MINOS new minimum message
  hist->Fit("DGfit","Q","",xmin,xmax); // Q suppresses fit results. Was QMLES
  gErrorIgnoreLevel = kInfo;
  gStyle->SetOptFit(1);
  return hist->GetFunction("DGfit");
  // Call output with, e.g. ->GetParameter(2)
}



void scatterPlot(float x[], float y[], float ex[], float ey[], int points, TH1F *hist[], string xTitle, string yTitle, string Title, string filename) {
  // Copy input arrays, so changes won't affect the initial arrays
  float X[points];
  float Y[points];
  float eX[points];
  float eY[points];
  std::copy(x, x+points, X);
  std::copy(y, y+points, Y);
  std::copy(ex, ex+points, eX);
  std::copy(ey, ey+points, eY);

  // Correct For Low-Event Histograms:
  int  nBad = 0;
  for (int i=0; i<points; i++) {
    if (hist[i]->GetEntries() < 1E3) {
      for (int j=i; j<(points-1); j++){
        X[j-nBad] = X[j-nBad+1]; // delete bad element and shift others left
        Y[j-nBad] = Y[j-nBad+1]; // Note last element gets duplicated but also stays where it is
        eX[j-nBad]=eX[j-nBad+1]; // This is ok, since we will plot the first (points-nBad) points
        eY[j-nBad]=eY[j-nBad+1];
      }
      nBad +=1;
    }
  }

  // Make the Plot:
  TCanvas *c = new TCanvas("c","c",200,10,600,400);
  c->SetLeftMargin(0.14);
  c->SetBottomMargin(0.15);
  c->SetTopMargin(0.12);
  c->cd();

  TGraphErrors* ge = new TGraphErrors(points-nBad, X, Y, eX, eY);
  ge->Draw("ap");

  TAxis* Xaxis = ge->GetXaxis();
  Xaxis->SetTitle( xTitle.c_str() );
  Xaxis->SetTitleSize(0.08);
  Xaxis->SetLabelSize(0.06);
  Xaxis->SetTitleOffset(0.8);

  TAxis* Yaxis = ge->GetYaxis();
  Yaxis->SetTitle( yTitle.c_str() );
  Yaxis->SetTitleSize(0.08);
  Yaxis->SetLabelSize(0.06);
  Yaxis->SetTitleOffset(0.85);

  gStyle->SetTitleFontSize(0.1);
  ge->SetTitle( Title.c_str() );
  ge->SetLineWidth(3);

  gErrorIgnoreLevel = kWarning; //Suppress info about created file
  c->SaveAs( ("plots/pdf/"+filename+".pdf").c_str() );
  c->SaveAs( ("plots/png/"+filename+".png").c_str() );
  c->SaveAs( ("plots/c/"+filename+".C").c_str() );
  delete c;
  gErrorIgnoreLevel = kInfo;
}



float GetDGSigma(TF1 *fit) {
  float p[6];
  for (int i=0; i<6; i++) p[i] = fit->GetParameter(i);
  return (p[0] * p[2] + p[3] * p[5]) / (p[0] + p[3]);
}



float GetDGMean(TF1 *fit) {
  float p[6];
  for (int i=0; i<6; i++) p[i] = fit->GetParameter(i);
  return (p[0] * p[1] + p[3] * p[4]) / (p[0] + p[3]);
}



float GetDGSigmaError(TF1 *fit) {
  float p[6];
  float e[6];
  for (int i=0; i<6; i++) {
    p[i] = fit->GetParameter(i);
    e[i] = fit->GetParError(i);
  }
  return sqrt(pow(p[0]*e[2], 2) + pow(p[3]*e[5], 2) + pow( (p[2]-p[5])/(p[0]+p[3]), 2) * ( pow(p[3]*e[0], 2) + pow(p[0]*e[3], 2) ))
         / (p[0] + p[3]);
}



float GetDGMeanError(TF1 *fit) {
  float p[6];
  float e[6];
  for (int i=0; i<6; i++) {
    p[i] = fit->GetParameter(i);
    e[i] = fit->GetParError(i);
  }
  return sqrt(pow(p[0]*e[1], 2) + pow(p[3]*e[4], 2) + pow( (p[1]-p[4])/(p[0]+p[3]), 2) * ( pow(p[3]*e[0], 2) + pow(p[0]*e[3], 2) ))
         / (p[0] + p[3]);
}



template <class T1, class T2>
void SigMeanTGraph(TFile *file, TH1F** hists, string histName, T1 cutArray[], T2 stepSize, int nSteps, string Time, string xTitle, string outFile) {
  TF1 *Fit[nSteps];
  float  Sigma[nSteps];
  float  SigEr[nSteps];
  float   Mean[nSteps];
  float MeanEr[nSteps];

  for (int i=0; i<nSteps; i++) {
    if ( !(hists[i]->GetEntries() != 0) ) continue;
    Fit[i] = Fitter(hists[i]);
    file->WriteTObject(hists[i], ( histName+Form("[%d]",i) ).c_str(),"WriteDelete");
    Sigma[i] = GetDGSigma( Fit[i] );      //Fit[i]->GetParameter(2);
    SigEr[i] = GetDGSigmaError( Fit[i] ); //Fit[i]->GetParError(2);
    Mean[i]  = GetDGMean( Fit[i] );       //Fit[i]->GetParameter(1);
    MeanEr[i]= GetDGMeanError( Fit[i] );  //Fit[i]->GetParError(1);
  }

  float xVals[nSteps];
  float xErr[nSteps];
  for(int i=0; i<nSteps; i++) xVals[i] = cutArray[i] + stepSize/2;
  std::fill(xErr, xErr+nSteps, float(stepSize)/2);

  scatterPlot(xVals, Sigma, xErr, SigEr, nSteps, hists, /*hists[0]->GetXaxis()->GetTitle()*/ xTitle.c_str(), "Time Resolution (ns)", ("#sigma of t_{1}-t_{2} "+Time).c_str() , (outFile+"_Sigma").c_str() );
  scatterPlot(xVals, Mean, xErr, MeanEr, nSteps, hists, /*hists[0]->GetXaxis()->GetTitle()*/ xTitle.c_str(), "Mean", ("Mean of t_{1}-t_{2} "+Time).c_str() , (outFile+"_Mean").c_str() );
}



void SaveHist(TH1 *hist) {
  TCanvas *c = new TCanvas("c","c",200,10,600,400);
  c->SetLeftMargin(0.14);
  c->SetBottomMargin(0.15);
  c->SetTopMargin(0.12);
  c->SetRightMargin(0.12);
  c->cd();
  
  hist->Draw();

  TAxis* Xaxis = hist->GetXaxis();
  Xaxis->SetTitleSize(0.08);
  Xaxis->SetLabelSize(0.06);
  Xaxis->SetTitleOffset(0.8);

  TAxis* Yaxis = hist->GetYaxis();
  Yaxis->SetTitleSize(0.08);
  Yaxis->SetLabelSize(0.06);
  Yaxis->SetTitleOffset(0.85);

  gStyle->SetTitleFontSize(0.1);
  gStyle->SetOptStat(0);

  gErrorIgnoreLevel = kWarning; //Suppress info about created file
  c->SaveAs( (string("plots/pdf/")+ hist->GetName() +".pdf").c_str() );
  c->SaveAs( (string("plots/png/")+ hist->GetName() +".png").c_str() );
  c->SaveAs( (string("plots/c/")  + hist->GetName() +".C").c_str() );
  delete c;
  gErrorIgnoreLevel = kInfo;
}



void EtaPhi1D( int etaChannels, int phiChannels, TH1F** etaHists, TH1F** phiHists, string timeTitle ) {
  TF1 *etaFit[etaChannels];
  TF1 *phiFit[phiChannels];
  float etaSigma[etaChannels];
  float etaSigEr[etaChannels];
  float phiSigma[phiChannels];
  float phiSigEr[phiChannels];
  float etaMean[etaChannels];
  float etaMeanEr[etaChannels];
  float phiMean[phiChannels];
  float phiMeanEr[phiChannels];

  for (int i=0; i<phiChannels; i++) {
    if ( !(phiHists[i]->GetEntries() != 0) ) continue;
    phiFit[i] = Fitter(phiHists[i]);
    phiSigma[i] = GetDGSigma( phiFit[i] );      //phiFit[i]->GetParameter(2);
    phiSigEr[i] = GetDGSigmaError( phiFit[i] ); //phiFit[i]->GetParError(2);
    phiMean[i]  = GetDGMean( phiFit[i] );       //phiFit[i]->GetParameter(1);
    phiMeanEr[i]= GetDGMeanError( phiFit[i] );  //phiFit[i]->GetParError(1);
  }
  for (int i=0; i<etaChannels; i++) {
    if ( !(etaHists[i]->GetEntries() != 0) ) continue;
    etaFit[i] = Fitter(etaHists[i]);
    etaSigma[i] = GetDGSigma( etaFit[i] );//etaFit[i]->GetParameter(2);
    etaSigEr[i] = GetDGSigmaError( etaFit[i] );//etaFit[i]->GetParError(2);
    etaMean[i]  = GetDGMean( etaFit[i] );//etaFit[i]->GetParameter(1);
    etaMeanEr[i]= GetDGMeanError( etaFit[i] );//etaFit[i]->GetParError(1);
  }
  TH1D* etaSigmaHist = new TH1D( ("EtaSigma_"+timeTitle).c_str() , (timeTitle+";i#eta;Time Resolution #sigma (ns)").c_str() ,etaChannels, -85.5, 85.5);
  TH1D* etaMeanHist  = new TH1D( ("EtaMean_"+timeTitle ).c_str() , (timeTitle+";i#eta;Mean").c_str() ,etaChannels, -85.5, 85.5);
  TH1D* phiSigmaHist = new TH1D( ("PhiSigma_"+timeTitle).c_str() , (timeTitle+";i#phi;Time Resolution #sigma (ns)").c_str() ,phiChannels, 0, 360);
  TH1D* phiMeanHist  = new TH1D( ("PhiMean_"+timeTitle ).c_str() , (timeTitle+";i#phi;Mean").c_str() ,phiChannels, 0, 360);
  for (int i=0; i<etaChannels; i++) {
    if (etaHists[i]->GetEntries() > 100 ) {
      etaSigmaHist->SetBinContent(i, etaSigma[i]);
      etaSigmaHist->SetBinError(i, etaSigEr[i]);
      etaMeanHist->SetBinContent(i, etaMean[i]);
      etaMeanHist->SetBinError(i, etaMeanEr[i]);
    }
  }
  for (int i=0; i<phiChannels; i++) {
    if (phiHists[i]->GetEntries() > 100 ) {
      phiSigmaHist->SetBinContent(i, phiSigma[i]);
      phiSigmaHist->SetBinError(i, phiSigEr[i]);
      phiMeanHist->SetBinContent(i, phiMean[i]);
      phiMeanHist->SetBinError(i, phiMeanEr[i]);
    }
  }

  etaSigmaHist->SetOption("ehist");
  etaMeanHist->SetOption("ehist");
  phiSigmaHist->SetOption("ehist");
  phiMeanHist->SetOption("ehist");

  SaveHist(etaSigmaHist);
  SaveHist(etaMeanHist);
  SaveHist(phiSigmaHist);
  SaveHist(phiMeanHist);
}



template <class U1>
U1 CutArray(U1 cuts[], int nSteps, U1 min, U1 max) {
  U1 stepSize = (max - min) / nSteps;
  for (int i=0; i<nSteps; i++) cuts[i] = min + i * stepSize;
  cuts[nSteps] = max; // Last element outside loop because StepSize might be rounded if int
  return stepSize;
}



void Zee_beta() {
  string filename = "All2016.root";

  TFile *inputfile = TFile::Open(filename.c_str(),"READ");
  TTree *tree = (TTree*)inputfile->Get("ZeeTiming");

  float mass;
  unsigned int run;
  unsigned int nPV;
  unsigned int eventTime;

  float t1;
  float t1_seed;
  float t1raw_seed;
  float t1calib_seed;
  float t1calib_seed_sept;
  int   ele1SeedIEta;
  int   ele1SeedIPhi;
  float ele1Eta;
  float ele1Phi;
  bool  ele1IsEB;
  float ele1Pt;
  float seed1_transpCorr;

  float t2;
  float t2_seed;
  float t2raw_seed;
  float t2calib_seed;
  float t2calib_seed_sept;
  int   ele2SeedIEta;
  int   ele2SeedIPhi;
  float ele2Eta;
  float ele2Phi;
  bool  ele2IsEB;
  float ele2Pt;
  float seed2_transpCorr;

  tree->SetBranchStatus("*",1);

  tree->SetBranchAddress("mass",&mass);
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("nPV",&nPV);
  tree->SetBranchAddress("eventTime",&eventTime);

  tree->SetBranchAddress("t1",&t1);
  tree->SetBranchAddress("t1_seed",&t1_seed);
  tree->SetBranchAddress("t1raw_seed",&t1raw_seed);
  tree->SetBranchAddress("t1calib_seed",&t1calib_seed);
  tree->SetBranchAddress("t1calib_seed_sept",&t1calib_seed_sept);
  tree->SetBranchAddress("ele1SeedIEta",&ele1SeedIEta);
  tree->SetBranchAddress("ele1SeedIPhi",&ele1SeedIPhi);
  tree->SetBranchAddress("ele1Eta",&ele1Eta);
  tree->SetBranchAddress("ele1Phi",&ele1Phi);
  tree->SetBranchAddress("ele1IsEB",&ele1IsEB);
  tree->SetBranchAddress("ele1Pt",&ele1Pt);
  tree->SetBranchAddress("seed1_transpCorr",&seed1_transpCorr);

  tree->SetBranchAddress("t2",&t2);
  tree->SetBranchAddress("t2_seed",&t2_seed);
  tree->SetBranchAddress("t2raw_seed",&t2raw_seed);
  tree->SetBranchAddress("t2calib_seed",&t2calib_seed);
  tree->SetBranchAddress("t2calib_seed_sept",&t2calib_seed_sept);
  tree->SetBranchAddress("ele2SeedIEta",&ele2SeedIEta);
  tree->SetBranchAddress("ele2SeedIPhi",&ele2SeedIPhi);
  tree->SetBranchAddress("ele2Eta",&ele2Eta);
  tree->SetBranchAddress("ele2Phi",&ele2Phi);
  tree->SetBranchAddress("ele2IsEB",&ele2IsEB);
  tree->SetBranchAddress("ele2Pt",&ele2Pt);
  tree->SetBranchAddress("seed2_transpCorr",&seed2_transpCorr);

  Long64_t nentries = tree->GetEntries();
  std::cout<<"Number of Events in Sample: "<<nentries<<std::endl;

  // Find values of min and max
  float transpMin   = 0.4;
  float transpMax   = 3.4;
  unsigned int runMin = 10E7;
  unsigned int runMax = 0;
  unsigned int nPVMin = 1;
  unsigned int nPVMax = nPVMin;
  unsigned int eventTimeMin = 4E9;
  unsigned int eventTimeMax = 0;
  float t1Min = 10;
  float t1Max = -10;
  float t1seedMin = 10;
  float t1seedMax = -10;
  float t1rawseedMin = 10;
  float t1rawseedMax = -10;
  float t1calibseedMin = 10;
  float t1calibseedMax = -10;
  float t1calibseedseptMin = 10;
  float t1calibseedseptMax = -10;
  float t2Min = 10;
  float t2Max = -10;
  float t2seedMin = 10;
  float t2seedMax = -10;
  float t2rawseedMin = 10;
  float t2rawseedMax = -10;
  float t2calibseedMin = 10;
  float t2calibseedMax = -10;
  float t2calibseedseptMin = 10;
  float t2calibseedseptMax = -10;

  for (Long64_t iEntry=0; iEntry<nentries; iEntry++) {
    tree->GetEntry(iEntry);
    
    if      (run < runMin) runMin = run;
    else if (run > runMax) runMax = run; // the "else" should speed it up. If the data order is monotonic, remove it

    if      (eventTime < eventTimeMin) eventTimeMin = eventTime;
    else if (eventTime > eventTimeMax) eventTimeMax = eventTime;

    if (nPV > nPVMax) nPVMax = nPV;

    if (t1 > -5 && t1 < 5) {
      if      (t1 < t1Min) t1Min = t1;
      else if (t1 > t1Max) t1Max = t1;
      if      (t1_seed < t1seedMin) t1seedMin = t1_seed;
      else if (t1_seed > t1seedMax) t1seedMax = t1_seed;
      if      (t1raw_seed < t1rawseedMin) t1rawseedMin = t1raw_seed;
      else if (t1raw_seed > t1rawseedMax) t1rawseedMax = t1raw_seed;
      if      (t1calib_seed < t1calibseedMin) t1calibseedMin = t1calib_seed;
      else if (t1calib_seed > t1calibseedMax) t1calibseedMax = t1calib_seed;
      if      (t1calib_seed_sept < t1calibseedseptMin) t1calibseedseptMin = t1calib_seed_sept;
      else if (t1calib_seed_sept > t1calibseedseptMax) t1calibseedseptMax = t1calib_seed_sept;
    }

    if (t2 > -5 && t2 < 5) {
      if      (t2 < t2Min) t2Min = t2;
      else if (t2 > t2Max) t2Max = t2;
      if      (t2_seed < t2seedMin) t2seedMin = t2_seed;
      else if (t2_seed > t2seedMax) t2seedMax = t2_seed;
      if      (t2raw_seed < t2rawseedMin) t2rawseedMin = t2raw_seed;
      else if (t2raw_seed > t2rawseedMax) t2rawseedMax = t2raw_seed;
      if      (t2calib_seed < t2calibseedMin) t2calibseedMin = t2calib_seed;
      else if (t2calib_seed > t2calibseedMax) t2calibseedMax = t2calib_seed;
      if      (t2calib_seed_sept < t2calibseedseptMin) t2calibseedseptMin = t2calib_seed_sept;
      else if (t2calib_seed_sept > t2calibseedseptMax) t2calibseedseptMax = t2calib_seed_sept;
    }
  }

  // Construct (nSteps) intervals of: Run, Transparency, nPV
  int nSteps = 50;
  unsigned int runCuts[nSteps + 1];
  int stepSize = CutArray(runCuts, nSteps, runMin, runMax);

  unsigned int nPVCuts[nSteps + 1];
  int nPVStepSize = CutArray(nPVCuts, nSteps, nPVMin, nPVMax);

  unsigned int eventTimeCuts[nSteps + 1];
  int eventTimeStepSize = CutArray(eventTimeCuts, nSteps, eventTimeMin, eventTimeMax);

  float transpCuts[nSteps+1];
  float transpStepSize = CutArray(transpCuts, nSteps, transpMin, transpMax);

  float t1Cuts[nSteps + 1];
  float t1StepSize = CutArray(t1Cuts, nSteps, t1Min, t1Max);

  // Declare Histograms
  TH1F *histRun_t[nSteps];
  TH1F *histRun_tseed[nSteps];
  TH1F *histRun_trawseed[nSteps];
  TH1F *histRun_tcalibseed[nSteps];
  TH1F *histRun_tcalibseedsept[nSteps];

  TH1F *histnPV_t[nSteps];

  TH1F *histEventTime_t[nSteps];

  TH1F *histTransp1_t[nSteps];
  TH1F *histTransp2_t[nSteps];
  TH1F *histTransp1_tseed[nSteps];
  TH1F *histTransp2_tseed[nSteps];
  TH1F *histTransp1_trawseed[nSteps];
  TH1F *histTransp2_trawseed[nSteps];
  TH1F *histTransp1_tcalibseed[nSteps];
  TH1F *histTransp2_tcalibseed[nSteps];
  TH1F *histTransp1_tcalibseedsept[nSteps];
  TH1F *histTransp2_tcalibseedsept[nSteps];

  int BinSize = 3;
  int phiChannels = 360;
  int etaChannels = 171;
  int phiBins = phiChannels/BinSize;
  int etaBins = etaChannels/BinSize;
  TH1F *histEtaPhi_t[phiBins][etaBins]; // Phi from 0 to 360 inclusive. Eta from -85 to 85 inclusive. Bins of 3.
  TH1F *histEta_t[etaChannels];
  TH1F *histPhi_t[phiChannels];
  TH1F *histEta_tseed[etaChannels];
  TH1F *histPhi_tseed[phiChannels];
  TH1F *histEta_trawseed[etaChannels];
  TH1F *histPhi_trawseed[phiChannels];
  TH1F *histEta_tcalibseed[etaChannels];
  TH1F *histPhi_tcalibseed[phiChannels];
  TH1F *histEta_tcalibseedsept[etaChannels];
  TH1F *histPhi_tcalibseedsept[phiChannels];

  // Initialize Histograms
  for(int i=0; i<nSteps; i++) {
    histRun_t[i]     = new TH1F( Form("histRun_t[%d]",i) ,";t_{1}-t_{2};Entries", 120, -3, 3);
    histRun_tseed[i] = new TH1F( Form("histRun_tseed[%d]",i) ,";t_{1}-t_{2} seed;Entries", 120, -3, 3);
    histRun_trawseed[i] = new TH1F( Form("histRun_trawseed[%d]",i) ,";t_{1}-t_{2} raw seed;Entries", 120, -3, 3);
    histRun_tcalibseed[i] = new TH1F( Form("histRun_tcalibseed[%d]",i) ,";t_{1}-t_{2} calib seed;Entries", 120, -3, 3);
    histRun_tcalibseedsept[i] = new TH1F( Form("histRun_tcalibseedsept[%d]",i) ,";t_{1}-t_{2} calib seed sept;Entries", 120, -3, 3);

    histnPV_t[i] = new TH1F( Form("histnPV_t[%d]",i), ";t_{1}-t_{2};Entries", 120, -3, 3);

    histEventTime_t[i] = new TH1F( Form("histEventTime[%d]",i), ";t_{1}-t_{2};Entries", 120, -3, 3);

    histTransp1_t[i] = new TH1F( Form("histTransp1_t[%d]",i), ";t_{1}-t_{2};Entries", 120, -3, 3);
    histTransp2_t[i] = new TH1F( Form("histTransp2_t[%d]",i), ";t_{1}-t_{2};Entries", 120, -3, 3);
    histTransp1_tseed[i] = new TH1F( Form("histTransp1_tseed[%d]",i), ";t_{1}-t_{2} seed;Entries", 120, -3, 3);
    histTransp2_tseed[i] = new TH1F( Form("histTransp2_tseed[%d]",i), ";t_{1}-t_{2} seed;Entries", 120, -3, 3);
    histTransp1_trawseed[i] = new TH1F( Form("histTransp1_trawseed[%d]",i), ";t_{1}-t_{2} raw seed;Entries", 120, -3, 3);
    histTransp2_trawseed[i] = new TH1F( Form("histTransp2_trawseed[%d]",i), ";t_{1}-t_{2} raw seed;Entries", 120, -3, 3);
    histTransp1_tcalibseed[i] = new TH1F( Form("histTransp1_tcalibseed[%d]",i), ";t_{1}-t_{2} calib seed;Entries", 120, -3, 3);
    histTransp2_tcalibseed[i] = new TH1F( Form("histTransp2_tcalibseed[%d]",i), ";t_{1}-t_{2} calib seed;Entries", 120, -3, 3);
    histTransp1_tcalibseedsept[i] = new TH1F( Form("histTransp1_tcalibseedsept[%d]",i), ";t_{1}-t_{2} calib seed sept;Entries", 120, -3, 3);
    histTransp2_tcalibseedsept[i] = new TH1F( Form("histTransp2_tcalibseedsept[%d]",i), ";t_{1}-t_{2} calib seed sept;Entries", 120, -3, 3);
  }
  for (int i=0; i<phiBins; i++){
    for (int j=0; j<etaBins; j++){
      histEtaPhi_t[i][j] = new TH1F( Form("histEtaPhi_t[%d][%d]",i,j) ,";t_{1}-t_{2};Entries", 120, -3, 3);
    }
  }
  for (int i=0; i<etaChannels; i++) {
    histEta_t[i] = new TH1F( Form("histEta_t[%d]",i),";t_{1}-t_{2};Entries", 120, -3, 3);
    histEta_tseed[i] = new TH1F( Form("histEta_tseed[%d]",i),";t_{1}-t_{2} seed;Entries", 120, -3, 3);
    histEta_trawseed[i] = new TH1F( Form("histEta_trawseed[%d]",i),";t_{1}-t_{2} raw seed;Entries", 120, -3, 3);
    histEta_tcalibseed[i] = new TH1F( Form("histEtacalibseed_t[%d]",i),";t_{1}-t_{2} calib seed;Entries", 120, -3, 3);
    histEta_tcalibseedsept[i] = new TH1F( Form("histEta_tcalibseedsept[%d]",i),";t_{1}-t_{2} calib seed sept;Entries", 120, -3, 3);
  }
  for (int i=0; i<phiChannels; i++) {
    histPhi_t[i] = new TH1F( Form("histPhi_t[%d]",i),";t_{1}-t_{2};Entries", 120, -3, 3);
    histPhi_tseed[i] = new TH1F( Form("histPhi_tseed[%d]",i),";t_{1}-t_{2} seed;Entries", 120, -3, 3);
    histPhi_trawseed[i] = new TH1F( Form("histPhi_trawseed[%d]",i),";t_{1}-t_{2} raw seed;Entries", 120, -3, 3);
    histPhi_tcalibseed[i] = new TH1F( Form("histPhi_tcalibseed[%d]",i),";t_{1}-t_{2 calib seed};Entries", 120, -3, 3);
    histPhi_tcalibseedsept[i] = new TH1F( Form("histPhi_tcalibseedsept[%d]",i),";t_{1}-t_{2} calib seed sept;Entries", 120, -3, 3);
  }


  // LOOP THROUGH EVERY EVENT
  for (Long64_t iEntry=0; iEntry<nentries; iEntry++) {
    if (iEntry %500000 == 0) cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);

    if( !(ele1Pt>30 && ele2Pt>30 && mass>75 && mass<105 && ele1IsEB && ele2IsEB) ) continue; //Prelim cuts for everything

    //Time resolution vs Run
    for (int i=0; i<nSteps; i++) {
      if( !(run>=runCuts[i] && run<runCuts[i+1]) ) continue;
      histRun_t[i]->Fill( t1-t2 ); 
      histRun_tseed[i]->Fill( t1_seed-t2_seed );
      histRun_trawseed[i]->Fill( t1raw_seed-t2raw_seed );
      histRun_tcalibseed[i]->Fill( t1calib_seed-t2calib_seed );
      histRun_tcalibseedsept[i]->Fill( t1calib_seed_sept-t2calib_seed_sept );
      break; // Don't bother to check other cases in the for loop
    }

    for (int i=0; i<nSteps; i++) {
      if( !(nPV>=nPVCuts[i] && nPV<nPVCuts[i+1]) ) continue;
        histnPV_t[i]->Fill( t1-t2 ); 
    }

    for (int i=0; i<nSteps; i++) {
      if( !(eventTime>=eventTimeCuts[i] && eventTime<eventTimeCuts[i+1]) ) continue;
        histEventTime_t[i]->Fill( t1-t2 );
    }

    //Time resolution vs Transparency
    for (int i=0; i<nSteps; i++) {
      if( !(seed1_transpCorr>=transpCuts[i] && seed1_transpCorr<transpCuts[i+1]) ) continue;
      histTransp1_t[i]->Fill( t1-t2 );
      histTransp1_tseed[i]->Fill( t1_seed-t2_seed );
      histTransp1_trawseed[i]->Fill( t1raw_seed-t2raw_seed );
      histTransp1_tcalibseed[i]->Fill( t1calib_seed-t2calib_seed );
      histTransp1_tcalibseedsept[i]->Fill( t1calib_seed_sept-t2calib_seed_sept );
      break; // Don't bother to check other cases in the for loop
    }
    for (int i=0; i<nSteps; i++) {
      if( !(seed2_transpCorr>=transpCuts[i] && seed2_transpCorr<transpCuts[i+1]) ) continue;
      histTransp2_t[i]->Fill( t1-t2 );
      histTransp2_tseed[i]->Fill( t1_seed-t2_seed );
      histTransp2_trawseed[i]->Fill( t1raw_seed-t2raw_seed );
      histTransp2_tcalibseed[i]->Fill( t1calib_seed-t2calib_seed );
      histTransp2_tcalibseedsept[i]->Fill( t1calib_seed_sept-t2calib_seed_sept );
      break; // Don't bother to check other cases in the for loop
    }

    // Cut on Eta, Phi
    if( abs(ele1Eta)<1.47 && abs(ele2Eta)<1.47 ) { // Electrons in barrel
      histEtaPhi_t[ (ele1SeedIPhi-1)/BinSize ][ (ele1SeedIEta+85)/BinSize ]->Fill(t1-t2);
      histEta_t[ (ele1SeedIEta+85) ]->Fill( t1-t2 );
      histPhi_t[ (ele1SeedIPhi-1)  ]->Fill( t1-t2 );
      histEta_tseed[ (ele1SeedIEta+85) ]->Fill( t1_seed-t2_seed );
      histPhi_tseed[ (ele1SeedIPhi-1)  ]->Fill( t1_seed-t2_seed );
      histEta_trawseed[ (ele1SeedIEta+85) ]->Fill( t1raw_seed-t2raw_seed );
      histPhi_trawseed[ (ele1SeedIPhi-1)  ]->Fill( t1raw_seed-t2raw_seed );
      histEta_tcalibseed[ (ele1SeedIEta+85) ]->Fill( t1calib_seed-t2calib_seed );
      histPhi_tcalibseed[ (ele1SeedIPhi-1)  ]->Fill( t1calib_seed-t2calib_seed );
      histEta_tcalibseedsept[ (ele1SeedIEta+85) ]->Fill( t1calib_seed_sept-t2calib_seed_sept );
      histPhi_tcalibseedsept[ (ele1SeedIPhi-1)  ]->Fill( t1calib_seed_sept-t2calib_seed_sept );
    }
  }

  TFile *file = TFile::Open(("output"+filename).c_str(), "RECREATE");
  file->cd();

  std::cout<<"Generating Run Plots..."<<std::endl;
  SigMeanTGraph(file, histRun_t, "histRun_t", runCuts, stepSize, nSteps, "", "Run Number", "run_t");
  SigMeanTGraph(file, histRun_tseed, "histRun_tseed", runCuts, stepSize, nSteps, "Seed", "Run Number", "run_tseed");
  SigMeanTGraph(file, histRun_trawseed, "histRun_trawseed", runCuts, stepSize, nSteps, "Raw Seed", "Run Number", "run_trawseed");
  SigMeanTGraph(file, histRun_tcalibseed, "histRun_tcalibseed", runCuts, stepSize, nSteps, "Calib Seed", "Run Number", "run_tcalibseed");
  SigMeanTGraph(file, histRun_tcalibseedsept, "histRun_tcalibseedsept", runCuts, stepSize, nSteps, "Calib Seed Sept", "Run Number", "run_tcalibseedsept");


  std::cout<<"Generating nPV Plots..."<<std::endl;
  SigMeanTGraph(file, histnPV_t, "histnPV_t", nPVCuts, nPVStepSize, nSteps, "", "nPV", "nPV_t");


  std::cout<<"Generating EventTime Plots..."<<std::endl;
  SigMeanTGraph(file, histEventTime_t, "histEventTime_t", eventTimeCuts, eventTimeStepSize, nSteps, "", "Event Time", "eventTime_t");


  std::cout<<"Generating Transparency Plots..."<<std::endl;
  SigMeanTGraph(file, histTransp1_t, "histTransp1_t", transpCuts, transpStepSize, nSteps, "", "Seed 1 Transparency", "transp1_t");
  SigMeanTGraph(file, histTransp2_t, "histTransp2_t", transpCuts, transpStepSize, nSteps, "", "Seed 2 Transparency", "transp2_t");
  SigMeanTGraph(file, histTransp1_tseed, "histTransp1_tseed", transpCuts, transpStepSize, nSteps, "Seed", "Seed 1 Transparency", "transp1_tseed");
  SigMeanTGraph(file, histTransp2_tseed, "histTransp2_tseed", transpCuts, transpStepSize, nSteps, "Seed", "Seed 2 Transparency", "transp2_tseed");
  SigMeanTGraph(file, histTransp1_trawseed, "histTransp1_trawseed", transpCuts, transpStepSize, nSteps, "Raw Seed", "Seed 1 Transparency", "transp1_trawseed");
  SigMeanTGraph(file, histTransp2_trawseed, "histTransp2_trawseed", transpCuts, transpStepSize, nSteps, "Raw Seed", "Seed 2 Transparency", "transp2_trawseed");
  SigMeanTGraph(file, histTransp1_tcalibseed, "histTransp1_tcalibseed", transpCuts, transpStepSize, nSteps, "Calib Seed", "Seed 1 Transparency", "transp1_tcalibseed");
  SigMeanTGraph(file, histTransp2_tcalibseed, "histTransp2_tcalibseed", transpCuts, transpStepSize, nSteps, "Calib Seed", "Seed 2 Transparency", "transp2_tcalibseed");
  SigMeanTGraph(file, histTransp1_tcalibseedsept, "histTransp1_tcalibseedsept", transpCuts, transpStepSize, nSteps, "Calib Seed Sept", "Seed 1 Transparency", "transp1_tcalibseedsept");
  SigMeanTGraph(file, histTransp2_tcalibseedsept, "histTransp2_tcalibseedsept", transpCuts, transpStepSize, nSteps, "Calib Seed Sept", "Seed 2 Transparency", "transp2_tcalibseedsept");


  std::cout<<"Generating 1D Eta and Phi Plots..."<<std::endl;
  EtaPhi1D( etaChannels, phiChannels, histEta_t,              histPhi_t,              "t" );
  EtaPhi1D( etaChannels, phiChannels, histEta_tseed,          histPhi_tseed,          "tseed" );
  EtaPhi1D( etaChannels, phiChannels, histEta_trawseed,       histPhi_trawseed,       "trawseed" );
  EtaPhi1D( etaChannels, phiChannels, histEta_tcalibseed,     histPhi_tcalibseed,     "tcalibseed" );
  EtaPhi1D( etaChannels, phiChannels, histEta_tcalibseedsept, histPhi_tcalibseedsept, "tcalibseedsept" );

  // Get gauss fit and save histograms to file
  TF1 *EtaPhiFit_t[phiBins][etaBins];
  float EtaPhiSigma_t[phiBins][etaBins];
  float EtaPhiMean_t[phiBins][etaBins];


  std::cout<<"Starting Gauss Fits for Eta vs Phi Plots..."<<std::endl;
  for (int i=0; i<phiBins; i++){
    for (int j=0; j<etaBins; j++){
      if (i % 10 == 0 && j==0) std::cout<<"Processing Phi Value of "<<i*BinSize<<std::endl;
      if ( !(histEtaPhi_t[i][j]->GetEntries() != 0) ) continue;
      EtaPhiFit_t[i][j] = Fitter(histEtaPhi_t[i][j]);
      EtaPhiSigma_t[i][j] = GetDGSigma( EtaPhiFit_t[i][j] );//EtaPhiFit_t[i][j]->GetParameter(2);
      EtaPhiMean_t[i][j]  = GetDGMean( EtaPhiFit_t[i][j] );//EtaPhiFit_t[i][j]->GetParameter(1);
    }
  }
  
  // Create Eta vs Phi Histograms
  TH2D *histEtaPhiSigma_t = new TH2D("histEtaPhiSigma_t","#sigma; i#phi; i#eta", phiBins,0,360, etaBins,-85.5,85.5); 
  TH2D *histEtaPhiMean_t  = new TH2D("histEtaPhiMean_t","Mean; i#phi; i#eta", phiBins,0,360, etaBins,-85.5,85.5);
  double minSigma = 5E7;
  double minMean = 5E7;
  for (int i=0; i<phiBins; i++) {
    for (int j=0; j<etaBins; j++) {
      // Sigma
      if ( histEtaPhi_t[i][j]->GetEntries() > 150 && EtaPhiSigma_t[i][j]<1) {
        histEtaPhiSigma_t->SetBinContent(i+1, j+1, EtaPhiSigma_t[i][j]);
        if ( EtaPhiSigma_t[i][j]<minSigma ) minSigma = EtaPhiSigma_t[i][j];
      }
      else histEtaPhiSigma_t->SetBinContent(i+1, j+1, -999);

      // Mean
      if ( histEtaPhi_t[i][j]->GetEntries() > 150 && abs(EtaPhiMean_t[i][j])<.2) {
        histEtaPhiMean_t->SetBinContent(i+1, j+1, EtaPhiMean_t[i][j]);
        if ( EtaPhiMean_t[i][j]<minMean ) minMean = EtaPhiMean_t[i][j];
      }
      else histEtaPhiMean_t->SetBinContent(i+1, j+1, -999);
    }
  }
  histEtaPhiSigma_t->SetMinimum(minSigma);
  histEtaPhiMean_t->SetMinimum(minMean);
  histEtaPhiSigma_t->SetOption("colz");
  histEtaPhiMean_t->SetOption("colz");
  SaveHist(histEtaPhiSigma_t);
  SaveHist(histEtaPhiMean_t);
}
