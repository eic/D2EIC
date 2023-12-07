///////////////////////////////////////////
// Check roman pot detector acceptance
// Use single particle simulation
// 12/07/2023 Jihee Kim (jkim11@bnl.gov)
///////////////////////////////////////////
#include <string>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include <stdio.h>
#include <stdlib.h>
// Defining constant parameter
#define IP8_CROSSING_ANGLE 0.035 // unit [rad]
#define RPSFSt1_X 100.850 // [cm]
#define RPSFSt2_X 104.100 // [cm]

void plot_ip8_RPSF_acceptance_edm4hep_output()
{
  // Setting for figures
  TStyle* kStyle = new TStyle("kStyle","Kim's Style");
  kStyle->SetOptStat(0); kStyle->SetOptTitle(0); kStyle->SetOptFit(1);
  kStyle->SetStatColor(0); kStyle->SetStatW(0.15); kStyle->SetStatH(0.10); kStyle->SetStatX(0.85); kStyle->SetStatY(0.9);
  kStyle->SetStatBorderSize(1); kStyle->SetTitleBorderSize(0); kStyle->SetTitleFillColor(-1); kStyle->SetTitleStyle(0);
  kStyle->SetLabelSize(0.045,"xyz"); kStyle->SetTitleSize(0.050,"xyz");
  kStyle->SetTitleOffset(1.2,"x"); kStyle->SetTitleOffset(1.3,"y"); kStyle->SetTitleOffset(1.2,"z");
  kStyle->SetLineWidth(2);
  kStyle->SetTitleFont(42,"xyz"); kStyle->SetLabelFont(42,"xyz");
  kStyle->SetCanvasDefW(1600); kStyle->SetCanvasDefH(800); kStyle->SetCanvasColor(0);
  kStyle->SetPadTickX(1); kStyle->SetPadTickY(1); kStyle->SetPadGridX(1); kStyle->SetPadGridY(1);
  kStyle->SetPadLeftMargin(0.15); kStyle->SetPadRightMargin(0.15); kStyle->SetPadTopMargin(0.1); kStyle->SetPadBottomMargin(0.15);
  TGaxis::SetMaxDigits(3);
  gStyle->SetPalette(1);
  gROOT->SetStyle("kStyle");
  // Reading input files
  TChain *chain = new TChain("events");
  chain->Add("/gpfs02/eic/jkim/condorFullSim/simfiles_ff_rpsf/*.root");
  TTreeReader myReader(chain);
  // Accessing variables
  TTreeReaderArray<Int_t>    pdg        = {myReader, "MCParticles.PDG"}; 
  TTreeReaderArray<Float_t>  px_mc      = {myReader, "MCParticles.momentum.x"}; 
  TTreeReaderArray<Float_t>  py_mc      = {myReader, "MCParticles.momentum.y"}; 
  TTreeReaderArray<Float_t>  pz_mc      = {myReader, "MCParticles.momentum.z"}; 
  TTreeReaderArray<Double_t> rpsf_pos_x = {myReader, "IP8SimpleRomanPot_Hits.position.x"};
  TTreeReaderArray<Double_t> rpsf_pos_y = {myReader, "IP8SimpleRomanPot_Hits.position.y"};
  TTreeReaderArray<Double_t> rpsf_pos_z = {myReader, "IP8SimpleRomanPot_Hits.position.z"};
  // Setting ranges for figures
  Int_t    nbins      = 100;
  Double_t etarange   = 15.;
  Double_t ptrange    = 2.;
  Double_t thetarange = 6.5;
  Double_t phirange   = 3.14;
  // Defining histograms
  TH1D *hThetaMC       = new TH1D("hThetaMC",       "; #theta_{mc} [mrad]; Events",                    nbins, 0., thetarange);
  TH1D *hPhiMC         = new TH1D("hPhiMC",         "; #phi_{mc} [mrad]; Events",                      nbins, -phirange, phirange);
  TH1D *hEtaMC         = new TH1D("hEtaMC",         "; #eta_{mc}; Events",                             nbins, 5., etarange);
  TH1D *hPtMC          = new TH1D("hPtMC",          "; P_{Tmc} [GeV]; Events",                         nbins, 0., ptrange);
  TH2D *hThetaVsPhiMC  = new TH2D("hThetaVsPhiMC",  "; #theta_{mc} [mrad]; #phi_{mc} [rad]",           nbins, 0., thetarange, nbins,-phirange, phirange);
  TH1D *hAccRPSFTheta  = new TH1D("hAccRPSFTheta",  "; #theta_{mc} [mrad]; Events",                    nbins, 0., thetarange);
  TH1D *hAccRPSFPhi    = new TH1D("hAccRPSFPhi",    "; #phi_{mc} [mrad]; Events",                      nbins, -phirange, phirange);
  TH1D *hAccRPSFEta    = new TH1D("hAccRPSFEta",    "; #eta_{mc}; Events",                             nbins, 5., etarange);
  TH1D *hAccRPSFPt     = new TH1D("hAccRPSFPt",     "; P_{Tmc} [GeV]; Events",                         nbins, 0., ptrange);
  TH2D *hAccRPSF       = new TH2D("hAccRPSF",       "; RPSF #theta_{mc} [mrad]; RPSF #phi_{mc} [rad]", nbins, 0., thetarange, nbins,-phirange, phirange);
  TH2D *hAccRPSFPosSt1 = new TH2D("hAccRPSFPosSt1", "; RPSF1 X position [cm]; RPSF1 Y position [cm]",  700, 85., 120., 320, -8., 8.);
  TH2D *hAccRPSFPosSt2 = new TH2D("hAccRPSFPosSt2", "; RPSF2 X position [cm]; RPSF2 Y position [cm]",  700, 85., 120., 320, -8., 8.);
  // Calculating errors  
  hThetaMC->Sumw2(); hPhiMC->Sumw2(); hEtaMC->Sumw2(); hPtMC->Sumw2();
  hAccRPSF->Sumw2(); hAccRPSFTheta->Sumw2(); hAccRPSFPhi->Sumw2(); hAccRPSFEta->Sumw2(); hAccRPSFPt->Sumw2(); 
  // Counting number of events
  Int_t ievt = 0; 
  Int_t iaccevt = 0; 
  // Reading event by event
  while (myReader.Next()) 
  {
    // Showing processing status          
    ievt++;
    if (ievt % 10000 == 0)
      std::cout << "Processing " << ievt << std::endl;
    // Transforming vector properly with crossing angle
    TVector3 mc(px_mc[0], py_mc[0], pz_mc[0]);
    Double_t s = std::sin(IP8_CROSSING_ANGLE);
    Double_t c = std::cos(IP8_CROSSING_ANGLE);
    Double_t pz = (c * mc.Z()) - (s * mc.X());
    Double_t px = (s * mc.Z()) - (c * mc.X());
    TVector3 transmc(px, mc.Y(), pz);
    Double_t theta = transmc.Theta() * 1000.; // unit [mrad]
    Double_t eta = transmc.Eta();
    Double_t phi = transmc.Phi();
    Double_t pt = transmc.Pt();
    // Filling histogram for generated events
    hThetaMC->Fill(theta);
    hPhiMC->Fill(phi);
    hEtaMC->Fill(eta);
    hPtMC->Fill(pt);
    hThetaVsPhiMC->Fill(theta,phi);
    // Checking for accepted events
    if (rpsf_pos_x.GetSize() != 0)
    {
      vector<double> rpsf;
      Int_t rpsf_count = 0;
      for (Int_t i = 0; i < rpsf_pos_x.GetSize(); i++)
      {
	rpsf.push_back(TMath::Floor(rpsf_pos_z[i]/10.));
        if (rpsf_pos_z[i] < 45000.)
          hAccRPSFPosSt1->Fill(rpsf_pos_x[i]/10.,rpsf_pos_y[i]/10.);
	else
          hAccRPSFPosSt2->Fill(rpsf_pos_x[i]/10.,rpsf_pos_y[i]/10.);
      }
      Double_t rpsf_1st = 4392.;
      Double_t rpsf_2nd = 4542.;
      if (std::find(rpsf.begin(), rpsf.end(), rpsf_1st) != rpsf.end()) rpsf_count++;
      if (std::find(rpsf.begin(), rpsf.end(), rpsf_2nd) != rpsf.end()) rpsf_count++;
      // Requring both two stations have hits to be accepted for reconstruction purpose
      if (rpsf_count == 2)
      {
        iaccevt++;
        hAccRPSF->Fill(theta,phi);
        hAccRPSFTheta->Fill(theta);
        hAccRPSFPhi->Fill(phi);
        hAccRPSFEta->Fill(eta);
        hAccRPSFPt->Fill(pt);
      }
    }
  }
  // Printing number of events
  std::cout << "*** Total Events: " << ievt << std::endl;
  std::cout << "*** Total Accepted Events: " << iaccevt << std::endl;
  std::cout << "*** Total Accepted Events [%]: " << ((double)iaccevt/(double)ievt) * 100. << std::endl;
  // Making figures
  // Figures of generated and accepted 2D histograms
  TCanvas* c1 = new TCanvas("c1","c1");
  c1->Divide(2,1);
  c1->cd(1);
  hThetaVsPhiMC->GetXaxis()->CenterTitle(true);
  hThetaVsPhiMC->GetYaxis()->CenterTitle(true);
  hThetaVsPhiMC->Draw("COLZ");
  c1->cd(2);
  hAccRPSF->GetXaxis()->CenterTitle(true);
  hAccRPSF->GetYaxis()->CenterTitle(true);
  hAccRPSF->Draw("COLZ");
  c1->SaveAs("./figures_RPSF/hAccRPSF.png"); 
  c1->SaveAs("./figures_RPSF/hAccRPSF.pdf");
  // Figures of generated accepted 1D histograms in theta and phi
  TCanvas* c2 = new TCanvas("c2","c2");
  c2->Divide(2,1);
  c2->cd(1);
  hThetaMC->GetXaxis()->CenterTitle(true);
  hThetaMC->GetYaxis()->CenterTitle(true);
  hThetaMC->SetLineWidth(2);  
  hThetaMC->Draw("HIST");
  hAccRPSFTheta->SetMarkerStyle(24);
  hAccRPSFTheta->Draw("PE SAME");
  TLegend* L1 = new TLegend(0.23, 0.58, 0.33, 0.73);
  L1->SetLineColor(kWhite); L1->SetFillColor(0); L1->SetTextSize(25); L1->SetTextFont(45);
  L1->AddEntry((TObject*)0, "RPSF", "");
  L1->AddEntry(hThetaMC, "MC  ", "L");
  L1->AddEntry(hAccRPSFTheta, "ACCEPTED", "PE");
  L1->Draw("SAME");  
  c2->cd(2);
  hPhiMC->GetXaxis()->CenterTitle(true);
  hPhiMC->GetYaxis()->CenterTitle(true);
  hPhiMC->GetYaxis()->SetRangeUser(0., hPhiMC->GetMaximum() * 1.1);
  hPhiMC->SetLineWidth(2);  
  hPhiMC->Draw("HIST");
  hAccRPSFPhi->SetMarkerStyle(24);
  hAccRPSFPhi->Draw("PE SAME");
  TLegend* L2 = new TLegend(0.23, 0.58, 0.33, 0.73);
  L2->SetLineColor(kWhite); L2->SetFillColor(0); L2->SetTextSize(25); L2->SetTextFont(45);
  L2->AddEntry((TObject*)0, "RPSF", "");
  L2->AddEntry(hPhiMC, "MC  ", "L");
  L2->AddEntry(hAccRPSFPhi, "ACCEPTED", "PE");
  L2->Draw("SAME");    
  c2->SaveAs("./figures_RPSF/hAccRPSFThetaPhi.png"); 
  c2->SaveAs("./figures_RPSF/hAccRPSFThetaPhi.pdf");
  // Figures of generated accepted 1D histograms in eta and pT
  TCanvas* c3 = new TCanvas("c3","c3");
  c3->Divide(2,1);
  c3->cd(1);
  gPad->SetLogy(1);
  hEtaMC->GetXaxis()->CenterTitle(true);
  hEtaMC->GetYaxis()->CenterTitle(true);
  hEtaMC->GetYaxis()->SetRangeUser(0.1, hEtaMC->GetMaximum() * 2.1);
  hEtaMC->SetLineWidth(2);
  hEtaMC->Draw("HIST");
  hAccRPSFEta->SetMarkerStyle(24);
  hAccRPSFEta->Draw("PE SAME");
  TLegend* L3 = new TLegend(0.6, 0.60, 0.76, 0.75);
  L3->SetLineColor(kWhite); L3->SetFillColor(0); L3->SetTextSize(25); L3->SetTextFont(45);
  L3->AddEntry((TObject*)0, "RPSF", "");
  L3->AddEntry(hEtaMC, "MC  ", "L");
  L3->AddEntry(hAccRPSFEta, "ACCEPTED", "PE");
  L3->Draw("SAME");    
  c3->cd(2);
  hPtMC->GetXaxis()->CenterTitle(true);
  hPtMC->GetYaxis()->CenterTitle(true);
  hPtMC->SetLineWidth(2);
  hPtMC->Draw("HIST");
  hAccRPSFPt->SetMarkerStyle(24);
  hAccRPSFPt->Draw("PE SAME");
  TLegend* L4 = new TLegend(0.23, 0.58, 0.33, 0.73);
  L4->SetLineColor(kWhite); L4->SetFillColor(0); L4->SetTextSize(25); L4->SetTextFont(45);
  L4->AddEntry((TObject*)0, "RPSF", "");
  L4->AddEntry(hPtMC, "MC  ", "L");
  L4->AddEntry(hAccRPSFPt, "ACCEPTED", "PE");
  L4->Draw("SAME");    
  c3->SaveAs("./figures_RPSF/hAccRPSFEtaPt.png"); 
  c3->SaveAs("./figures_RPSF/hAccRPSFEtaPt.pdf"); 
  // Figures of hit positions in each roman pot 
  TCanvas* c4 = new TCanvas("c4","c4");
  c4->Divide(2,1);
  c4->cd(1);
  gPad->SetLogz(1);
  hAccRPSFPosSt1->GetXaxis()->CenterTitle(true);
  hAccRPSFPosSt1->GetYaxis()->CenterTitle(true);
  hAccRPSFPosSt1->GetXaxis()->SetNdivisions(505);
  hAccRPSFPosSt1->Draw("COLZ");
  c4->cd(2);
  gPad->SetLogz(1);
  hAccRPSFPosSt2->GetXaxis()->CenterTitle(true);
  hAccRPSFPosSt2->GetYaxis()->CenterTitle(true);
  hAccRPSFPosSt2->GetXaxis()->SetNdivisions(505);
  hAccRPSFPosSt2->Draw("COLZ");
  c4->SaveAs("./figures_RPSF/hAccRPSFPos.png"); 
  c4->SaveAs("./figures_RPSF/hAccRPSFPos.pdf"); 
} 
