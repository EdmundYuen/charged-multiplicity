#ifndef HISTO_PLOTTING_CPP_INCLUDED
#define HISTO_PLOTTING_CPP_INCLUDED



#endif // HISTO_PLOTTING_CPP_INCLUDED

#include "TROOT.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <algorithm>
#include "math.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TLegend.h"
#include "TText.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TAttAxis.h"
#include "TString.h"
#include "TChain.h"
#include "TChainElement.h"

using namespace std;
using namespace ROOT::Math;

void histo_plotting()
{
    TFile *data_plot = TFile::Open("Histos/data_trees.root");
    TFile *herwig_plot = TFile::Open("Histos/data_treesCUETP8M1_66.root");
    gStyle->SetOptStat(0);
    //gStyle->SetTitleAlign(33);

    //============================================= Eta ===============================================================

    TH1F *data_eta_histo = (TH1F*)data_plot->Get("data eta");
    TH1F *herwig_eta_histo = (TH1F*)herwig_plot->Get("reco eta");
    TCanvas *eta_canvas = new TCanvas ("eta canvas", "Eta Working Plots", 2);
    TLegend *legend_eta = new TLegend (0.6, 0.6, 0.8, 0.8);
    eta_canvas->cd();
    data_eta_histo->GetXaxis()->SetRangeUser(-2.4, 2.4);
    data_eta_histo->SetMinimum(0);
    data_eta_histo->SetMaximum(0.04);
    data_eta_histo->SetTitle("Ongoing Analysis #sqrt{s} = 13TeV");
    data_eta_histo->GetXaxis()->SetLabelSize(0.02);
    data_eta_histo->GetYaxis()->SetTitleOffset(2.0);
    data_eta_histo->GetYaxis()->SetTitleSize(0.02);
    data_eta_histo->GetYaxis()->SetLabelOffset(0.001);
    data_eta_histo->GetYaxis()->SetLabelSize(0.02);
    data_eta_histo->SetLineStyle(0);
    data_eta_histo->SetLineColorAlpha(kBlack,1);
    data_eta_histo->SetMarkerStyle(20);
    data_eta_histo->Draw("E1");

    herwig_eta_histo->SetMinimum(0);
    herwig_eta_histo->SetMaximum(0.04);
    herwig_eta_histo->SetTitle("Ongoing Analysis #sqrt{s} = 13TeV");
    herwig_eta_histo->SetLineWidth(3);
    herwig_eta_histo->GetXaxis()->SetRangeUser(-2.4, 2.4);
    herwig_eta_histo->GetXaxis()->SetLabelSize(0.02);
    herwig_eta_histo->GetYaxis()->SetTitleOffset(2.0);
    herwig_eta_histo->GetYaxis()->SetTitleSize(0.02);
    herwig_eta_histo->GetYaxis()->SetLabelOffset(0.001);
    herwig_eta_histo->GetYaxis()->SetLabelSize(0.02);
    //herwig_eta_histo->SetFillColor(kRed);
    herwig_eta_histo->SetLineColorAlpha(kRed,0.8);
    //herwig_eta_histo->SetMarkerStyle(20);
    herwig_eta_histo->Draw("SAME HIST");

    legend_eta->SetFillColor(0);
    legend_eta->SetFillStyle(0);
    legend_eta->SetBorderSize(0);
    legend_eta->SetTextSize(0.02);
    legend_eta->AddEntry(data_eta_histo, "Uncorrected Data #eta", "pe");
    legend_eta->AddEntry(herwig_eta_histo, "Herwig #eta", "l");
    legend_eta->Draw();

    eta_canvas->Update();
    eta_canvas->SaveAs("eta.png");

    //============================================= Phi ===============================================================

    TH1F *data_phi_histo = (TH1F*)data_plot->Get("data phi");
    TH1F *herwig_phi_histo = (TH1F*)herwig_plot->Get("reco phi");
    TCanvas *phi_canvas = new TCanvas ("phi canvas", "Phi Working Plot", 2);
    TLegend *leg_phi = new TLegend (0.6, 0.6, 0.8, 0.8);
    phi_canvas->cd();
    data_phi_histo->SetMinimum(0);
    data_phi_histo->SetMaximum(0.04);
    data_phi_histo->SetTitle("Ongoing Analysis #sqrt{s} = 13TeV");
    data_phi_histo->GetXaxis()->SetLabelSize(0.02);
    data_phi_histo->GetYaxis()->SetTitleOffset(1.9);
    data_phi_histo->GetYaxis()->SetTitleSize(0.02);
    data_phi_histo->GetYaxis()->SetLabelOffset(0.001);
    data_phi_histo->GetYaxis()->SetLabelSize(0.02);
    data_phi_histo->SetLineStyle(0);
    data_phi_histo->SetLineColorAlpha(kBlack,1);
    data_phi_histo->SetMarkerStyle(20);
    data_phi_histo->Draw("E1");

    herwig_phi_histo->SetMinimum(0);
    herwig_phi_histo->SetMaximum(0.04);
    herwig_phi_histo->SetTitle("Ongoing Analysis #sqrt{s} = 13TeV");
    herwig_phi_histo->SetLineWidth(3);
    herwig_phi_histo->GetXaxis()->SetRangeUser(-3, 3);
    herwig_phi_histo->GetXaxis()->SetLabelSize(0.02);
    herwig_phi_histo->GetYaxis()->SetTitleOffset(1.9);
    herwig_phi_histo->GetYaxis()->SetTitleSize(0.02);
    herwig_phi_histo->GetYaxis()->SetLabelOffset(0.001);
    herwig_phi_histo->GetYaxis()->SetLabelSize(0.02);
    herwig_phi_histo->SetLineStyle(0);
    herwig_phi_histo->SetLineColorAlpha(kRed,0.8);
    herwig_phi_histo->Draw("SAME HIST");

    leg_phi->SetFillColor(0);
    leg_phi->SetFillStyle(0);
    leg_phi->SetBorderSize(0);
    leg_phi->SetTextSize(0.02);
    leg_phi->AddEntry(data_phi_histo, "Uncorrected Data #phi", "pe");
    leg_phi->AddEntry(herwig_phi_histo, "Herwig #phi", "l");
    leg_phi->Draw();

    phi_canvas->Update();
    phi_canvas->SaveAs("phi.png");

    //============================================= pT ===============================================================

    TH1F *data_pt_histo = (TH1F*)data_plot->Get("data pT");
    TH1F *herwig_pt_histo = (TH1F*)herwig_plot->Get("reco pT");
    TCanvas *pt_canvas = new TCanvas ("pt canvas", "pt Working Plot", 4);
    TLegend *leg_pt = new TLegend (0.6, 0.6, 0.8, 0.8);

    pt_canvas->cd();
    gPad->SetLogy();
    data_pt_histo->SetTitle("Ongoing Analysis #sqrt{s} = 13TeV");
    data_pt_histo->GetXaxis()->SetTitle("p_{T} (GeV)");
    data_pt_histo->GetXaxis()->SetLabelSize(0.02);
    data_pt_histo->GetYaxis()->SetTitleOffset(1.3);
    data_pt_histo->GetYaxis()->SetTitleSize(0.03);
    data_pt_histo->GetYaxis()->SetLabelOffset(0.001);
    data_pt_histo->GetYaxis()->SetLabelSize(0.02);
    data_pt_histo->SetLineStyle(0);
    data_pt_histo->SetLineColorAlpha(kBlack,1);
    data_pt_histo->SetMarkerStyle(20);
    data_pt_histo->Draw("E1");

    herwig_pt_histo->SetTitle("Ongoing Analysis #sqrt{s} = 13TeV");
    herwig_pt_histo->SetLineWidth(3);
    herwig_pt_histo->GetXaxis()->SetTitle("p_{T} (GeV)");
    herwig_pt_histo->GetXaxis()->SetLabelSize(0.02);
    herwig_pt_histo->GetYaxis()->SetTitleOffset(1.3);
    herwig_pt_histo->GetYaxis()->SetTitleSize(0.03);
    herwig_pt_histo->GetYaxis()->SetLabelOffset(0.001);
    herwig_pt_histo->GetYaxis()->SetLabelSize(0.02);
    herwig_pt_histo->SetLineColorAlpha(kRed,0.8);
    herwig_pt_histo->Draw("SAME HIST");

    leg_pt->SetFillColor(0);
    leg_pt->SetFillStyle(0);
    leg_pt->SetBorderSize(0);
    leg_pt->SetTextSize(0.02);
    leg_pt->AddEntry(data_pt_histo, "Uncorrected Data p_{T}", "pe");
    leg_pt->AddEntry(herwig_pt_histo, "Herwig p_{T}", "l");
    leg_pt->Draw();

    pt_canvas->Update();
    pt_canvas->SaveAs("pt.png");

    //============================================= Multiplicity ===============================================================

    TH1F *data_multiplicity = (TH1F*)data_plot->Get("Normalized_Multiplicity");
    TH1F *herwig_multiplicity = (TH1F*)herwig_plot->Get("Normalized_Multiplicity");
    TCanvas *multiplicity_canvas = new TCanvas ("multiplicity canvas", "Multiplicity Working Plot", 2);
    TLegend *leg_multiplicity = new TLegend (0.6, 0.6, 0.8, 0.8);

    multiplicity_canvas->cd();
    gPad->SetLogy();
    cout << "Getting data multiplicity" << endl;
    data_multiplicity->SetTitle("Ongoing Analysis #sqrt{s} = 13TeV");
    data_multiplicity->GetXaxis()->SetLabelSize(0.02);
    data_multiplicity->GetXaxis()->SetRangeUser(0, 140);
    data_multiplicity->GetYaxis()->SetTitleOffset(1.3);
    data_multiplicity->GetYaxis()->SetTitleSize(0.03);
    data_multiplicity->GetYaxis()->SetLabelOffset(0.001);
    data_multiplicity->GetYaxis()->SetLabelSize(0.02);
    data_multiplicity->GetYaxis()->SetTitle("#frac{1}{N_{ev}} #frac{dN_{ch}}{dN_{ev}}");
    data_multiplicity->SetLineStyle(0);
    data_multiplicity->SetLineColorAlpha(kBlack,1);
    data_multiplicity->SetMarkerStyle(20);
    data_multiplicity->SetMarkerColor(kBlack);
    data_multiplicity->Draw("E1");

    herwig_multiplicity->SetTitle("Ongoing Analysis #sqrt{s} = 13TeV");
    cout << "Getting Herwig multiplicity" << endl;
    herwig_multiplicity->SetLineWidth(3);
    herwig_multiplicity->GetXaxis()->SetLabelSize(0.02);
    herwig_multiplicity->GetYaxis()->SetTitleOffset(1.3);
    herwig_multiplicity->GetYaxis()->SetTitleSize(0.03);
    herwig_multiplicity->GetYaxis()->SetLabelOffset(0.001);
    herwig_multiplicity->GetYaxis()->SetLabelSize(0.02);
    herwig_multiplicity->GetYaxis()->SetTitle("#frac{1}{N_{ev}} #frac{dN_{ch}}{dN_{ev}}");
    herwig_multiplicity->SetLineColorAlpha(kRed,0.8);
    herwig_multiplicity->Draw("SAME HIST");

    leg_multiplicity->SetFillColor(0);
    leg_multiplicity->SetFillStyle(0);
    leg_multiplicity->SetBorderSize(0);
    leg_multiplicity->SetTextSize(0.02);
    leg_multiplicity->AddEntry(data_multiplicity, "Uncorrected Data", "pe");
    leg_multiplicity->AddEntry(herwig_multiplicity, "Herwig" "l");
    leg_multiplicity->Draw();

    multiplicity_canvas->Update();
    multiplicity_canvas->SaveAs("multiplicity.png");

}
