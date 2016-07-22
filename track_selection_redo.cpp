#include "TROOT.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include "TLorentzVector.h"
#include "TSystem.h"
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

using namespace std;
using namespace ROOT::Math;

void track_selection_redo()
{
    //===========================implement cuts================================
    const double eta_cut = 2.;
    const double vtx_number_cut = 1.0;
    const double dist_to_BS = 10;
    const double pt_cut = 0.5;
    const int lumi_cut = 90;
    const float ptErr_pt_cut = 0.05;
    const float dz_dzErr_cut = 3;
    const float d0_d0Err_cut = 3;

    //===========================retrieve ROOT file============================

    TFile *datafile = TFile::Open("treesCUETP8M1_66.root", "READ");
    TTree *datatree = (TTree*)datafile->Get("UETree/data");

    //===========================define variables to read TTree================

    vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *data_tracks = 0;
    datatree->SetBranchAddress("recoTracksp4", &data_tracks);

    int ndata_lumi;
    datatree->SetBranchAddress("lumi", &ndata_lumi);

    int ndata_zerobias;
    datatree->SetBranchAddress("trgZeroBias", &ndata_zerobias);

    vector<float> *fvecdata_vtxx = 0;
    datatree->SetBranchAddress("vtxx", &fvecdata_vtxx);
    int ndata_numberofvtxx = fvecdata_vtxx->size();

    vector<float> *fvecdata_vtxy = 0;
    datatree->SetBranchAddress("vtxy", &fvecdata_vtxy);
    int ndata_numberofvtxy = fvecdata_vtxy->size();

    vector<float> *fvecdata_vtxz = 0;
    datatree->SetBranchAddress("vtxz", &fvecdata_vtxz);
    int ndata_numberofvtxz = fvecdata_vtxz->size();

    vector<float> *fvecdata_vtxxBS = 0;
    datatree->SetBranchAddress("vtxxBS", &fvecdata_vtxxBS);

    vector<float> *fvecdata_vtxyBS = 0;
    datatree->SetBranchAddress("vtxyBS", &fvecdata_vtxyBS);

    vector<float> *fvecdata_vtxzBS = 0;
    datatree->SetBranchAddress("vtxzBS", &fvecdata_vtxzBS);

    vector<int> *nvecdata_highpurity = 0;
    datatree->SetBranchAddress("recoTrackshighPurity", &nvecdata_highpurity);

    vector<int> *nvecdata_vtxndof = 0;
    datatree->SetBranchAddress("vtxndof", &nvecdata_vtxndof);

    vector<float> *fvecdata_dz = 0;
    datatree->SetBranchAddress("recoTracksdz", &fvecdata_dz);

    vector<float> *fvecdata_d0 = 0;
    datatree->SetBranchAddress("recoTracksd0", &fvecdata_d0);

    vector<float> *fvecdata_dzerr = 0;
    datatree->SetBranchAddress("recoTracksdzErr", &fvecdata_dzerr);

    vector<float> *fvecdata_d0err = 0;
    datatree->SetBranchAddress("recoTracksd0Err", &fvecdata_d0err);

    vector<float> *fvecdata_pterr = 0;
    datatree->SetBranchAddress("recoTracksptErr", &fvecdata_pterr);


    //==============================================Histograms==============================================================

    TH1F *data_pt_histo = new TH1F ("data p_{T}", "Normalized data p_{T}", 200, 0, 10);
    TH1F *data_eta_histo = new TH1F("data eta", "Normalised data #eta", 50, -2.5, 2.5);
    TH1F *data_phi_histo = new TH1F ("data phi", "Normalized_data #phi", 60, -3, 3);
    TH1F *data_dz_sigmadz = new TH1F ("data_dz_sigmadz", "data d_{z}/#sigma_{z}", 160, -20, 20);
    TH1F *data_d0_sigmad0 = new TH1F ("data_d0_sigmad0", "data d_{0}/#sigma_{xy}", 160, -20, 20);
    TH1F *data_sigmapt_pt = new TH1F ("data_sigmapt_pt", "data #sigma_{p_{T}}/p_{T}", 20, 0, 0.2);

    //==============================================Variables==============================================================

    float fdata_evt = 0;
    float fdata_trk = 0;
    float fdata_trketa = 0;
    float fdata_trkd0 = 0;
    float fdata_trkpt = 0;
    float fdata_dz_sigmadz, fdata_d0_sigmad0, fdata_sigmapt_pt;
    //float fdata_dz_sigmadz = 0;
    //float fdata_d0_sigmad0 = 0;
    //float fdata_sigmapt_pt = 0;
    int ndata_totaltrk;

    Int_t ndata_totalEvt = (Int_t)datatree->GetEntries();
    cout << "There is a total of " << ndata_totalEvt << " events." << endl;

    for (Int_t i = 0; i < ndata_totalEvt; ++i)
    {
        datatree->GetEntry(i);
        //cout << "At entry " << i << endl;

        if (ndata_lumi >= lumi_cut)
        {
            if (ndata_zerobias == 1)
            {
                ndata_totaltrk = data_tracks->size();
                ++fdata_evt;

                //looping through tracks
                for (int t = 0; t != ndata_totaltrk; ++t)
                {
                    XYZTVector data_vec = (*data_tracks)[t];
                    //cout << "Within track " << endl;

                    fdata_dz_sigmadz = ((*fvecdata_dz)[t])/((*fvecdata_dzerr)[t]);
                    fdata_d0_sigmad0 = ((*fvecdata_d0)[t])/((*fvecdata_d0err)[t]);
                    fdata_sigmapt_pt = (((*fvecdata_pterr)[t])/(data_vec.Pt()));

                    if (abs(data_vec.Eta()) <= eta_cut)
                    {
                        data_pt_histo->Fill(data_vec.Pt());
                        ++fdata_trkpt;
                    }

                    if (data_vec.Pt() >= pt_cut)
                    {
                        data_eta_histo->Fill(data_vec.Eta());
                        ++fdata_trketa;

                    }

                    if (data_vec.Pt() >= pt_cut && abs(data_vec.Eta()) <= eta_cut)
                    {
                        data_phi_histo->Fill(data_vec.Phi());
                        data_d0_sigmad0->Fill(fdata_d0_sigmad0);
                        data_dz_sigmadz->Fill(fdata_dz_sigmadz);
                        ++fdata_trkd0;
                        data_sigmapt_pt->Fill(fdata_sigmapt_pt);
                        ++fdata_trk;
                    }

                }
            }
        }
    }
    cout << "Before plotting." << endl;

    TFile data_plot("no_cut_data_histo.root", "recreate");
    TCanvas *canvas = new TCanvas ("data eta", "data #eta");
    //canvas->Divide (2,2);

    data_eta_histo->Scale(1/fdata_trketa);
    data_eta_histo->SetMaximum(0.04);
    data_eta_histo->Draw();
    data_eta_histo->Write();

    //canvas->cd(3);
    data_phi_histo->Scale(1/fdata_trk);
    data_phi_histo->SetMinimum(0.014);
    data_phi_histo->SetMaximum(0.024);
    data_phi_histo->Draw();
    data_phi_histo->Write();

    //canvas->cd(1);
    gPad->SetLogy();
    data_pt_histo->Scale(1/fdata_evt);
    data_pt_histo->Draw();
    data_pt_histo->Write();

    data_dz_sigmadz->Scale(1/fdata_trkd0);
    data_dz_sigmadz->Draw();
    data_dz_sigmadz->Write();
    data_dz_sigmadz->SetMinimum(1E-5);
    data_dz_sigmadz->SetMaximum(1E-1);

    data_d0_sigmad0->SetMinimum(1E-5);
    data_d0_sigmad0->SetMaximum(0.1);
    data_d0_sigmad0->Scale(1/fdata_trkd0);
    //data_d0_sigmad0->SetLineColor(kGreen+3);
    data_d0_sigmad0->Draw();
    data_d0_sigmad0->Write();

    data_sigmapt_pt->Scale(1/fdata_trkd0);
    data_sigmapt_pt->SetMinimum(1E-8);
    data_sigmapt_pt->Draw();
    data_sigmapt_pt->Write();

    //canvas->cd(2);


    //canvas->cd(4);

    data_plot.Write();

    //data_eta_histo->DrawNormalized("", 1);
    //data_eta_histo->Write();
    //data_eta_histo->Scale(1/fdata_trketa);
    //data_eta_histo->GetYaxis()->SetRange(0,0.04);
    //data_eta_histo->Draw();
    //data_eta_histo->Write();

    //canvas->cd(2);
    //gPad->SetLogy();
    //data_eta_histo->DrawNormalized("", 1);
    //data_eta_histo->Draw();

    /*canvas->cd(3);
    gPad->SetLogy();
    data_dz_sigmadz->Scale(1/fdata_trk);
    data_dz_sigmadz->Draw();

    canvas->cd(4);
    gPad->SetLogy();
    data_pt_histo->Scale(1/fdata_trk);
    data_pt_histo->Draw();*/

}
