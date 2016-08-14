#include "TROOT.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <algorithm>
#include "math.h"
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
//#include "TUnfold.h"
//#include "TUnfoldSys.h"
//#include "TUnfoldDensity.h"
//#include "TUnfoldBinning.h"
#include "TNamed.h"
//#include "TVectorDfwd.h"
//#include "RooUnfold.h"
//#include "RooUnfoldBayes.h"
//#include "RooUnfoldResponse.h"
//#include "RooUnfoldSvd.h"

using namespace std;
using namespace ROOT::Math;

void track_unfoldMC()
{
    //===========================implement cuts================================
    const double eta_cut = 2.;
    const double vtx_number_cut = 1.0;
    const float vtxxysize = 2;
    const float vtxzsize = 10;
    const double pt_cut = 0.5;
    const int lumi_cut = 90;
    const float ptErr_pt_cut = 0.05;
    const float dz_dzErr_cut = 3;
    //const float dz_dzErr_cut = 4;
    const float d0_d0Err_cut = 3;
    //const float d0_d0Err_cut = 4;
    const float dof_cut = 4;

    //===========================retrieve ROOT file============================

    TFile *datafile = TFile::Open("treesCUETP8M1_66.root", "READ");
    TTree *datatree = (TTree*)datafile->Get("UETree/data");

    //==============================================Histograms==============================================================

    TH1F *data_pt_histo = new TH1F ("data p_{T}", "Normalized data p_{T}", 200, 0, 10);
    TH1F *data_eta_histo = new TH1F("data eta", "Normalised data #eta", 50, -2.5, 2.5);
    TH1F *data_phi_histo = new TH1F ("data phi", "Normalized_data #phi", 60, -3, 3);
    TH1F *gen_pt_histo = new TH1F ("gen pT", "Normalized gen p_{T}", 200, 0, 10);
    TH1F *gen_eta_histo = new TH1F("gen eta", "Normalised gen #eta", 50, -2.5, 2.5);
    TH1F *gen_phi_histo = new TH1F ("gen phi", "Normalized gen #phi", 60, -3, 3);

    TH1F *data_d0leaf = new TH1F ("reco d0 leaf", "reco d_{0} leaf", 200, -10, 10);
    TH1F *data_sigmad0_wxwy = new TH1F ("reco sigmad0", "reco #sigma_{0}", 200, 0, 10);
    TH1F *data_d0_sigmad0 = new TH1F ("reco d0/sigmad0", "reco #frac{d_{0}}{#sigma_{0}}", 400, -20, 20);
    TH1F *data_d0_sigmad0_wxwy = new TH1F("reco d0/sigmad0 wxwy", "reco #frac{d_{0}}{#sigma_{0}} wxwy", 400, -20, 20);

    TH1F *data_dz_sigmadz = new TH1F ("reco dz/sigmadz", "reco #frac{d_{z}}{#sigma_{z}}", 400, -20, 20);
    TH1F *data_dz_sigmadz_vtxzerr = new TH1F ("reco dz/sigmadz vtxzerr", "reco #frac{d_{z}}{#sigma_{z}} vtxzerr", 400, -20, 20);
    TH1F *data_dz_sigmadz_wz = new TH1F ("reco dz/sigmadz wz", "reco #frac{d_{z}}{#sigma_{z}}", -200, -10, 10);

    TH1F *data_sigmapt_pt = new TH1F ("reco sigma_pt/pt", "reco #frac{#sigma_{p_{T}}}{p_{T}}", 20, 0, 0.2);

    TH1F *reco_multiplicity = new TH1F ("reco_multiplicity", "Normalized_reco_Multiplicity", 50, 0, 50);
    TH1F *gen_multiplicity = new TH1F ("gen_multiplicity", "Normalized_gen_Multiplicity", 50, 0, 50);

    //==============================================Variables==============================================================

    float fdata_evt = 0;
    float fdata_trk = 0;
    float fdata_trketa = 0;
    float fdata_trkd0 = 0;
    float fdata_trkpt = 0;
    float fdata_trkdz = 0;
    float fdata_trkdpt = 0;
    float fdata_dz, fdata_sigmadz, fdata_dz_sigmadz, fdata_dz_sigmadz_vtxzerr, fdata_dz_sigmadz_wz;
    float fdata_d0, fdata_sigmad0, fdata_sigmad0_wxwy, fdata_d0_sigmad0, fdata_d0_sigmad0calc, fdata_d0_sigmad0_wxwy;
    float fdata_sigmapt_pt, fdata_trkdx, fdata_trkdy;
    float fdata_wx, fdata_wy,fdata_wz, fdata_vtxxysize, fdata_vtxzsize, fdata_sigmad0calc;
    float fdata_numberofvtxxBS, fdata_vtxxBSvalue, fdata_vtxxBSlower, fdata_vtxxBSupper;
    float fdata_numberofvtxyBS, fdata_vtxyBSvalue, fdata_vtxyBSlower, fdata_vtxyBSupper;
    float fdata_numberofvtxzBS, fdata_vtxzBSvalue, fdata_vtxzBSlower, fdata_vtxzBSupper;
    float fdata_vtxzminusvtxz;
    float fdata_sqvtxx = 0;
    float fdata_sqvtxxnumber = 0;
    float fdata_sqvtxy = 0;
    float fdata_sqvtxynumber = 0;
    float fdata_sqvtxz = 0;
    float fdata_sqvtxznumber = 0;
    float fdata_numberoftrkdx = 0;
    float fdata_numberoftrkdy = 0;
    float fdata_numselectedvtxz = 0;
    float fdata_numvtxzminusvtxz = 0;
    float freco_multiplicity = 0;
    float fgen_multiplicity = 0;
    int ndata_numberofvtxx, ndata_numberofvtxy, ndata_numberofvtxz, ndata_numberofvtxxBS, ndata_numberofvtxyBS, ndata_numberofvtxzBS;
    int ngen_numberofvtxx, ngen_numberofvtxy, ngen_numberofvtxz;
    float freco_unfold, fMC_unfold;
    //float fdata_dz_sigmadz = 0;
    //float fdata_d0_sigmad0 = 0;
    //float fdata_sigmapt_pt = 0;
    int ndata_totaltrk = 0;
    int ngen_totaltrk = 0;
    bool is_gen, is_reco;

        //===========================define variables to read TTree================


    vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *reco_tracks = 0;
    datatree->SetBranchAddress("recoTracksp4", &reco_tracks);

    int ndata_lumi;
    datatree->SetBranchAddress("lumi", &ndata_lumi);

    int ndata_zerobias;
    datatree->SetBranchAddress("trgZeroBias", &ndata_zerobias);

    vector<float> *fvecdata_vtxx = 0;
    datatree->SetBranchAddress("vtxx", &fvecdata_vtxx);

    vector<float> *fvecdata_vtxy = 0;
    datatree->SetBranchAddress("vtxy", &fvecdata_vtxy);

    vector<float> *fvecdata_vtxz = 0;
    datatree->SetBranchAddress("vtxz", &fvecdata_vtxz);

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

    vector<float> *fvecdata_d0 = 0;
    datatree->SetBranchAddress("recoTracksd0", &fvecdata_d0);

    vector<float> *fvecdata_d0err = 0;
    datatree->SetBranchAddress("recoTracksd0Err", &fvecdata_d0err);

    vector<float> *fvecdata_dz = 0;
    datatree->SetBranchAddress("recoTracksdz", &fvecdata_dz);

    vector<float> *fvecdata_dzerr = 0;
    datatree->SetBranchAddress("recoTracksdzErr", &fvecdata_dzerr);

    vector<float> *fvecdata_pterr = 0;
    datatree->SetBranchAddress ("recoTracksptErr", &fvecdata_pterr);

//================================================= Read Gen Branches ================================================================================

    vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *gen_tracks = 0;
    datatree->SetBranchAddress("genParticlesp4", &gen_tracks);

    vector<float> *fvecgen_simvtxx = 0;
    datatree->SetBranchAddress("simvtxx", &fvecgen_simvtxx);

    vector<float> *fvecgen_simvtxy = 0;
    datatree->SetBranchAddress("simvtxy", &fvecgen_simvtxy;

    vector<float> *fvecgen_simvtxz = 0;
    datatree->SetBranchAddress("simvtxz", &fvecgen_simvtxz);

//================================================= Start of RMS Loop ================================================================================

    Int_t ndata_totalEvt = (Int_t)datatree->GetEntries();
    cout << "There is a total of " << ndata_totalEvt << " events." << endl;

    for (int ev = 0; ev < ndata_totalEvt; ++ev)
    {
        datatree->GetEntry(ev);

        if (ndata_lumi >= lumi_cut)
        {
            if (ndata_zerobias == 1)
            {
                ndata_numberofvtxx = fvecdata_vtxx->size();
                ndata_numberofvtxy = fvecdata_vtxy->size();
                ndata_numberofvtxz = fvecdata_vtxz->size();
                ndata_numberofvtxxBS = fvecdata_vtxxBS->size();
                ndata_numberofvtxyBS = fvecdata_vtxyBS->size();
                ndata_numberofvtxzBS = fvecdata_vtxzBS->size();

                fdata_vtxxBSlower = *min_element(fvecdata_vtxxBS->begin(), fvecdata_vtxxBS->end());
                fdata_vtxxBSupper = *max_element(fvecdata_vtxxBS->begin(), fvecdata_vtxxBS->end());
                fdata_vtxyBSlower = *min_element(fvecdata_vtxyBS->begin(), fvecdata_vtxyBS->end());
                fdata_vtxyBSupper = *max_element(fvecdata_vtxyBS->begin(), fvecdata_vtxyBS->end());
                fdata_vtxzBSlower = *min_element(fvecdata_vtxzBS->begin(), fvecdata_vtxzBS->end());
                fdata_vtxzBSupper = *max_element(fvecdata_vtxzBS->begin(), fvecdata_vtxzBS->end());

                for(int vtxx = 0; vtxx != ndata_numberofvtxx; ++vtxx)
                {
                    //fdata_sqvtxx += pow(((*fvecdata_vtxx)[vtxx]) - ((*fvecdata_vtxxBS)[vtxx]), 2);
                    fdata_sqvtxx += pow((((*fvecdata_vtxx)[vtxx]) - (fdata_vtxxBSlower)), 2);
                    ++fdata_sqvtxxnumber;
                }

                for(int vtxy = 0; vtxy != ndata_numberofvtxy; ++vtxy)
                {
                    //fdata_sqvtxy += pow(((*fvecdata_vtxy)[vtxy]) - ((*fvecdata_vtxyBS)[vtxy]), 2);
                    fdata_sqvtxy += pow((((*fvecdata_vtxy)[vtxy]) - (fdata_vtxyBSupper)), 2);
                    ++fdata_sqvtxynumber;
                }

                for(int vtxz = 0; vtxz != ndata_numberofvtxz; ++vtxz)
                {
                    //fdata_sqvtxz += pow(((*fvecdata_vtxz)[vtxz]) - ((*fvecdata_vtxzBS)[vtxz]), 2);
                    fdata_sqvtxz += pow((((*fvecdata_vtxz)[vtxz]) - (fdata_vtxzBSupper)), 2);
                    ++fdata_sqvtxznumber;
                }

                for (int vtxzBS = 0; vtxzBS != ndata_numberofvtxzBS; ++vtxzBS)
                {
                    if ((*fvecdata_vtxzBS)[vtxzBS] == fdata_vtxzBSupper)
                    {
                        ++fdata_numberofvtxzBS;
                    }

                }
            }
        }
    }
    cout << "Largest BS z-coordinate is " << fdata_vtxzBSupper << endl;
    fdata_wx = sqrt((fdata_sqvtxx) / (fdata_sqvtxxnumber)); //RMS(?) of position of vtx x-coordinate. Averaged over number of vtx
    fdata_wy = sqrt((fdata_sqvtxy) / (fdata_sqvtxynumber)); //RMS(?) of position of vtx y-coordinate Averaged over number of vtx
    fdata_wz = sqrt((fdata_sqvtxz) / (fdata_sqvtxznumber)); //

//========================================================= End of RMS Loop ============================================================

    //========================================================= Start of Event Loop ============================================================

    for (Int_t i = 0; i < ndata_totalEvt; ++i)
    {
        datatree->GetEntry(i);
        //cout << "At entry " << i << endl;
        //based on http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html#use arguments are the dimension of the reco and gen multiplicity histo
        //RooUnfoldResponse response (hist_measured, hist_truth)
        RooUnfoldResponse response (1, 1);

        if (ndata_lumi >= lumi_cut)
        {
            if (ndata_zerobias == 1)
            {
                ndata_totaltrk = reco_tracks->size();
                ngen_totaltrk = gen_tracks->size();
                ++fdata_evt;

                int vtxdof = 0;
                ndata_numberofvtxx = fvecdata_vtxx->size();
                ndata_numberofvtxy = fvecdata_vtxy->size();
                ndata_numberofvtxz = fvecdata_vtxz->size();
                ndata_numberofvtxxBS = fvecdata_vtxxBS->size();
                ndata_numberofvtxyBS = fvecdata_vtxyBS->size();
                ndata_numberofvtxzBS = fvecdata_vtxzBS->size();

                ngen_numberofvtxx = fvecgen_simvtxx->size();
                ngen_numberofvtxy = fvecgen_simvtxy->size();
                ngen_numberofvtxz = fvecgen_simvtxz->size();

//========================================================= Start of Gen Vertex Loop ============================================================

                is_gen = false;
                for (int gen_vtxnumber = 0; gen_vtxnumber != ngen_numberofvtxx; ++gen_vyxnumber)
                {
                    if (ngen_numberofvtxx == vtx_number_cut && ngen_numberofvtxy == vtx_number_cut && ngen_numberofvtxz == vtx_number_cut)
                    {
                        is_gen = true;
                        fgen_multiplicity = 0;
//========================================================= Start of Gen Track Loop ============================================================

                        for (int gen_trk = 0; gen_trk != ngen_totaltrk; ++gen_trk)
                        {
                            XYZTVector gen_vec = (*gen_tracks)[t];

                            if(fabs(gen_vec.Eta()) < eta_cut && fabs(gen_vec.Pt()) >= pt_cut)
                            {
                                ++fgen_multiplicity;
                                gen_phi_histo->Fill(MC_vec.Phi());
                            }

                            if(fabs(MC_vec.Eta()) <= eta_cut)
                            {
                                gen_pt_histo->Fill(MC_vec.Pt());
                            }

                            if(fabs(MC_vec.Pt()) >= pt_cut)
                            {
                                gen_eta_histo->Fill(MC_vec.Eta());
                            }
                        }

//========================================================= End of Gen Vertex Loop ============================================================
                    gen_multiplicity->Fill(fgen_multiplicity);
                    }
                }
//========================================================= End of Gen Vertex Loop ============================================================


//========================================================= Start of Reco Vertex Loop ============================================================

                is_reco = false;

                for (int reco_vtxnumber = 0; reco_vtxnumber != ndata_numberofvtxx; ++reco_vtxnumber)
                {

                    if((*nvecdata_vtxndof)[vtxdof] > dof_cut)
                    {
                        fdata_vtxxysize = sqrt(pow(((*fvecdata_vtxx)[reco_vtxnumber]) - ((*fvecdata_vtxxBS)[reco_vtxnumber]), 2) + pow(((*fvecdata_vtxy)[reco_vtxnumber]) - ((*fvecdata_vtxyBS)[reco_vtxnumber]), 2));
                        fdata_vtxzsize = fabs(((*fvecdata_vtxz)[reco_vtxnumber]) - ((*fvecdata_vtxzBS)[reco_vtxnumber]));

                        if (ndata_numberofvtxx == vtx_number_cut && ndata_numberofvtxy == vtx_number_cut && ndata_numberofvtxz == vtx_number_cut)
                        {
                            if (fdata_vtxxysize < vtxxysize && fdata_vtxzsize < vtxzsize)
                            {
                                freco_multiplicity = 0;
                                is_reco = true;

//========================================================= Start of Reco Track Loop ============================================================

                                for (int t = 0; t != ndata_totaltrk; ++t)
                                {
                                    XYZTVector data_vec = (*reco_tracks)[t];
                                    //cout << "Within track " << endl;

//========================================================= dz variales ============================================================

                                    //using formula from paper of run 1 result, dz is leaf value
                                    fdata_dz = (*fvecdata_dz)[t];
                                    fdata_sigmadz = sqrt(pow(((*fvecdata_dzerr)[t]), 2) + pow(fdata_wz, 2));
                                    fdata_dz_sigmadz = ((*fvecdata_dz)[t])/((*fvecdata_dzerr)[t]);
                                    fdata_dz_sigmadz_wz = ((fdata_dz) / (fdata_sigmadz));

                                    //fdata_dz_sigmadz = ((*fvecdata_dz)[t])/((*fvecdata_vtxzerr)[t]);
                                    //fdata_dz_sigmadz_vtxzerr = (((*fvecdata_dz)[t])/sqrt(pow(((*fvecdata_dzerr)[t]),2)+pow((*fvecdata_vtxzerr)[t],2)));

//========================================================= d0 variales ============================================================

                                    //using formula from paper of run 1 result, d0 is leaf value
                                    fdata_d0_sigmad0 = ((*fvecdata_d0)[t])/((*fvecdata_d0err)[t]);
                                    fdata_d0 = (*fvecdata_d0)[t];
                                    fdata_sigmad0_wxwy = sqrt(pow(((*fvecdata_d0err)[t]), 2) + (fdata_wx)*(fdata_wy));
                                    //fdata_sigmad0 = fmin(0.05,(fdata_sigmad0_wxwy));
                                    data_d0leaf->Fill(fdata_d0);
                                    fdata_d0_sigmad0_wxwy = ((fdata_d0) / (fdata_sigmad0_wxwy));
                                    //2 differeent versions of d0/sigma_d0
                                    //fdata_d0_sigmad0calc = (((*fvecdata_d0)[t])/sqrt(pow(((*fvecdata_d0err)[t]),2)+pow((*fvecdata_vtxxerr)[t],2)+pow((*fvecdata_vtxyerr)[t],2)));

//========================================================= pT variales ============================================================

                                    fdata_sigmapt_pt = (((*fvecdata_pterr)[t])/(data_vec.Pt()));

//========================================================= Cuts on Tracks ============================================================

                                    if ((*nvecdata_highpurity)[t] == 1)
                                    {

                                        if (fabs(data_vec.Eta()) <= eta_cut)
                                        {
                                            data_pt_histo->Fill(data_vec.Pt());
                                            ++fdata_trkpt;
                                        }

                                        //if (data_vec.Pt() >= pt_cut && fabs(fdata_sigmapt_pt) < ptErr_pt_cut && fabs(fdata_d0_sigmad0run1) < d0_d0Err_cut && fabs(fdata_dz_sigmadzrun1 < dz_dzErr_cut))
                                        if (data_vec.Pt() >= pt_cut && fabs(fdata_sigmapt_pt) < ptErr_pt_cut && fabs(fdata_dz_sigmadz_wz < dz_dzErr_cut))
                                        {
                                            data_eta_histo->Fill(data_vec.Eta());
                                            ++fdata_trketa;
                                        }

                                        if (data_vec.Pt() >= pt_cut && fabs(data_vec.Eta()) <= eta_cut)
                                        {
                                            data_phi_histo->Fill(data_vec.Phi());
                                            ++fdata_trk;

                                            //if (fabs(fdata_d0_sigmad0calc) < d0_d0Err_cut && fabs(fdata_sigmapt_pt) < ptErr_pt_cut && fdata_d0 < fmin(0.05, fdata_sigmad0calc)))
                                            if (fabs(fdata_d0_sigmad0calc) < d0_d0Err_cut && fabs(fdata_sigmapt_pt) < ptErr_pt_cut )
                                            {
                                                data_dz_sigmadz->Fill(fdata_dz_sigmadz);
                                                data_dz_sigmadz_vtxzerr->Fill(fdata_dz_sigmadz_vtxzerr);
                                                data_dz_sigmadz_wz->Fill(fdata_dz_sigmadz_wz);
                                                ++fdata_trkdz;
                                            }

                                            //if (fabs(fdata_dz_sigmadz_wz) < dz_dzErr_cut && fabs(fdata_sigmapt_pt) < ptErr_pt_cut)
                                            if (fabs(fdata_dz_sigmadz) < dz_dzErr_cut && fabs(fdata_sigmapt_pt) < ptErr_pt_cut)
                                            {
                                                data_d0_sigmad0->Fill(fdata_d0_sigmad0);
                                                data_d0_sigmad0_wxwy->Fill(fdata_d0_sigmad0_wxwy);
                                                //fdata_sigmad0 = sqrt(pow(((*fvecdata_d0)[t]), 2) + (fdata_wx)*(fdata_wy)); //check the formula again
                                                //fdata_sigmad0 = sqrt(pow((*fvecdata_d0)[t], 2) + 0.01429*0.01727);
                                                //data_d0leaf->Fill(fdata_d0);
                                                data_sigmad0_wxwy->Fill(fdata_sigmad0_wxwy);
                                                ++fdata_trkd0;
                                                ++freco_multiplicity;
                                            }

                                            if (fabs(fdata_dz_sigmadz) < dz_dzErr_cut && fabs(fdata_d0_sigmad0calc) < d0_d0Err_cut)
                                            {
                                                data_sigmapt_pt->Fill(fdata_sigmapt_pt);
                                                ++fdata_trkdpt;
                                            }

                                        }
                                    }
                                }
                                reco_multiplicity->Fill(freco_multiplicity);
//========================================================= End of Track Loop ============================================================

                            }
                        }
                    }
                }

//========================================================= Train ========================================================================

                if (is_reco == true && is_gen == true)
                {
                    response.Fill (reco_multiplicity, gen_multiplicity);
                }

                if (is_reco == true && is_gen == false)
                {
                    response.Fake (reco_multiplicity);
                }

                if (is_reco == false && is_gen == true)
                {
                    response.Miss(gen_multiplicity);
                }

//========================================================= End of Vertex Loop ============================================================
            }

//========================================================= End of ZeroBias ================================================================

        }
    }
    cout << "Before plotting." << endl;
    cout << "wx is " << fdata_wx << endl;
    cout << "wy is " << fdata_wy << endl;
    cout << "wz is " << fdata_wz << endl;
    cout << "cut on d0 is " << sqrt(0.0025 - fdata_wx*fdata_wy) << endl;
    cout << "Largest z-coordinate of BS is " << fdata_vtxzBSupper << endl;

    TFile data_plot("no_cut_reco_histo.root", "recreate");
    TCanvas *canvas = new TCanvas ("data eta", "data #eta");
    //canvas->Divide (2,2);

//=================================================Output 3 Main Variables==============================================================

    data_eta_histo->Scale(1/fdata_trketa);
    data_eta_histo->SetMaximum(0.04);
    data_eta_histo->SetLineColor(kGreen);
    data_eta_histo->SetMarkerStyle(kFullDotMedium);
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
    data_pt_histo->Scale(1/fdata_trkpt);
    data_pt_histo->Draw();
    data_pt_histo->Write();

    gen_pt_histo->Scale(1/ngen_totaltrk);
    gen_pt_histo->Draw();
    gen_pt_histo->Write();

    gen_eta_histo->Scale(1/ngen_totaltrk);
    gen_eta_histo->Draw();
    gen_eta_histo->Write();

    gen_phi_histo->Scale(1/ngen_totaltrk);
    gen_phi_histo->SetMinimum(0.014);
    gen_phi_histo->SetMaximum(0.024);
    gen_phi_histo->Draw();
    gen_phi_histo->Write();

//=================================================Output dz/sigma_dz==============================================================

    data_dz_sigmadz->Scale(1/fdata_trkdz);
    data_dz_sigmadz->SetMinimum(1E-5);
    data_dz_sigmadz->SetMaximum(1E-1);
    data_dz_sigmadz->Draw();
    data_dz_sigmadz->Write();

    data_dz_sigmadz_vtxzerr->Scale(1/fdata_trkdz);
    data_dz_sigmadz_vtxzerr->Draw();
    data_dz_sigmadz_vtxzerr->Write();

    data_dz_sigmadz_wz->Scale(1/fdata_trkdz);
    data_dz_sigmadz_wz->Write();

//=================================================Output d0/sigma_d0==============================================================

    data_d0_sigmad0->SetMinimum(1E-5);
    data_d0_sigmad0->SetMaximum(0.1);
    data_d0_sigmad0->Scale(1/fdata_trkd0);
    //data_d0_sigmad0->SetLineColor(kGreen+3);
    data_d0_sigmad0->Draw();
    data_d0_sigmad0->Write();

    data_d0leaf->Scale(1/fdata_trkd0);
    data_d0leaf->Write();
    data_sigmad0_wxwy->Scale(1/fdata_trkd0);
    data_sigmad0_wxwy->Write();

    data_d0_sigmad0_wxwy->Scale(1/fdata_trkd0);
    data_d0_sigmad0_wxwy->Draw();
    data_d0_sigmad0_wxwy->Write();

//=================================================Output sigma_pT/pT==============================================================

    data_sigmapt_pt->Scale(1/fdata_trkdpt);
    data_sigmapt_pt->SetMinimum(1E-8);
    data_sigmapt_pt->Draw();
    data_sigmapt_pt->Write();

//=================================================Output Multiplicity==============================================================

    reco_multiplicity->Scale(1/ndata_totaltrk);
    reco_multiplicity->Draw();
    //reco_multiplicity->DrawNormalized("", 1);
    reco_multiplicity->Write();

    gen_multiplicity->DrawNormalized("", 1);
    gen_multiplicity->Write();

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

    TFile unfold_histo("unfold.root", "new");
    RooUnfoldBayes unfold (&response, reco_multiplicity, 4);

    TH1F *hReco = (TH1F*) unfold.Hreco();
    hReco->Draw();
    hReco->Write();
    reco_multiplicity->SetLineColor(kGreen);
    reco_multiplicity->Draw("SAME");
    gen_multiplicity->SetLineColor(8);
    gen_multiplicity->Draw("SAME");

    unfold_histo.Write();


}
