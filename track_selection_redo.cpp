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

using namespace std;
using namespace ROOT::Math;

void track_selection_redo()
{
//===========================implement cuts================================
    const double eta_cut = 2.;
    const double vtx_number_cut = 1.0;
    const float vtxxysize = 2;
    const float vtxzsize = 10;
    const double pt_cut = 0.5;
    const int lumi_cut = 90;
    const float ptErr_pt_cut = 0.05;
    //const float dz_dzErr_cut = 3;
    const float dz_dzErr_cut = 4;
    //const float d0_d0Err_cut = 3;
    const float d0_d0Err_cut = 4;
    const float dof_cut = 4;
    float fdata_multiplicity_norm = 0;
    float freco_multiplicity_norm = 0;

    TH1F *data_multiplicity = new TH1F("Normalized_Data_Multiplicity", "Normalized Data Multiplicity", 200, 0, 200);
    TH1F *reco_multiplicity = new TH1F("Normalized_Reco_Multiplicity", "Normalized Reco Multiplicity", 200, 0, 200);

    int nselect = 0;
    cout << "Enter 0 for data, 1 for reco and 2 for both: ";
    cin >> nselect;

    while (nselect > 2)
    {
        cout << "Please enter 0, 1 or 2 only: ";
        cin >> nselect;
    }
//==================================================== Data Loop ===================================================================

    if ((nselect == 0) || (nselect == 2))
    {
//==============================================Variables==============================================================

        float fdata_evt = 0;
        float fdata_trkphi = 0;
        float fdata_trketa = 0;
        float fdata_trkd0 = 0;
        float fdata_trkpt = 0;
        float fdata_trkdz = 0;
        float fdata_trkdpt = 0;
        float fdata_dz_sigmadz, fdata_d0_sigmad0, fdata_sigmapt_pt, fdata_dz_sigmadzcalc, fdata_d0_sigmad0calc, fdata_sigmad0, fdata_sigmad0run1, fdata_d0, fdata_d0_sigmad0run1, fdata_trkdx, fdata_trkdy;
        float fdata_wx, fdata_wy,fdata_wz, fdata_vtxxysize, fdata_vtxzsize, fdata_sigmad0calc, fdata_dz, fdata_sigmadz, fdata_dz_sigmadzrun1;
        float fdata_numberofvtxxBS, fdata_vtxxBSvalue, fdata_vtxxBSlower, fdata_vtxxBSupper;
        float fdata_numberofvtxyBS, fdata_vtxyBSvalue, fdata_vtxyBSlower, fdata_vtxyBSupper;
        float fdata_numberofvtxzBS, fdata_vtxzBSvalue, fdata_vtxzBSlower, fdata_vtxzBSupper;
        float fdata_vtxzminusvtxz, fdata_multiplicity;
        //float fdata_multiplicity_norm = 0;
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
        float fdata_trkvalidhits = 0;
        float fdata_trkchi2n = 0;
        int ndata_numberofvtxx, ndata_numberofvtxy, ndata_numberofvtxz, ndata_numberofvtxxBS, ndata_numberofvtxyBS, ndata_numberofvtxzBS;
        //float fdata_dz_sigmadz = 0;
        //float fdata_d0_sigmad0 = 0;
        //float fdata_sigmapt_pt = 0;
        int ndata_totaltrk;

//===========================retrieve ROOT file============================

        TFile *datafile = TFile::Open("tree1.root", "READ");
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

        vector<float> *fvecdata_vtxxerr = 0;
        datatree->SetBranchAddress("vtxxErr", &fvecdata_vtxxerr);

        vector<float> *fvecdata_vtxyerr = 0;
        datatree->SetBranchAddress("vtxyErr", &fvecdata_vtxyerr);

        vector<float> *fvecdata_vtxzerr = 0;
        datatree->SetBranchAddress("vtxzErr", &fvecdata_vtxzerr);

        vector<int> *nvecdata_validhits = 0;
        datatree->SetBranchAddress("recoTracksnValidHits", &nvecdata_validhits);

        vector<float> *fvecdata_trackschi2n = 0;
        datatree->SetBranchAddress("recoTrackschi2n", &fvecdata_trackschi2n);


    //============================================== Histos for pT, eta, phi ==============================================================

        TH1F *data_pt_histo = new TH1F ("data pT", "Normalized data p_{T}", 200, 0, 10);
        TH1F *data_eta_histo = new TH1F("data eta", "Normalised data #eta", 50, -2.5, 2.5);
        TH1F *data_phi_histo = new TH1F ("data phi", "Normalized_data #phi", 60, -3, 3);

    //================================================== Histos for dz and sigma_dz =========================================================================
        TH1F *data_dzleaf = new TH1F ("data_dz", "data dz", 400, -10, 10);
        TH1F *data_sigmadzcalc = new TH1F ("data_sigmadzcalc", "data #sigma_{z} calc", 200, 0, 4);//plot of sigmadz using formula from WY
        TH1F *data_sigmadzrun1 = new TH1F ("data_sigmadzrun1", "data #sigma_{z} run 1", 200, 0, 4);//plot of sigmadz using formula from run 1
        TH1F *data_dz_sigmadz = new TH1F ("data_dz_sigmadz", "data d_{z}/#sigma_{z}", 160, -20, 20);
        TH1F *data_dz_sigmadzcalc = new TH1F ("data_dz_sigmadzcalc", "data d_{z}/#sigma_{z} calc", 160, -20, 20);
        TH1F *data_dz_sigmadzrun1 = new TH1F ("data_dz_sigmadz run 1", "data d_{z}/#sigma_{z} run 1", 160, -20, 20);
        //TH1F *data_dz_sigmadzcalcb4cut = new TH1F ("data_dz_sigmadzcalcb4cut", "data d_{z}/#sigma_{z} calc", 160, -20, 20);

    //================================================== Histos for d0 and sigma_d0 =========================================================================

        TH1F *data_d0leaf = new TH1F ("data_d0_calcrun1", "data d_{0} calc run 1", 100, -1, 1); //plot of d0 using leaf value to determine how much data is cut
        TH1F *data_sigmad0calc = new TH1F ("data_sigmad0calc", "data #sigma_{xy} calc", 200, 0, 4);//plot of sigmad0 using formula from WY
        TH1F *data_sigmad0run1 = new TH1F ("data_sigmad0run1", "data #sigma_{xy} run 1", 200, 0, 4);//plot of sigmad0 using run 1 formula

        TH1F *data_d0_sigmad0 = new TH1F ("data_d0_sigmad0", "data d_{0}/#sigma_{xy}", 160, -20, 20); //both leaf values
        TH1F *data_d0_sigmad0run1 = new TH1F ("data_d0_sigmad0calcrun1", "data d_{0}/#sigma_{xy} calc run 1", 160, -20, 20); //plot using run 1 formula
        TH1F *data_d0_sigmad0calc = new TH1F ("data_d0_sigmad0calc", "data d_{0}/#sigma_{xy} calc", 160, -20, 20); //using formula from WY
        //TH1F *data_d0_sigmad0calcb4cut = new TH1F ("data_d0_sigmad0calcb4cut", "data d_{0}/#sigma_{xy} calc", 160, -20, 20);
        //TH1F *data_d0_sigmad0calc = new TH1F ("data_d0_sigmad0calc", "data d_{0}/#sigma_{xy} calc", 150, 0, 0.5);

    //================================================== Histos for pT and sigma_pT =========================================================================

        TH1F *data_sigmapt_pt = new TH1F ("data_sigmapt_pt", "data #sigma_{p_{T}}/p_{T}", 20, 0, 0.2);

    //================================================== Histos for chi2n and ValidHits =========================================================================

        TH1F *data_validhits = new TH1F ("Tracks_vs_validhits", "Tracks vs validhits", 50, 0, 50);
        TH1F *data_chi2n = new TH1F ("Tracks_vs_chi2n", "Tracks vs #chi^{2/ndof}", 50, 0, 5);

    //================================================== Histos for Multiplicity =========================================================================

        TH1F *data_multiplicity = new TH1F("Normalized_Multiplicity", "Normalized Multiplicity", 200, 0, 200);

        TH1F *data_vtxzminusvtxz = new TH1F ("vtxzminusvtxz", "vtxzminusvtxz", 800, -4, 4);
        TH1F *data_vtxzposn = new TH1F ("vtxzposn", "vtxzpos^{n}", 400, -20, 20);


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

    //========================================================= Start of Evt Loop ================================================================

        for (Int_t i = 0; i < ndata_totalEvt; ++i)
        {
            datatree->GetEntry(i);
            //cout << "At entry " << i << endl;

            if (ndata_lumi >= lumi_cut)
            {
                if (ndata_zerobias == 1)
                {
                    ndata_totaltrk = data_tracks->size();

                    int vtxdof = 0;

                    ndata_numberofvtxxBS = fvecdata_vtxxBS->size();
                    ndata_numberofvtxyBS = fvecdata_vtxyBS->size();
                    ndata_numberofvtxzBS = fvecdata_vtxzBS->size();
                    fdata_multiplicity = 0;
                    ++fdata_evt;

    //========================================================= Start of Vertex Loop ================================================================

                    for (int vtxnumber = 0; vtxnumber != ndata_numberofvtxxBS; ++vtxnumber)
                    {
                        if((*nvecdata_vtxndof)[vtxdof] > dof_cut)
                        {
                            fdata_vtxxysize = sqrt(pow(fabs(((*fvecdata_vtxx)[vtxnumber]) - ((*fvecdata_vtxxBS)[vtxnumber])), 2) + pow(fabs((*fvecdata_vtxy)[vtxnumber] - (*fvecdata_vtxyBS)[vtxnumber]), 2));
                            fdata_vtxzsize = fabs(((*fvecdata_vtxz)[vtxnumber]) - ((*fvecdata_vtxzBS)[vtxnumber]));

                            //if (fdata_vtxxysize <= vtxxysize && fdata_vtxzsize <= vtxzsize)
                            //{
                            if (ndata_numberofvtxzBS == vtx_number_cut)
                            {
                                data_vtxzposn->Fill((*fvecdata_vtxzBS)[vtxnumber]);
                                ++fdata_numselectedvtxz;

    //========================================================= Start of Trk Loop ================================================================

                                for (int t = 0; t != ndata_totaltrk; ++t)
                                {
                                    XYZTVector data_vec = (*data_tracks)[t];
                                    //cout << "Within track " << endl;
                                    //using formula from paper of run 1 result, dz is leaf value

                        //======================================= dz ===================================================

                                    fdata_dz = (*fvecdata_dz)[t];
                                    data_dzleaf->Fill((*fvecdata_dz)[t]);
                                    fdata_sigmadz = sqrt(pow(((*fvecdata_dzerr)[t]), 2) + pow(fdata_wz, 2));

                                    fdata_dz_sigmadz = ((*fvecdata_dz)[t])/((*fvecdata_dzerr)[t]); //leaf
                                    fdata_dz_sigmadzcalc = ((fdata_dz)/sqrt(pow(((*fvecdata_dzerr)[t]),2)+pow((*fvecdata_vtxzerr)[t],2))); //WY
                                    fdata_dz_sigmadzrun1 = ((fdata_dz) / (fdata_sigmadz)); //run 1 formula

                        //======================================= d0 ===================================================

                                    fdata_d0 = (*fvecdata_d0)[t];
                                    fdata_sigmad0run1 = sqrt(pow(((*fvecdata_d0err)[t]), 2) + (fdata_wx)*(fdata_wy));
                                    fdata_sigmad0calc = sqrt(pow(((*fvecdata_d0err)[t]),2)+pow((*fvecdata_vtxxerr)[t],2)+pow((*fvecdata_vtxyerr)[t],2)); //WY

                                    fdata_d0_sigmad0 = ((*fvecdata_d0)[t])/((*fvecdata_d0err)[t]); //both leaf values
                                    fdata_d0_sigmad0run1 = (((*fvecdata_d0)[t]) / (fdata_sigmad0run1)); //run 1 formula
                                    fdata_d0_sigmad0calc = (((*fvecdata_d0)[t]) / (fdata_sigmad0calc)); //with sigmad0 from WY

                                    //fdata_sigmad0 = fmin(0.05,(fdata_sigmad0run1));

                                    //fdata_d0_sigmad0 = ((*fvecdata_d0)[t])/sqrt((pow((*fvecdata_vtxxerr)[t],2)+pow((*fvecdata_vtxyerr)[t],2)));
                                    //fdata_d0_sigmad0run1 = ((fdata_d0) / (fdata_sigmad0));

                        //======================================= pT ===================================================

                                    fdata_sigmapt_pt = (((*fvecdata_pterr)[t])/(data_vec.Pt())); //both leaf values

                        //======================================= Trk cuts ===================================================

                                    if ((*nvecdata_highpurity)[t] == 1)
                                    {
                                        if (fabs(data_vec.Eta()) <= eta_cut)
                                        {
                                            data_pt_histo->Fill(data_vec.Pt());
                                            ++fdata_trkpt;
                                        }

                                        //if (data_vec.Pt() >= pt_cut && fabs(fdata_sigmapt_pt) < ptErr_pt_cut && fabs(fdata_dz_sigmadzrun1 < dz_dzErr_cut))
                                        //if ((data_vec.Pt() >= pt_cut) && (fabs(fdata_sigmapt_pt) < ptErr_pt_cut))
                                        if (data_vec.Pt() >= pt_cut)
                                        {
                                            data_eta_histo->Fill(data_vec.Eta());
                                            ++fdata_trketa;
                                        }

                                        if ((data_vec.Pt() >= pt_cut) && (fabs(data_vec.Eta()) <= eta_cut))
                                        {
                                            data_phi_histo->Fill(data_vec.Phi());
                                            ++fdata_trkphi;

                                            data_sigmapt_pt->Fill(fdata_sigmapt_pt);
                                            ++fdata_trkdpt;

                                            data_validhits->Fill((*nvecdata_validhits)[t]);
                                            ++fdata_trkvalidhits;

                                            data_chi2n->Fill((*fvecdata_trackschi2n)[t]);
                                            ++fdata_trkchi2n;

                                            //if (fabs(fdata_d0_sigmad0calc) < d0_d0Err_cut && fabs(fdata_sigmapt_pt) < ptErr_pt_cut && fdata_d0 < fmin(0.05, fdata_sigmad0calc)))
                                            //if (fabs(fdata_d0_sigmad0calc) < d0_d0Err_cut && fabs(fdata_sigmapt_pt) < ptErr_pt_cut)
                                            if (fabs(fdata_sigmapt_pt) < ptErr_pt_cut)
                                            {
                                                data_dz_sigmadz->Fill(fdata_dz_sigmadz);
                                                data_dz_sigmadzcalc->Fill(fdata_dz_sigmadzcalc);
                                                data_dz_sigmadzrun1->Fill(fdata_dz_sigmadzrun1);
                                                ++fdata_trkdz;
                                            }

                                            //if (fabs(fdata_dz_sigmadzrun1) < dz_dzErr_cut && fabs(fdata_sigmapt_pt) < ptErr_pt_cut)
                                            if (fabs(fdata_sigmapt_pt) < ptErr_pt_cut)
                                            {
                                                data_d0leaf->Fill((*fvecdata_d0)[t]);
                                                data_sigmad0run1->Fill(fdata_sigmad0run1);
                                                data_d0_sigmad0->Fill(fdata_d0_sigmad0);
                                                data_d0_sigmad0run1->Fill(fdata_d0_sigmad0run1);
                                                data_d0_sigmad0calc->Fill(fdata_d0_sigmad0calc);

                                                //fdata_sigmad0 = sqrt(pow(((*fvecdata_d0)[t]), 2) + (fdata_wx)*(fdata_wy)); //check the formula again
                                                //fdata_sigmad0 = sqrt(pow((*fvecdata_d0)[t], 2) + 0.01429*0.01727);

                                                //data_d0leaf->Fill(fdata_d0);
                                                ++fdata_trkd0;
                                                ++fdata_multiplicity;
                                                fdata_multiplicity_norm += fdata_multiplicity;
                                            }

                                            /*if ((fabs(fdata_dz_sigmadz) < dz_dzErr_cut) && (fabs(fdata_d0_sigmad0run1) < d0_d0Err_cut))
                                            {
                                                data_sigmapt_pt->Fill(fdata_sigmapt_pt);
                                                ++fdata_trkdpt;
                                            }*/

                                        }
                                    }
                                }
    //========================================================= End of Trk Loop ================================================================
                            }
                            //}
                        }
                    }
    //========================================================= End of Vertex Loop ================================================================
                }
                data_multiplicity->Fill(fdata_multiplicity);
            }
        }
    //========================================================= End of Evt Loop ================================================================

        cout << "Before plotting." << endl;
        cout << "wx is " << fdata_wx << endl;
        cout << "wy is " << fdata_wy << endl;
        cout << "wz is " << fdata_wz << endl;
        cout << "cut on d0 is " << sqrt(0.0025 - fdata_wx*fdata_wy) << endl;
        cout << "Largest z-coordinate of BS is " << fdata_vtxzBSupper << endl;

        TFile data_plot("Histos/data_tree1.root", "recreate");
        //gStyle->SetOptLogy();

        //canvas->Divide (2,2);

        data_eta_histo->Scale(1/fdata_trketa);
        data_eta_histo->SetMaximum(0.04);
        data_eta_histo->GetXaxis()->SetTitle("#eta");
        data_eta_histo->GetYaxis()->SetTitleOffset(1.3);
        data_eta_histo->GetYaxis()->SetTitle("Fraction of Tracks");
        data_eta_histo->Draw();
        data_eta_histo->Write();

        //canvas->cd(3);
        data_phi_histo->Scale(1/fdata_trkphi);
        data_phi_histo->SetMinimum(0.014);
        data_phi_histo->SetMaximum(0.024);
        data_phi_histo->GetXaxis()->SetTitle("#phi");
        data_phi_histo->GetYaxis()->SetTitleOffset(1.3);
        data_phi_histo->GetYaxis()->SetTitle("Fraction of Tracks");
        data_phi_histo->Write();

        //canvas->cd(1);
        gPad->SetLogy();
        data_pt_histo->Scale(1/fdata_trkpt);
        data_pt_histo->GetXaxis()->SetTitle("p_{T}");
        data_pt_histo->GetYaxis()->SetTitleOffset(1.3);
        data_pt_histo->GetYaxis()->SetTitle("Fraction of Tracks");
        data_pt_histo->Write();

        data_dzleaf->DrawNormalized("", 1);
        data_dzleaf->Write();

        data_dz_sigmadz->Scale(1/fdata_trkdz);
        data_dz_sigmadz->SetMinimum(1E-5);
        data_dz_sigmadz->SetMaximum(1E-1);
        data_dz_sigmadz->Draw();
        data_dz_sigmadz->Write();

        data_dz_sigmadzrun1->Scale(1/fdata_trkdz);
        data_dz_sigmadzrun1->Write();

        data_dz_sigmadzcalc->Scale(1/fdata_trkdz);
        data_dz_sigmadzcalc->Write();

        //both with leaf values

        data_d0leaf->Scale(1/fdata_trkd0);
        data_d0leaf->GetXaxis()->SetTitle("d_{0}");
        data_d0leaf->GetYaxis()->SetTitleOffset(1.3);
        data_d0leaf->GetYaxis()->SetTitle("Fraction of Tracks");
        data_d0leaf->Write();

        data_sigmad0run1->Scale(1/fdata_trkd0);
        data_sigmad0run1->GetXaxis()->SetTitle("#sigma_{xy} with Run 1 formula");
        data_sigmad0run1->GetYaxis()->SetTitleOffset(1.25);
        data_sigmad0run1->GetYaxis()->SetTitle("Fraction of Tracks");
        data_sigmad0run1->Write();

        data_d0_sigmad0->SetMinimum(1E-5); //leaf values
        data_d0_sigmad0->SetMaximum(0.1);
        data_d0_sigmad0->Scale(1/fdata_trkd0);
        data_d0_sigmad0->GetXaxis()->SetTitle("d_{0}/#sigma_{0} leaf values");
        data_d0_sigmad0->GetYaxis()->SetTitleOffset(1.2);
        data_d0_sigmad0->GetYaxis()->SetTitle("Fraction of Tracks");
        data_d0_sigmad0->Write();

        data_d0_sigmad0run1->Scale(1/fdata_trkd0); //= d0leaf/sigma_d0run1
        data_d0_sigmad0run1->GetXaxis()->SetTitle("d_{0}/#sigma_{0} Run 1 formula");
        data_d0_sigmad0run1->GetYaxis()->SetTitleOffset(1.2);
        data_d0_sigmad0run1->GetYaxis()->SetTitle("Fraction of Tracks");
        data_d0_sigmad0run1->Write();

        //data_d0_sigmad0calc->Scale(1/fdata_trkd0);
        //data_d0_sigmad0calc->Draw();
        //data_d0_sigmad0calc->Write();

        data_sigmapt_pt->Scale(1/fdata_trkdpt);
        data_sigmapt_pt->SetMinimum(1E-8);
        data_sigmapt_pt->GetXaxis()->SetTitle("#sigma{p_{T}}/p_{T}");
        data_sigmapt_pt->GetYaxis()->SetTitleOffset(1.2);
        data_sigmapt_pt->GetYaxis()->SetTitle("Fraction of Tracks");
        data_sigmapt_pt->Write();

        data_validhits->Scale(1/fdata_trkvalidhits);
        data_validhits->GetXaxis()->SetTitle("# Valid Hits");
        data_validhits->GetYaxis()->SetTitleOffset(1.3);
        data_validhits->GetYaxis()->SetTitle("Fraction of Tracks");
        data_validhits->Write();

        data_chi2n->Scale(1/fdata_trkchi2n);
        data_chi2n->GetXaxis()->SetTitle("#chi^2/dof");
        data_chi2n->GetYaxis()->SetTitleOffset(1.3);
        data_chi2n->GetYaxis()->SetTitle("Fraction of Tracks");
        data_chi2n->Write();

        TCanvas *multiplicity = new TCanvas ("multiplicity", "Multiplicity");
        TPad *data_pad = new TPad ("data_pad", "Data Pad", 0, 0.3, 1.0, 1.0);
        //TPad *ratio_pad = new TPad ("ratio_pad", "Ratio Pad", 0, 0, 1.0, 0.3);

        data_pad->Draw();
        //ratio_pad->Draw();

        data_pad->cd();
        data_multiplicity->Scale(1/fdata_evt);
        data_multiplicity->GetXaxis()->SetTitle("Multiplicity");
        data_multiplicity->GetYaxis()->SetTitleOffset(1.1);
        data_multiplicity->GetYaxis()->SetTitle("Fraction of Tracks");
        data_multiplicity->SetMarkerStyle(8);
        data_multiplicity->SetMarkerColor(kViolet);
        data_multiplicity->SetMarkerSize(1.5);
        data_multiplicity->SetLineWidth(2);
        data_multiplicity->SetLineColor(kGreen);
        gPad->SetTickx();
        gPad->SetTicky();
        data_multiplicity->Draw();
        data_multiplicity->Write();

        TLegend *dataleg_multiplicity = new TLegend (0.6, 0.7, 0.9, 0.9);
        dataleg_multiplicity->SetFillColor(0);
        dataleg_multiplicity->SetBorderSize(0);
        dataleg_multiplicity->SetTextSize(0.02);
        dataleg_multiplicity->AddEntry(data_multiplicity, "data", "lf");
        dataleg_multiplicity->Draw();
        dataleg_multiplicity->Write();

        data_pad->Modified();
        data_pad->Update();
        data_plot.Write();

        /*data_vtxzminusvtxz->Scale(1/fdata_numvtxzminusvtxz);
        data_vtxzminusvtxz->Write();

        data_vtxzposn->Scale(1/fdata_numselectedvtxz);
        data_vtxzposn->Write();*/

        //canvas->cd(2);


        //canvas->cd(4);


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

//==================================================== Reco Loop ===================================================================

    if ((nselect == 1) || (nselect == 2))
    {

        //==============================================Variables==============================================================

        float freco_evt = 0;
        float freco_trkphi = 0;
        float freco_trketa = 0;
        float freco_trkd0 = 0;
        float freco_trkpt = 0;
        float freco_trkdz = 0;
        float freco_trkdpt = 0;
        float freco_dz_sigmadz, freco_d0_sigmad0, freco_sigmapt_pt, freco_dz_sigmadzcalc, freco_d0_sigmad0calc, freco_sigmad0, freco_sigmad0run1, freco_d0, freco_d0_sigmad0run1, freco_trkdx, freco_trkdy;
        float freco_wx, freco_wy,freco_wz, freco_vtxxysize, freco_vtxzsize, freco_sigmad0calc, freco_dz, freco_sigmadz, freco_dz_sigmadzrun1;
        float freco_numberofvtxxBS, freco_vtxxBSvalue, freco_vtxxBSlower, freco_vtxxBSupper;
        float freco_numberofvtxyBS, freco_vtxyBSvalue, freco_vtxyBSlower, freco_vtxyBSupper;
        float freco_numberofvtxzBS, freco_vtxzBSvalue, freco_vtxzBSlower, freco_vtxzBSupper;
        float freco_vtxzminusvtxz, freco_multiplicity;
        //float freco_multiplicity_norm = 0;
        float freco_sqvtxx = 0;
        float freco_sqvtxxnumber = 0;
        float freco_sqvtxy = 0;
        float freco_sqvtxynumber = 0;
        float freco_sqvtxz = 0;
        float freco_sqvtxznumber = 0;
        float freco_numberoftrkdx = 0;
        float freco_numberoftrkdy = 0;
        float freco_numselectedvtxz = 0;
        float freco_numvtxzminusvtxz = 0;
        float freco_trkvalidhits = 0;
        float freco_trkchi2n = 0;
        int nreco_numberofvtxx, nreco_numberofvtxy, nreco_numberofvtxz, nreco_numberofvtxxBS, nreco_numberofvtxyBS, nreco_numberofvtxzBS;
        //float freco_dz_sigmadz = 0;
        //float freco_d0_sigmad0 = 0;
        //float freco_sigmapt_pt = 0;
        int nreco_totaltrk;

//=========================== Retrieve ROOT file ============================

        TFile *recofile = TFile::Open("treesCUETP8M1_66.root", "READ");
        TTree *recotree = (TTree*)recofile->Get("UETree/data");

//=========================== Define variables to read TTree ================

        vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *reco_tracks = 0;
        recotree->SetBranchAddress("recoTracksp4", &reco_tracks);

        int nreco_lumi;
        recotree->SetBranchAddress("lumi", &nreco_lumi);

        int nreco_zerobias;
        recotree->SetBranchAddress("trgZeroBias", &nreco_zerobias);

        vector<float> *fvecreco_vtxx = 0;
        recotree->SetBranchAddress("vtxx", &fvecreco_vtxx);

        vector<float> *fvecreco_vtxy = 0;
        recotree->SetBranchAddress("vtxy", &fvecreco_vtxy);

        vector<float> *fvecreco_vtxz = 0;
        recotree->SetBranchAddress("vtxz", &fvecreco_vtxz);

        vector<float> *fvecreco_vtxxBS = 0;
        recotree->SetBranchAddress("vtxxBS", &fvecreco_vtxxBS);

        vector<float> *fvecreco_vtxyBS = 0;
        recotree->SetBranchAddress("vtxyBS", &fvecreco_vtxyBS);

        vector<float> *fvecreco_vtxzBS = 0;
        recotree->SetBranchAddress("vtxzBS", &fvecreco_vtxzBS);

        vector<int> *nvecreco_highpurity = 0;
        recotree->SetBranchAddress("recoTrackshighPurity", &nvecreco_highpurity);

        vector<int> *nvecreco_vtxndof = 0;
        recotree->SetBranchAddress("vtxndof", &nvecreco_vtxndof);

        vector<float> *fvecreco_dz = 0;
        recotree->SetBranchAddress("recoTracksdz", &fvecreco_dz);

        vector<float> *fvecreco_d0 = 0;
        recotree->SetBranchAddress("recoTracksd0", &fvecreco_d0);

        vector<float> *fvecreco_dzerr = 0;
        recotree->SetBranchAddress("recoTracksdzErr", &fvecreco_dzerr);

        vector<float> *fvecreco_d0err = 0;
        recotree->SetBranchAddress("recoTracksd0Err", &fvecreco_d0err);

        vector<float> *fvecreco_pterr = 0;
        recotree->SetBranchAddress("recoTracksptErr", &fvecreco_pterr);

        vector<float> *fvecreco_vtxxerr = 0;
        recotree->SetBranchAddress("vtxxErr", &fvecreco_vtxxerr);

        vector<float> *fvecreco_vtxyerr = 0;
        recotree->SetBranchAddress("vtxyErr", &fvecreco_vtxyerr);

        vector<float> *fvecreco_vtxzerr = 0;
        recotree->SetBranchAddress("vtxzErr", &fvecreco_vtxzerr);

        vector<int> *nvecreco_validhits = 0;
        recotree->SetBranchAddress("recoTracksnValidHits", &nvecreco_validhits);

        vector<float> *fvecreco_trackschi2n = 0;
        recotree->SetBranchAddress("recoTrackschi2n", &fvecreco_trackschi2n);


    //============================================== Histos for pT, eta, phi ==============================================================

        TH1F *reco_pt_histo = new TH1F ("reco pT", "Normalized reco p_{T}", 200, 0, 10);
        TH1F *reco_eta_histo = new TH1F("reco eta", "Normalised reco #eta", 50, -2.5, 2.5);
        TH1F *reco_phi_histo = new TH1F ("reco phi", "Normalized_reco #phi", 60, -3, 3);

    //================================================== Histos for dz and sigma_dz =========================================================================
        TH1F *reco_dzleaf = new TH1F ("reco_dz", "reco dz", 400, -10, 10);
        TH1F *reco_sigmadzcalc = new TH1F ("reco_sigmadzcalc", "reco #sigma_{z} calc", 200, 0, 4);//plot of sigmadz using formula from WY
        TH1F *reco_sigmadzrun1 = new TH1F ("reco_sigmadzrun1", "reco #sigma_{z} run 1", 200, 0, 4);//plot of sigmadz using formula from run 1
        TH1F *reco_dz_sigmadz = new TH1F ("reco_dz_sigmadz", "reco d_{z}/#sigma_{z}", 160, -20, 20);
        TH1F *reco_dz_sigmadzcalc = new TH1F ("reco_dz_sigmadzcalc", "reco d_{z}/#sigma_{z} calc", 160, -20, 20);
        TH1F *reco_dz_sigmadzrun1 = new TH1F ("reco_dz_sigmadz run 1", "reco d_{z}/#sigma_{z} run 1", 160, -20, 20);
        //TH1F *reco_dz_sigmadzcalcb4cut = new TH1F ("reco_dz_sigmadzcalcb4cut", "reco d_{z}/#sigma_{z} calc", 160, -20, 20);

    //================================================== Histos for d0 and sigma_d0 =========================================================================

        TH1F *reco_d0leaf = new TH1F ("reco_d0_calcrun1", "reco d_{0} calc run 1", 100, -1, 1); //plot of d0 using leaf value to determine how much reco is cut
        TH1F *reco_sigmad0calc = new TH1F ("reco_sigmad0calc", "reco #sigma_{xy} calc", 200, 0, 4);//plot of sigmad0 using formula from WY
        TH1F *reco_sigmad0run1 = new TH1F ("reco_sigmad0run1", "reco #sigma_{xy} run 1", 200, 0, 4);//plot of sigmad0 using run 1 formula

        TH1F *reco_d0_sigmad0 = new TH1F ("reco_d0_sigmad0", "reco d_{0}/#sigma_{xy}", 160, -20, 20); //both leaf values
        TH1F *reco_d0_sigmad0run1 = new TH1F ("reco_d0_sigmad0calcrun1", "reco d_{0}/#sigma_{xy} calc run 1", 160, -20, 20); //plot using run 1 formula
        TH1F *reco_d0_sigmad0calc = new TH1F ("reco_d0_sigmad0calc", "reco d_{0}/#sigma_{xy} calc", 160, -20, 20); //using formula from WY
        //TH1F *reco_d0_sigmad0calcb4cut = new TH1F ("reco_d0_sigmad0calcb4cut", "reco d_{0}/#sigma_{xy} calc", 160, -20, 20);
        //TH1F *reco_d0_sigmad0calc = new TH1F ("reco_d0_sigmad0calc", "reco d_{0}/#sigma_{xy} calc", 150, 0, 0.5);

    //================================================== Histos for pT and sigma_pT =========================================================================

        TH1F *reco_sigmapt_pt = new TH1F ("reco_sigmapt_pt", "reco #sigma_{p_{T}}/p_{T}", 20, 0, 0.2);

    //================================================== Histos for chi2n and ValidHits =========================================================================

        TH1F *reco_validhits = new TH1F ("Tracks_vs_validhits", "Tracks vs validhits", 50, 0, 50);
        TH1F *reco_chi2n = new TH1F ("Tracks_vs_chi2n", "Tracks vs #chi^{2/ndof}", 50, 0, 5);

    //================================================== Histos for Multiplicity =========================================================================

        //TH1F *reco_multiplicity = new TH1F("Normalized_Multiplicity", "Normalized Multiplicity", 200, 0, 200);
        TH1F *reco_vtxzminusvtxz = new TH1F ("vtxzminusvtxz", "vtxzminusvtxz", 800, -4, 4);
        TH1F *reco_vtxzposn = new TH1F ("vtxzposn", "vtxzpos^{n}", 400, -20, 20);


        Int_t nreco_totalEvt = (Int_t)recotree->GetEntries();
        cout << "There is a total of " << nreco_totalEvt << " events." << endl;

        for (int ev = 0; ev < nreco_totalEvt; ++ev)
        {
            recotree->GetEntry(ev);

            if (nreco_lumi >= lumi_cut)
            {
                if (nreco_zerobias == 1)
                {
                    nreco_numberofvtxx = fvecreco_vtxx->size();
                    nreco_numberofvtxy = fvecreco_vtxy->size();
                    nreco_numberofvtxz = fvecreco_vtxz->size();
                    nreco_numberofvtxxBS = fvecreco_vtxxBS->size();
                    nreco_numberofvtxyBS = fvecreco_vtxyBS->size();
                    nreco_numberofvtxzBS = fvecreco_vtxzBS->size();

                    freco_vtxxBSlower = *min_element(fvecreco_vtxxBS->begin(), fvecreco_vtxxBS->end());
                    freco_vtxxBSupper = *max_element(fvecreco_vtxxBS->begin(), fvecreco_vtxxBS->end());
                    freco_vtxyBSlower = *min_element(fvecreco_vtxyBS->begin(), fvecreco_vtxyBS->end());
                    freco_vtxyBSupper = *max_element(fvecreco_vtxyBS->begin(), fvecreco_vtxyBS->end());
                    freco_vtxzBSlower = *min_element(fvecreco_vtxzBS->begin(), fvecreco_vtxzBS->end());
                    freco_vtxzBSupper = *max_element(fvecreco_vtxzBS->begin(), fvecreco_vtxzBS->end());

                    for(int vtxx = 0; vtxx != nreco_numberofvtxx; ++vtxx)
                    {
                        //freco_sqvtxx += pow(((*fvecreco_vtxx)[vtxx]) - ((*fvecreco_vtxxBS)[vtxx]), 2);
                        freco_sqvtxx += pow((((*fvecreco_vtxx)[vtxx]) - (freco_vtxxBSlower)), 2);
                        ++freco_sqvtxxnumber;
                    }

                    for(int vtxy = 0; vtxy != nreco_numberofvtxy; ++vtxy)
                    {
                        //freco_sqvtxy += pow(((*fvecreco_vtxy)[vtxy]) - ((*fvecreco_vtxyBS)[vtxy]), 2);
                        freco_sqvtxy += pow((((*fvecreco_vtxy)[vtxy]) - (freco_vtxyBSupper)), 2);
                        ++freco_sqvtxynumber;
                    }

                    for(int vtxz = 0; vtxz != nreco_numberofvtxz; ++vtxz)
                    {
                        //freco_sqvtxz += pow(((*fvecreco_vtxz)[vtxz]) - ((*fvecreco_vtxzBS)[vtxz]), 2);
                        freco_sqvtxz += pow((((*fvecreco_vtxz)[vtxz]) - (freco_vtxzBSupper)), 2);
                        ++freco_sqvtxznumber;
                    }

                    for (int vtxzBS = 0; vtxzBS != nreco_numberofvtxzBS; ++vtxzBS)
                    {
                        if ((*fvecreco_vtxzBS)[vtxzBS] == freco_vtxzBSupper)
                        {
                            ++freco_numberofvtxzBS;
                        }

                    }
                }
            }
        }
        cout << "Largest BS z-coordinate is " << freco_vtxzBSupper << endl;
        freco_wx = sqrt((freco_sqvtxx) / (freco_sqvtxxnumber)); //RMS(?) of position of vtx x-coordinate. Averaged over number of vtx
        freco_wy = sqrt((freco_sqvtxy) / (freco_sqvtxynumber)); //RMS(?) of position of vtx y-coordinate Averaged over number of vtx
        freco_wz = sqrt((freco_sqvtxz) / (freco_sqvtxznumber)); //

    //========================================================= Start of Evt Loop ================================================================

        for (Int_t i = 0; i < nreco_totalEvt; ++i)
        {
            recotree->GetEntry(i);
            //cout << "At entry " << i << endl;

            if (nreco_lumi >= lumi_cut)
            {
                if (nreco_zerobias == 1)
                {
                    nreco_totaltrk = reco_tracks->size();
                    ++freco_evt;
                    freco_multiplicity = 0;

                    int vtxdof = 0;

                    nreco_numberofvtxxBS = fvecreco_vtxxBS->size();
                    nreco_numberofvtxyBS = fvecreco_vtxyBS->size();
                    nreco_numberofvtxzBS = fvecreco_vtxzBS->size();

    //========================================================= Start of Vertex Loop ================================================================

                    for (int vtxnumber = 0; vtxnumber != nreco_numberofvtxxBS; ++vtxnumber)
                    {
                        if((*nvecreco_vtxndof)[vtxdof] > dof_cut)
                        {
                            freco_vtxxysize = sqrt(pow(fabs(((*fvecreco_vtxx)[vtxnumber]) - ((*fvecreco_vtxxBS)[vtxnumber])), 2) + pow(fabs((*fvecreco_vtxy)[vtxnumber] - (*fvecreco_vtxyBS)[vtxnumber]), 2));
                            freco_vtxzsize = fabs(((*fvecreco_vtxz)[vtxnumber]) - ((*fvecreco_vtxzBS)[vtxnumber]));

                            //if (freco_vtxxysize <= vtxxysize && freco_vtxzsize <= vtxzsize)
                            //{
                                if (nreco_numberofvtxzBS == vtx_number_cut)
                                {
                                    reco_vtxzposn->Fill((*fvecreco_vtxzBS)[vtxnumber]);
                                    ++freco_numselectedvtxz;

    //========================================================= Start of Trk Loop ================================================================

                                    for (int t = 0; t != nreco_totaltrk; ++t)
                                    {
                                        XYZTVector reco_vec = (*reco_tracks)[t];
                                        //cout << "Within track " << endl;
                                        //using formula from paper of run 1 result, dz is leaf value

                            //======================================= dz ===================================================

                                        freco_dz = (*fvecreco_dz)[t];
                                        reco_dzleaf->Fill((*fvecreco_dz)[t]);
                                        freco_sigmadz = sqrt(pow(((*fvecreco_dzerr)[t]), 2) + pow(freco_wz, 2));

                                        freco_dz_sigmadz = ((*fvecreco_dz)[t])/((*fvecreco_dzerr)[t]); //leaf
                                        freco_dz_sigmadzcalc = ((freco_dz)/sqrt(pow(((*fvecreco_dzerr)[t]),2)+pow((*fvecreco_vtxzerr)[t],2))); //WY
                                        freco_dz_sigmadzrun1 = ((freco_dz) / (freco_sigmadz)); //run 1 formula

                            //======================================= d0 ===================================================

                                        freco_d0 = (*fvecreco_d0)[t];
                                        freco_sigmad0run1 = sqrt(pow(((*fvecreco_d0err)[t]), 2) + (freco_wx)*(freco_wy));
                                        freco_sigmad0calc = sqrt(pow(((*fvecreco_d0err)[t]),2)+pow((*fvecreco_vtxxerr)[t],2)+pow((*fvecreco_vtxyerr)[t],2)); //WY

                                        freco_d0_sigmad0 = ((*fvecreco_d0)[t])/((*fvecreco_d0err)[t]); //both leaf values
                                        freco_d0_sigmad0run1 = (((*fvecreco_d0)[t]) / (freco_sigmad0run1)); //run 1 formula
                                        freco_d0_sigmad0calc = (((*fvecreco_d0)[t]) / (freco_sigmad0calc)); //with sigmad0 from WY

                                        //freco_sigmad0 = fmin(0.05,(freco_sigmad0run1));

                                        //freco_d0_sigmad0 = ((*fvecreco_d0)[t])/sqrt((pow((*fvecreco_vtxxerr)[t],2)+pow((*fvecreco_vtxyerr)[t],2)));
                                        //freco_d0_sigmad0run1 = ((freco_d0) / (freco_sigmad0));

                            //======================================= pT ===================================================

                                        freco_sigmapt_pt = (((*fvecreco_pterr)[t])/(reco_vec.Pt())); //both leaf values

                            //======================================= Trk cuts ===================================================

                                        if ((*nvecreco_highpurity)[t] == 1)
                                        {
                                            if (fabs(reco_vec.Eta()) <= eta_cut)
                                            {
                                                reco_pt_histo->Fill(reco_vec.Pt());
                                                ++freco_trkpt;
                                            }

                                            //if (reco_vec.Pt() >= pt_cut && fabs(freco_sigmapt_pt) < ptErr_pt_cut && fabs(freco_dz_sigmadzrun1 < dz_dzErr_cut))
                                            //if ((reco_vec.Pt() >= pt_cut) && (fabs(freco_sigmapt_pt) < ptErr_pt_cut))
                                            if (reco_vec.Pt() >= pt_cut)
                                            {
                                                reco_eta_histo->Fill(reco_vec.Eta());
                                                ++freco_trketa;
                                            }

                                            if ((reco_vec.Pt() >= pt_cut) && (fabs(reco_vec.Eta()) <= eta_cut))
                                            {
                                                reco_phi_histo->Fill(reco_vec.Phi());
                                                ++freco_trkphi;

                                                reco_sigmapt_pt->Fill(freco_sigmapt_pt);
                                                ++freco_trkdpt;

                                                reco_validhits->Fill((*nvecreco_validhits)[t]);
                                                ++freco_trkvalidhits;

                                                reco_chi2n->Fill((*fvecreco_trackschi2n)[t]);
                                                ++freco_trkchi2n;

                                                //if (fabs(freco_d0_sigmad0calc) < d0_d0Err_cut && fabs(freco_sigmapt_pt) < ptErr_pt_cut && freco_d0 < fmin(0.05, freco_sigmad0calc)))
                                                //if (fabs(freco_d0_sigmad0calc) < d0_d0Err_cut && fabs(freco_sigmapt_pt) < ptErr_pt_cut)
                                                if (fabs(freco_sigmapt_pt) < ptErr_pt_cut)
                                                {
                                                    reco_dz_sigmadz->Fill(freco_dz_sigmadz);
                                                    reco_dz_sigmadzcalc->Fill(freco_dz_sigmadzcalc);
                                                    reco_dz_sigmadzrun1->Fill(freco_dz_sigmadzrun1);
                                                    ++freco_trkdz;
                                                }

                                                //if (fabs(freco_dz_sigmadzrun1) < dz_dzErr_cut && fabs(freco_sigmapt_pt) < ptErr_pt_cut)
                                                if (fabs(freco_sigmapt_pt) < ptErr_pt_cut)
                                                {
                                                    reco_d0leaf->Fill((*fvecreco_d0)[t]);
                                                    reco_sigmad0run1->Fill(freco_sigmad0run1);
                                                    reco_d0_sigmad0->Fill(freco_d0_sigmad0);
                                                    reco_d0_sigmad0run1->Fill(freco_d0_sigmad0run1);
                                                    reco_d0_sigmad0calc->Fill(freco_d0_sigmad0calc);

                                                    //freco_sigmad0 = sqrt(pow(((*fvecreco_d0)[t]), 2) + (freco_wx)*(freco_wy)); //check the formula again
                                                    //freco_sigmad0 = sqrt(pow((*fvecreco_d0)[t], 2) + 0.01429*0.01727);

                                                    //reco_d0leaf->Fill(freco_d0);
                                                    ++freco_trkd0;
                                                    ++freco_multiplicity;
                                                    freco_multiplicity_norm += freco_multiplicity;
                                                }

                                                /*if ((fabs(freco_dz_sigmadz) < dz_dzErr_cut) && (fabs(freco_d0_sigmad0run1) < d0_d0Err_cut))
                                                {
                                                    reco_sigmapt_pt->Fill(freco_sigmapt_pt);
                                                    ++freco_trkdpt;
                                                }*/

                                            }
                                        }
                                    }

    //========================================================= End of Trk Loop ================================================================
                                }
                            //}
                        }
                    }
    //========================================================= End of Vertex Loop ================================================================
                    reco_multiplicity->Fill(freco_multiplicity);
                }
            }
        }
    //========================================================= End of Evt Loop ================================================================

        cout << "Before plotting." << endl;
        cout << "wx is " << freco_wx << endl;
        cout << "wy is " << freco_wy << endl;
        cout << "wz is " << freco_wz << endl;
        cout << "cut on d0 is " << sqrt(0.0025 - freco_wx*freco_wy) << endl;
        cout << "Largest z-coordinate of BS is " << freco_vtxzBSupper << endl;

        TFile reco_plot("Histos/data_treesCUETP8M1_66.root", "recreate");
        //gStyle->SetOptLogy();

        //canvas->Divide (2,2);

        reco_eta_histo->Scale(1/freco_trketa);
        reco_eta_histo->SetMaximum(0.04);
        reco_eta_histo->GetXaxis()->SetTitle("#eta");
        reco_eta_histo->GetYaxis()->SetTitleOffset(1.3);
        reco_eta_histo->GetYaxis()->SetTitle("Fraction of Tracks");
        reco_eta_histo->SetLineColor(4);
        reco_eta_histo->Draw();
        reco_eta_histo->Write();

        //canvas->cd(3);
        reco_phi_histo->Scale(1/freco_trkphi);
        reco_phi_histo->SetMinimum(0.014);
        reco_phi_histo->SetMaximum(0.024);
        reco_phi_histo->GetXaxis()->SetTitle("#phi");
        reco_phi_histo->GetYaxis()->SetTitleOffset(1.3);
        reco_phi_histo->GetYaxis()->SetTitle("Fraction of Tracks");
        reco_phi_histo->SetLineColor(4);
        reco_phi_histo->Write();

        //canvas->cd(1);
        gPad->SetLogy();
        reco_pt_histo->Scale(1/freco_trkpt);
        reco_pt_histo->GetXaxis()->SetTitle("p_{T}");
        reco_pt_histo->GetYaxis()->SetTitleOffset(1.3);
        reco_pt_histo->GetYaxis()->SetTitle("Fraction of Tracks");
        reco_pt_histo->SetLineColor(4);
        reco_pt_histo->Write();

        reco_dzleaf->DrawNormalized("", 1);
        reco_dzleaf->SetLineColor(4);
        reco_dzleaf->Write();

        reco_dz_sigmadz->Scale(1/freco_trkdz);
        reco_dz_sigmadz->SetMinimum(1E-5);
        reco_dz_sigmadz->SetMaximum(1E-1);
        reco_dz_sigmadz->SetLineColor(4);
        reco_dz_sigmadz->Draw();
        reco_dz_sigmadz->Write();

        reco_dz_sigmadzrun1->Scale(1/freco_trkdz);
        reco_dz_sigmadzrun1->SetLineColor(4);
        reco_dz_sigmadzrun1->Write();

        reco_dz_sigmadzcalc->Scale(1/freco_trkdz);
        reco_dz_sigmadzcalc->SetLineColor(4);
        reco_dz_sigmadzcalc->Write();

        //both with leaf values

        reco_d0leaf->Scale(1/freco_trkd0);
        reco_d0leaf->GetXaxis()->SetTitle("d_{0}");
        reco_d0leaf->GetYaxis()->SetTitleOffset(1.3);
        reco_d0leaf->GetYaxis()->SetTitle("Fraction of Tracks");
        reco_d0leaf->SetLineColor(4);
        reco_d0leaf->Write();

        reco_sigmad0run1->Scale(1/freco_trkd0);
        reco_sigmad0run1->GetXaxis()->SetTitle("#sigma_{xy} with Run 1 formula");
        reco_sigmad0run1->GetYaxis()->SetTitleOffset(1.25);
        reco_sigmad0run1->GetYaxis()->SetTitle("Fraction of Tracks");
        reco_sigmad0run1->SetLineColor(4);
        reco_sigmad0run1->Write();

        reco_d0_sigmad0->SetMinimum(1E-5); //leaf values
        reco_d0_sigmad0->SetMaximum(0.1);
        reco_d0_sigmad0->Scale(1/freco_trkd0);
        reco_d0_sigmad0->GetXaxis()->SetTitle("d_{0}/#sigma_{0} leaf values");
        reco_d0_sigmad0->GetYaxis()->SetTitleOffset(1.2);
        reco_d0_sigmad0->GetYaxis()->SetTitle("Fraction of Tracks");
        reco_d0_sigmad0->SetLineColor(4);
        reco_d0_sigmad0->Write();

        reco_d0_sigmad0run1->Scale(1/freco_trkd0); //= d0leaf/sigma_d0run1
        reco_d0_sigmad0run1->GetXaxis()->SetTitle("d_{0}/#sigma_{0} Run 1 formula");
        reco_d0_sigmad0run1->GetYaxis()->SetTitleOffset(1.2);
        reco_d0_sigmad0run1->GetYaxis()->SetTitle("Fraction of Tracks");
        reco_d0_sigmad0run1->SetLineColor(4);
        reco_d0_sigmad0run1->Write();

        //reco_d0_sigmad0calc->Scale(1/freco_trkd0);
        //reco_d0_sigmad0calc->Draw();
        //reco_d0_sigmad0calc->Write();

        reco_sigmapt_pt->Scale(1/freco_trkdpt);
        reco_sigmapt_pt->SetMinimum(1E-8);
        reco_sigmapt_pt->GetXaxis()->SetTitle("#sigma{p_{T}}/p_{T}");
        reco_sigmapt_pt->GetYaxis()->SetTitleOffset(1.2);
        reco_sigmapt_pt->GetYaxis()->SetTitle("Fraction of Tracks");
        reco_sigmapt_pt->SetLineColor(4);
        reco_sigmapt_pt->Write();

        reco_validhits->Scale(1/freco_trkvalidhits);
        reco_validhits->GetXaxis()->SetTitle("# Valid Hits");
        reco_validhits->GetYaxis()->SetTitleOffset(1.3);
        reco_validhits->GetYaxis()->SetTitle("Fraction of Tracks");
        reco_validhits->SetLineColor(4);
        reco_validhits->Write();

        reco_chi2n->Scale(1/freco_trkchi2n);
        reco_chi2n->GetXaxis()->SetTitle("#chi^2/dof");
        reco_chi2n->GetYaxis()->SetTitleOffset(1.3);
        reco_chi2n->GetYaxis()->SetTitle("Fraction of Tracks");
        reco_chi2n->SetLineColor(4);
        reco_chi2n->Write();

        TPad *reco_pad = new TPad ("reco_pad", "Data Pad", 0, 0.3, 1.0, 1.0);
        //TPad *ratio_pad = new TPad ("ratio_pad", "Ratio Pad", 0, 0, 1.0, 0.3);

        reco_pad->Draw();
        //ratio_pad->Draw();

        reco_pad->cd();
        reco_multiplicity->Scale(1/freco_multiplicity_norm);
        reco_multiplicity->GetXaxis()->SetTitle("Multiplicity");
        reco_multiplicity->GetYaxis()->SetTitleOffset(1.1);
        reco_multiplicity->GetYaxis()->SetTitle("Fraction of Tracks");
        reco_multiplicity->SetMarkerStyle(8);
        reco_multiplicity->SetMarkerColor(kViolet);
        reco_multiplicity->SetMarkerSize(1.5);
        reco_multiplicity->SetLineWidth(2);
        reco_multiplicity->SetLineColor(4);
        gPad->SetTickx();
        gPad->SetTicky();
        reco_multiplicity->Draw();
        reco_multiplicity->Write();

        TLegend *recoleg_multiplicity = new TLegend (0.6, 0.7, 0.9, 0.9);
        recoleg_multiplicity->SetFillColor(0);
        recoleg_multiplicity->SetBorderSize(0);
        recoleg_multiplicity->SetTextSize(0.02);
        recoleg_multiplicity->AddEntry(reco_multiplicity, "reco", "lf");
        recoleg_multiplicity->Draw();
        recoleg_multiplicity->Write();

        reco_pad->Modified();
        reco_pad->Update();
        reco_plot.Write();

        /*reco_vtxzminusvtxz->Scale(1/freco_numvtxzminusvtxz);
        reco_vtxzminusvtxz->Write();

        reco_vtxzposn->Scale(1/freco_numselectedvtxz);
        reco_vtxzposn->Write();*/

        //canvas->cd(2);


        //canvas->cd(4);


        //reco_eta_histo->DrawNormalized("", 1);
        //reco_eta_histo->Write();
        //reco_eta_histo->Scale(1/freco_trketa);
        //reco_eta_histo->GetYaxis()->SetRange(0,0.04);
        //reco_eta_histo->Draw();
        //reco_eta_histo->Write();

        //canvas->cd(2);
        //gPad->SetLogy();
        //reco_eta_histo->DrawNormalized("", 1);
        //reco_eta_histo->Draw();

        /*canvas->cd(3);
        gPad->SetLogy();
        reco_dz_sigmadz->Scale(1/freco_trk);
        reco_dz_sigmadz->Draw();

        canvas->cd(4);
        gPad->SetLogy();
        reco_pt_histo->Scale(1/freco_trk);
        reco_pt_histo->Draw();*/
    }

    TFile data_reco ("Histos/data_reco.root", "recreate");
    TCanvas *data_reco_canvas = new TCanvas;
    //data_reco_canvas->Divide(2,1);
    TPad *data_reco_pad = new TPad ("data_reco_pad", "Data and Reco Pad", 0, 0.3, 1.0, 1.0);
    TPad *divide_pad = new TPad ("ratio_pad", "Data/MC", 0, 0, 1, 0.3);

    data_reco_pad->cd();
    data_multiplicity->Scale(1/fdata_multiplicity_norm);
    data_multiplicity->GetXaxis()->SetTitle("Data Multiplicity");
    data_multiplicity->GetYaxis()->SetTitleOffset(1.1);
    data_multiplicity->GetYaxis()->SetTitle("Fraction of Tracks");
    data_multiplicity->SetMarkerStyle(8);
    data_multiplicity->SetMarkerColor(6);
    data_multiplicity->SetMarkerSize(1.5);
    data_multiplicity->SetLineWidth(2);
    data_multiplicity->SetLineColor(6);
    gPad->SetTickx();
    gPad->SetTicky();
    data_multiplicity->Draw();
    data_multiplicity->Write();

    reco_multiplicity->Scale(1/freco_multiplicity_norm);
    reco_multiplicity->GetXaxis()->SetTitle("Reco Multiplicity");
    reco_multiplicity->GetYaxis()->SetTitleOffset(1.1);
    reco_multiplicity->GetYaxis()->SetTitle("Fraction of Tracks");
    reco_multiplicity->SetMarkerStyle(8);
    reco_multiplicity->SetMarkerColor(4);
    reco_multiplicity->SetMarkerSize(1.5);
    reco_multiplicity->SetLineWidth(2);
    reco_multiplicity->SetLineColor(4);
    gPad->SetTickx();
    gPad->SetTicky();
    reco_multiplicity->Draw("same");
    reco_multiplicity->Write();

    TLegend *leg_multiplicity = new TLegend (0.6, 0.7, 0.9, 0.9);
    leg_multiplicity->SetFillStyle(0);
    leg_multiplicity->SetBorderSize(0);
    leg_multiplicity->SetTextSize(0.07);
    leg_multiplicity->AddEntry(data_multiplicity, "data", "lf");
    leg_multiplicity->AddEntry(reco_multiplicity, "reco", "lf");
    leg_multiplicity->Draw();
    leg_multiplicity->Write();
    data_reco_pad->Modified();
    data_reco_pad->Update();

    TH1F *data_reco_multiplicity = (TH1F*)data_multiplicity->Clone("data_reco_multiplicity");

    divide_pad->cd();
    data_reco_multiplicity->Divide(reco_multiplicity);
    data_reco_multiplicity->GetXaxis()->SetTitle("");
    data_reco_multiplicity->GetXaxis()->SetLabelSize(0.01);
    data_reco_multiplicity->SetMaximum(5);
    data_reco_multiplicity->SetMinimum(0);
    data_reco_multiplicity->GetYaxis()->SetTitle("Data/Reco");
    data_reco_multiplicity->GetYaxis()->SetTitleOffset(1.2);
    data_reco_multiplicity->GetYaxis()->SetLabelSize(0.1);
    data_reco_multiplicity->GetYaxis()->SetTitleSize(0.15);
    data_reco_multiplicity->SetLineColor(9);
    gPad->SetTickx();
    gPad->SetTicky();
    gPad->SetGrid();
    data_reco_multiplicity->Draw();
    data_reco_multiplicity->Write();
    divide_pad->Modified();
    divide_pad->Update();

    data_reco_canvas->Draw();
    data_reco_canvas->Write();
}





