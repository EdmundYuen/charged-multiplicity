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

void track_selection_data_MC()
{
   //gROOT->ProcessLine(".L Loader.C+");
   //TFile *myFile = TFile::Open("root://eoscms.cern.ch//eos/cms/store/user/wei/multiplicity/trees_10.root", "READ")
   //TTree* tree = (TTree*)myFile->Get("UETree/data");

    //---------------------Cuts-------------------------------------------

    const int lumi_cut = 90;
    const double eta_cut = 2.4;
    const double pt_cut = 0.5;
    const double vtxz_number = 1.;
    const double vtxz_size = 10;
    const float d0_d0Err_cut = 3;
    const float dz_dzErr_cut = 3;
    const float ptErr_pt_cut = 0.05;

    int ndata_MC = 0;
    cout << "Enter 0 for data, 1 for MC and 2 for both : ";
    cin >> ndata_MC;

    if (ndata_MC == 0 || ndata_MC == 2)
    {
        TFile *myFile = TFile::Open("tree1.root", "READ");
        TTree *tree = (TTree*)myFile->Get("UETree/data");

        //for this dataset we want lumisection of 90 and above
        //data from branch 62
        int ndata_lumi_section;
        tree->SetBranchAddress("lumi", &ndata_lumi_section);

        int ndata_ZeroBias; //data from Branch 60
        tree->SetBranchAddress("trgZeroBias",&ndata_ZeroBias);

    //---------------------Variables for data tracks-------------------------------------------

        vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *data_tracks = 0;
        tree->SetBranchAddress("recoTracksp4", &data_tracks);
        //vector<TLorentzVector> *tracks = new vector<TLorentzVector>;

        //variables to check if there is only 1 vertex
        vector<float> *fvec_data_vtxz = 0; //require initialisation to 0 to avoid crash
        tree->SetBranchAddress("vtxz",&fvec_data_vtxz);
        float fvec_data_vtxz_size;

        vector<float> *fvec_data_vtxzBS = 0; //from branch 54
        tree->SetBranchAddress("vtxzBS", &fvec_data_vtxzBS);
        float fdata_BeamSize;

        vector<float> *fvec_data_dz = 0;
        tree->SetBranchAddress ("recoTracksdz", &fvec_data_dz);

        vector<float> *fvec_data_dzErr = 0;
        tree->SetBranchAddress ("recoTracksdzErr", &fvec_data_dzErr);

        vector<float> *fvec_data_d0 = 0;
        tree->SetBranchAddress ("recoTracksd0", &fvec_data_d0);

        vector<float> *fvec_data_d0Err = 0;
        tree->SetBranchAddress ("recoTracksd0Err", &fvec_data_d0Err);

        vector<float> *fvec_data_ptErr = 0;
        tree->SetBranchAddress("recoTracksptErr", &fvec_data_ptErr);

        vector<int> *nvec_data_HighPurity = 0;
        tree->SetBranchAddress ("recoTrackshighPurity", &nvec_data_HighPurity);

        vector<float> *fvec_data_vtxy = 0; //require initialisation to 0 to avoid crash
        tree->SetBranchAddress("vtxy",&fvec_data_vtxy);

        vector<float> *fvec_data_vtxx = 0; //require initialisation to 0 to avoid crash
        tree->SetBranchAddress("vtxx",&fvec_data_vtxx);

        vector<float> *fvec_data_trackschi2n = 0;
        tree->SetBranchAddress("recoTrackschi2n", &fvec_data_trackschi2n);
        float fdata_chi2n;

        vector<int> *nvec_data_validhits = 0;
        tree->SetBranchAddress("recoTracksnValidHits", &nvec_data_validhits);
        int ndata_validhits;

        vector<int> *nvec_data_vtxndof = 0;
        tree->SetBranchAddress("vtxndof", &nvec_data_vtxndof);

        int ndata_totaltrk;
        double ddata_trk;
        double ddata_trkdzErr = 0;
        double ddata_trkdzErrcalc = 0;
        double ddata_trkd0Err = 0;
        double ddata_trkd0Errcalc = 0;
        double ddata_trkptErr = 0;
        double ddata_trkptErrcalc = 0;
        double ddata_trkpt = 0;
        double ddata_trketa = 0;
        double ddata_trkphi = 0;
        double ddata_evt = 0;
        int ndata_trk_normalized = 0;
        int ndata_HighPurity = 0;
        float fdata_dz_dzErr, fdata_d0_d0Err, fdata_ptErr_pt, fdata_dz_dzErrcalc, fdata_d0_d0Errcalc;

        double ddata_dzleaf, ddata_dzcalc, ddata_dzdiff;
        double ddata_d0leaf, ddata_d0calc, ddata_d0diff;

        //------------------------Histograms for data Tracks-------------------------------------------
        TCanvas *data_canvas = new TCanvas;
        //TCanvas *gen_canvas = new TCanvas;
        //TH1F *event_histo = new TH1F ("reco_evt", "reco_evt", 100, 0, 200);
        TH1F *data_pt_histo = new TH1F ("data p_{T}", "Normalized data p_{T}", 80, 0, 4);
        TH1F *data_eta_histo = new TH1F ("data #eta", "Normalized_data #eta", 60, -3, 3);
        TH1F *data_phi_histo = new TH1F ("data #phi", "Normalized_data #phi", 80, -4, 4);
        TH1F *data_lumi_histo = new TH1F ("lumi_section", "lumi_section", 160, 80, 230);
        TH1F *data_multiplicity = new TH1F ("data_multiplicity", "Normalized_data_Multiplicity", 200, 0, 200);
        TH1F *data_dz_sigmadz = new TH1F ("data_dz_sigmadz", "data d_{z}/#sigma_{z}", 160, -20, 20);
        TH1F *data_dz_sigmadzcalc = new TH1F ("data_dz_sigmadz calc", "data d_{z}/#sigma_{z} calc", 160, -20, 20);
        TH1F *data_d0_sigmad0 = new TH1F ("data_d0_sigmad0", "data d_{0}/#sigma_{xy}", 160, -20, 20);
        TH1F *data_d0_sigmad0calc = new TH1F ("data_d0_sigmad0calc", "d_{0}/#sigma_{xy} calc", 160, -20, 20);
        TH1F *data_sigmapt_pt = new TH1F ("data_sigmapt_pt", "data #sigma_{p_{T}}/p_{T}", 20, 0, 0.2);
        //TH1F *normalized_multiplicity_histo = new TH1F ("normalized_multiplicity", "normalized_multiplicity", 200, 0, 200);
        TH1F *data_validhits = new TH1F ("data_validhits", "data_validhits", 100, 0, 50);
        TH1F *data_trackschi2 = new TH1F ("data tracks #chi^{2}", "data tracks #chi^{2}", 80, 0, 40);

        Int_t ndata_totalEvt = (Int_t)tree->GetEntries();
        cout << "Tree Entries " << ndata_totalEvt << endl;

    //start loop over events
        for(Int_t i = 0; i != ndata_totalEvt; ++i)
        {
            //if (i%10000 == 0)
                //cout<< "Processing event number: " << i << endl; //1 TTree entry = 1 event
            cout << "Entry " << i << endl;
            tree -> GetEntry(i);

            //selects for events from lumisection >= 90
            if (ndata_lumi_section >= lumi_cut)
            {
                cout << "lumisection is " << ndata_lumi_section << endl;


                //we select only events that are triggered by ZeroBias Trigger. "True" in ZeroBias evaluate to 1
                if (ndata_ZeroBias == 1)
                {
                    fvec_data_vtxz_size = fvec_data_vtxz->size();
                    //vtxz_plot->Fill(fvecVtxz_size);
                    //ntrk_MC = MC_tracks->size();
                    //nMC_Multi = 0;

                    //---------------------Loops for data tracks-------------------------------------------

                    if (fvec_data_vtxz_size == vtxz_number)
                    {
                        //looping over vertices
                        for (int k = 0; k != fvec_data_vtxz_size; ++k)
                        {
                            fdata_BeamSize = fabs((*fvec_data_vtxz)[k] - (*fvec_data_vtxzBS)[0]);
                            cout << "Beam Size is " << fdata_BeamSize << endl;

                            if (fdata_BeamSize <= vtxz_size)
                            {
                                ndata_totaltrk = data_tracks->size();

                                //fill with the lumi sections that meet the above event-level cuts
                                data_lumi_histo->Fill(ndata_lumi_section);
                                ++ddata_evt;

                                if((*nvec_data_vtxndof)[0] > 4)
                                {
                                    //looping over tracks
                                    for (int j = 0; j != ndata_totaltrk; ++j)
                                    {
                                        XYZTVector data_vec = (*data_tracks)[j];
                                        ddata_dzleaf = (*fvec_data_dz)[j];
                                        ddata_dzcalc = ((data_vec.Z())-(*fvec_data_vtxz)[k])-((((data_vec.X())-(*fvec_data_vtxx)[k])*(data_vec.px())+((data_vec.Y())-(*fvec_data_vtxy)[k])*(data_vec.py()))/(data_vec.Pt())*(data_vec.pz()/data_vec.Pt()));
                                        ddata_dzdiff = ((data_vec.Z())-(*fvec_data_vtxz)[k]);
                                        //tr_d0= (- (tr_x-vtx_x)  track.py() + (tr_y-vtx_y)  track.px() ) / track.pt()
                                        ddata_d0leaf = (*fvec_data_d0)[j];
                                        ddata_d0calc = ((-(data_vec.X() - (*fvec_data_vtxx)[k])*data_vec.Py()) + ((data_vec.Y()-(*fvec_data_vtxy)[k])*data_vec.Px()))/data_vec.Pt();
                                        ddata_trk = 0;

                                        if ((*nvec_data_HighPurity)[j] == 1)
                                        {
                                            if (abs (data_vec.Eta()) <= eta_cut && data_vec.Pt() >= pt_cut)
                                            {
                                                fdata_dz_dzErr = ((*fvec_data_dz)[j])/((*fvec_data_dzErr)[j]);
                                                fdata_d0_d0Err = ((*fvec_data_d0)[j])/((*fvec_data_d0Err)[j]);
                                                fdata_ptErr_pt = ((*fvec_data_ptErr)[j]/(data_vec.Pt()));
                                                fdata_dz_dzErrcalc = ddata_dzcalc / ((*fvec_data_dzErr)[j]);
                                                fdata_d0_d0Errcalc = ddata_d0calc / ((*fvec_data_d0Err)[j]);

                                                if (abs(fdata_ptErr_pt) < ptErr_pt_cut && abs(fdata_d0_d0Err) < d0_d0Err_cut)
                                                {
                                                    data_dz_sigmadz->Fill(fdata_dz_dzErr);
                                                    ++ddata_trkdzErr;
                                                }

                                                if (abs(fdata_ptErr_pt) < ptErr_pt_cut && abs(fdata_d0_d0Errcalc) < d0_d0Err_cut)
                                                {
                                                    data_dz_sigmadzcalc->Fill(fdata_dz_dzErrcalc);
                                                    ++ddata_trkdzErrcalc;
                                                }

                                                if (abs(fdata_ptErr_pt) < ptErr_pt_cut && abs(fdata_dz_dzErr) < dz_dzErr_cut)
                                                {
                                                    data_d0_sigmad0->Fill(fdata_d0_d0Err);
                                                    ++ddata_trkd0Err;
                                                }

                                                if (abs(fdata_ptErr_pt) < ptErr_pt_cut && abs(fdata_dz_dzErrcalc) < dz_dzErr_cut)
                                                {
                                                    data_d0_sigmad0calc->Fill(fdata_d0_d0Errcalc);
                                                    ++ddata_trkd0Errcalc;
                                                }

                                                if (abs(fdata_d0_d0Err) < d0_d0Err_cut && abs(fdata_dz_dzErr) < dz_dzErr_cut)
                                                {
                                                    data_sigmapt_pt->Fill(fdata_ptErr_pt);
                                                    ++ddata_trkptErr;
                                                }

                                                if((abs(fdata_ptErr_pt) < ptErr_pt_cut && abs(fdata_d0_d0Err) < d0_d0Err_cut) && abs(fdata_dz_dzErr) < dz_dzErr_cut)
                                                {
                                                    data_phi_histo->Fill(data_vec.Phi());
                                                    data_validhits->Fill((*nvec_data_validhits)[j]);
                                                    data_trackschi2->Fill((*fvec_data_trackschi2n)[j]);
                                                    ++ddata_trk;
                                                }
                                            }

                                            if(abs(data_vec.Eta()) <= eta_cut)
                                            {
                                                if ((abs(fdata_ptErr_pt) < ptErr_pt_cut && abs(fdata_d0_d0Err) < d0_d0Err_cut) && abs(fdata_dz_dzErr) < dz_dzErr_cut)
                                                {
                                                    data_pt_histo->Fill(data_vec.Pt());
                                                    ++ddata_trkpt;
                                                }
                                            }

                                            if(abs(data_vec.Pt()) >= pt_cut)
                                            {
                                                if ((abs(fdata_ptErr_pt) < ptErr_pt_cut && abs(fdata_d0_d0Err) < d0_d0Err_cut) && abs(fdata_dz_dzErr) < dz_dzErr_cut)
                                                {
                                                    data_eta_histo->Fill(data_vec.Eta());
                                                    ++ddata_trketa;
                                                }
                                            }
                                            ++ndata_HighPurity;
                                        }
                                    }

                                    data_multiplicity->Fill(ddata_trk);
                                    ndata_trk_normalized += ddata_trk;
                                }


                            }


                        }

                    }

                }

            }

        }

    //-----------------------------output data histo------------------------------

    TFile data_histo("data_histo.root", "Recreate");
    gPad->SetLogy();
    data_multiplicity->Write();

    gPad->SetLogy();
    data_dz_sigmadzcalc->Scale(1/ddata_trkdzErrcalc);
    data_dz_sigmadzcalc->Write();
    data_dz_sigmadzcalc->Draw();

    gPad->SetLogy();
    data_d0_sigmad0calc->Scale(1/ddata_trkd0Errcalc);
    data_d0_sigmad0calc->Write();
    data_d0_sigmad0calc->Draw();

    TCanvas *cdata_validhits = new TCanvas ("valid hits", "normalized data valid hits");
    cdata_validhits->cd();
    gPad->SetLogy();
    data_validhits->Scale(1/ddata_evt);
    data_validhits->Draw();
    data_validhits->Write();

    TCanvas *cdata_chi2n = new TCanvas ("chi2n", "normalized data chi2n");
    TLegend *data_leg_validhits = new TLegend (0.1, 0.6, 0.4, 0.9);
    cdata_chi2n->cd();
    gPad->SetLogy();
    data_leg_validhits->AddEntry(data_trackschi2, "Reconstructed chi2", "L");
    data_trackschi2->Scale(1/ddata_evt);
    data_trackschi2->Draw();
    data_trackschi2->Write();

    data_canvas->Divide(3,2);

    data_canvas->cd(1);
    gPad->SetLogy();
    data_pt_histo->Scale(1/ddata_evt);
    data_pt_histo->Draw();
    data_pt_histo->Write();

    data_canvas->cd(2);
    gPad->SetLogy();
    data_eta_histo->Scale(1/ddata_evt);
    data_eta_histo->Draw();
    data_eta_histo->Write();

    data_canvas->cd(3);
    gPad->SetLogy();
    data_phi_histo->Scale(1/ddata_evt);
    data_phi_histo->Draw();
    data_phi_histo->Write();

    data_canvas->cd(4);
    gPad->SetLogy();
    data_dz_sigmadz->Scale(1/ddata_trkdzErr);
    data_dz_sigmadz->GetXaxis()->SetTitle("d_{z}/#sigma_{z}");
    data_dz_sigmadz->GetYaxis()->SetTitle("Fraction of Tracks");
    data_dz_sigmadz->SetLabelSize(0.01, "X");
    data_dz_sigmadz->SetLabelSize(0.01, "Y");
    data_dz_sigmadz->SetMarkerStyle(7);
    data_dz_sigmadz->Draw();
    data_dz_sigmadz->Write();

    data_canvas->cd(5);
    gPad->SetLogy();
    data_d0_sigmad0->Scale(1/ddata_trkd0Err);
    data_d0_sigmad0->GetXaxis()->SetTitle("d_{0}/#sigma_{xy}");
    data_d0_sigmad0->GetYaxis()->SetTitle("Fraction of Tracks");
    data_d0_sigmad0->Draw();
    data_d0_sigmad0->Write();

    data_canvas->cd(6);
    gPad->SetLogy();
    data_sigmapt_pt->Scale(1/ddata_trkpt);
    data_sigmapt_pt->GetXaxis()->SetTitle("#sigma_{p_{T}}/p_{T}");
    data_sigmapt_pt->GetYaxis()->SetTitle("Fraction of Tracks");
    data_sigmapt_pt->SetMarkerStyle(5);
    data_sigmapt_pt->Draw();
    data_sigmapt_pt->Write();

    data_canvas->Write();
    data_canvas->Update();
    data_histo.Write();
    //reco_canvas->SaveAs("data plots.pdf");
    cout << "End of data loop" << endl;
    }

    else if (ndata_MC == 1 || ndata_MC == 2)
    {
        //---------------------loop for MC-------------------------------------------------------
        TFile *myMCFile = TFile::Open("treesCUETP8M1_66.root", "READ");
        TTree *MCtree = (TTree*)myMCFile->Get("UETree/data");

        //for this dataset we want lumisection of 90 and above
        //data from branch 62
        int nMC_lumi_section;
        MCtree->SetBranchAddress("lumi", &nMC_lumi_section);

        //data from Branch 60
        int nMC_ZeroBias;
        MCtree->SetBranchAddress("trgZeroBias",&nMC_ZeroBias);

        //---------------------Variables for MC tracks-------------------------------------------

        vector<float> *fvec_MC_vtxz = 0;
        MCtree->SetBranchAddress("vtxz",&fvec_MC_vtxz);
        float fvec_MC_vtxz_size;

        vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *MC_tracks = 0;
        MCtree->SetBranchAddress("recoTracksp4", &MC_tracks);

        //vector<TLorentzVector> *tracks = new vector<TLorentzVector>;

        vector<float> *fvec_MC_vtxzBS = 0;
        MCtree->SetBranchAddress("vtxzBS", &fvec_MC_vtxzBS);
        float fMC_BeamSize;

        vector<float> *fvec_MC_dz = 0;
        MCtree->SetBranchAddress ("recoTracksdz", &fvec_MC_dz);

        vector<float> *fvec_MC_dzErr = 0;
        MCtree->SetBranchAddress ("recoTracksdzErr", &fvec_MC_dzErr);

        vector<float> *fvec_MC_d0 = 0;
        MCtree->SetBranchAddress ("recoTracksd0", &fvec_MC_d0);

        vector<float> *fvec_MC_d0Err = 0;
        MCtree->SetBranchAddress ("recoTracksd0Err", &fvec_MC_d0Err);

        vector<float> *fvec_MC_ptErr = 0;
        MCtree->SetBranchAddress("recoTracksptErr", &fvec_MC_ptErr);

        vector<int> *nvec_MC_HighPurity = 0;
        MCtree->SetBranchAddress ("recoTrackshighPurity", &nvec_MC_HighPurity);

        vector<float> *fvec_MC_vtxy = 0; //require initialisation to 0 to avoid crash
        MCtree->SetBranchAddress("vtxy",&fvec_MC_vtxy);

        vector<float> *fvec_MC_vtxx = 0; //require initialisation to 0 to avoid crash
        MCtree->SetBranchAddress("vtxx",&fvec_MC_vtxx);

        vector<float> *fvec_MC_trackschi2n = 0;
        MCtree->SetBranchAddress("recoTrackschi2n", &fvec_MC_trackschi2n);

        vector<int> *nvec_MC_validhits = 0;
        MCtree->SetBranchAddress("recoTracksnValidHits", &nvec_MC_validhits);

        vector<int> *nvec_MC_vtxndof = 0;
        MCtree->SetBranchAddress("vtxndof", &nvec_MC_vtxndof);

        int nMC_trk = 0;
        int nMC_Multi = 0;
        int nMC_HighPurity = 0;
        int nMC_Evt = 0;
        int nMC_totaltrk;
        float fMC_dz_dzErr, fMC_d0_d0Err, fMC_ptErr_pt, fMC_dz_dzErrcalc, fMC_d0_d0Errcalc;
        float fMC_dzleaf, fMC_dzcalc, fMC_dzdiff;
        float fMC_d0leaf, fMC_d0calc, fMC_d0diff;
        double dMC_evt = 0;
        double dMC_trk = 0;
        double dMC_trkd0Err = 0;
        double dMC_trkd0Errcalc = 0;
        double dMC_trkdzErr = 0;
        double dMC_trkdzErrcalc = 0;
        double dMC_trkptErr = 0;
        double dMC_trkptErrcalc = 0;

    //---------------------Histograms for MC tracks-------------------------------------------

        TCanvas *MC_canvas = new TCanvas;
        TH1F *MC_pt_histo = new TH1F ("MC p_{T}", "MC p{T}", 80, 0, 4);
        TH1F *MC_eta_histo = new TH1F ("MC #eta", "MC #eta", 60, -3, 3);
        TH1F *MC_phi_histo = new TH1F ("MC #phi", "MC #phi", 80, -4, 4);
        TH1F *MC_multiplicity = new TH1F ("MC_multiplicity", "MC_Multiplicity", 200, 0, 200);
        TH1F *MC_dz_sigmadz = new TH1F ("MC_dz_sigmadz", "MC d_{z}/#sigma_{z}", 160, -20, 20);
        TH1F *MC_dz_sigmadzcalc = new TH1F ("MC_dz_sigmadzcalc", "MC d_{z}/#sigma_{z} calc", 160, -20, 20);
        TH1F *MC_d0_sigmad0 = new TH1F ("MC_d0_sigmad0", "MC d_{0}/#sigma{xy}", 160, -20, 20);
        TH1F *MC_d0_sigmad0calc = new TH1F ("MC_d0_sigmad0calc", "MC d_{0}/#sigma_{xy} calc", 160, -20, 20);
        TH1F *MC_sigmapt_pt = new TH1F ("MC_sigmapt_pt", "MC #sigma_p_{T}/p_{T}", 20, 0, 0.2);
        //TH1F *event_histo = new TH1F ("reco_evt", "reco_evt", 100, 0, 200);
        TH1F *MC_lumi_histo = new TH1F ("MC lumi_section", "MC lumi_section", 160, 80, 230);
        TH1F *MC_validhits = new TH1F ("MC_validhits", "MC_validhits", 100, 0, 50);
        TH1F *MC_trackschi2 = new TH1F ("MC_trackschi2", "MC tracks #chi^{2}", 80, 0, 40);


        int nMC_totalEvt = (Int_t)MCtree->GetEntries();
        cout << "MC Tree Entries " << nMC_totalEvt << endl;

        for(Int_t i = 0; i != nMC_totalEvt; ++i)
        {
            //if (i%10000 == 0)
                //cout<< "Processing event number: " << i << endl; //1 TTree entry = 1 event
            cout << "Entry " << i << endl;
            MCtree -> GetEntry(i);

            //selects for events from lumisection >= 90
            if (nMC_lumi_section >= lumi_cut)
            {
                cout << "MC lumisection is " << nMC_lumi_section << endl;

                //we select only events that are triggered by ZeroBias Trigger. "True" in ZeroBias evaluate to 1
                if (nMC_ZeroBias == 1)
                {
                    fvec_MC_vtxz_size = fvec_MC_vtxz->size();
                    cout << "MC number of vertex is " << fvec_MC_vtxz_size << endl;

                    if (fvec_MC_vtxz_size == vtxz_number)
                    {
                        //looping over vertices
                        for (int k = 0; k != fvec_MC_vtxz_size; ++k)
                        {
                            fMC_BeamSize = fabs((*fvec_MC_vtxz)[k] - (*fvec_MC_vtxzBS)[0]);
                            cout << "Beam Size is " << fMC_BeamSize << endl;

                            if (fMC_BeamSize <= vtxz_size)
                            {
                                nMC_totaltrk = MC_tracks->size();
                                cout << "Total number of MC tracks is " << nMC_totaltrk << endl;

                                //fill with the lumi sections that meet the above event-level cuts
                                MC_lumi_histo->Fill(nMC_lumi_section);
                                ++dMC_evt;

                                if ((*nvec_MC_vtxndof)[0] > 4)
                                {
                                   //looping over tracks
                                    for (int j = 0; j != nMC_totaltrk; ++j)
                                    {
                                        XYZTVector MC_vec = (*MC_tracks)[j];
                                        fMC_dzleaf = (*fvec_MC_dz)[j];
                                        fMC_dzcalc = ((MC_vec.Z())-(*fvec_MC_vtxz)[k])-((((MC_vec.X())-(*fvec_MC_vtxx)[k])*(MC_vec.px())+((MC_vec.Y())-(*fvec_MC_vtxy)[k])*(MC_vec.py()))/(MC_vec.Pt())*(MC_vec.pz()/MC_vec.Pt()));
                                        fMC_dzdiff = ((MC_vec.Z())-(*fvec_MC_vtxz)[k]);
                                        //tr_d0= (- (tr_x-vtx_x)  track.py() + (tr_y-vtx_y)  track.px() ) / track.pt()
                                        fMC_d0leaf = (*fvec_MC_d0)[j];
                                        fMC_d0calc = ((-(MC_vec.X() - (*fvec_MC_vtxx)[k])*MC_vec.Py()) + ((MC_vec.Y()-(*fvec_MC_vtxy)[k])*MC_vec.Px()))/MC_vec.Pt();
                                        dMC_trk = 0;

                                         if ((*nvec_MC_HighPurity)[j] == 1)
                                         {
                                            if (abs (MC_vec.Eta()) <= eta_cut && MC_vec.Pt() >= pt_cut)
                                            {
                                                fMC_dz_dzErr = ((*fvec_MC_dz)[j]) / ((*fvec_MC_dzErr)[j]);
                                                fMC_dz_dzErrcalc = fMC_dzcalc / ((*fvec_MC_dzErr)[j]);
                                                fMC_d0_d0Err = ((*fvec_MC_d0)[j]) / ((*fvec_MC_d0Err)[j]);
                                                fMC_d0_d0Errcalc = fMC_d0calc / ((*fvec_MC_d0Err)[j]);
                                                fMC_ptErr_pt = ((*fvec_MC_ptErr)[j] / (MC_vec.Pt()));


                                                //fphi_eff = reco_vec.Phi()/gen_vec.Phi();
                                                //feta_eff = reco_vec.Eta()/gen_vec.Eta();
                                                //fpt_eff = reco_vec.Pt()/gen_vec.Pt();

                                                //pt_efficiency->Fill(fpt_eff);
                                                //eta_efficiency->Fill(feta_eff);
                                                //phi_efficiency->Fill(fphi_eff);

                                                if ((abs(fMC_ptErr_pt)) < ptErr_pt_cut && abs(fMC_d0_d0Err) < d0_d0Err_cut)
                                                {
                                                    MC_dz_sigmadz->Fill(fMC_dz_dzErr);
                                                    ++dMC_trkdzErr;
                                                }

                                                if ((abs(fMC_ptErr_pt)) < ptErr_pt_cut && abs(fMC_d0_d0Errcalc) < d0_d0Err_cut)
                                                {
                                                    MC_dz_sigmadzcalc->Fill(fMC_dz_dzErrcalc);
                                                    ++dMC_trkdzErrcalc;
                                                }

                                                if (abs(fMC_ptErr_pt) < ptErr_pt_cut && abs(fMC_dz_dzErr) < dz_dzErr_cut)
                                                {
                                                    MC_d0_sigmad0->Fill(fMC_d0_d0Err);
                                                    ++dMC_trkd0Err;
                                                }

                                                if (abs(fMC_ptErr_pt) < ptErr_pt_cut && abs(fMC_dz_dzErrcalc) < dz_dzErr_cut)
                                                {
                                                    MC_d0_sigmad0calc->Fill(fMC_d0_d0Errcalc);
                                                    ++dMC_trkd0Errcalc;
                                                }

                                                if (abs(fMC_d0_d0Err) < d0_d0Err_cut && abs(fMC_dz_dzErr) < dz_dzErr_cut)
                                                {
                                                    MC_sigmapt_pt->Fill(fMC_ptErr_pt);
                                                    ++dMC_trkptErr;
                                                }

                                                if((abs(fMC_ptErr_pt) < ptErr_pt_cut && abs(fMC_d0_d0Err) < d0_d0Err_cut) && abs(fMC_dz_dzErr) < dz_dzErr_cut)
                                                {
                                                    MC_phi_histo->Fill(MC_vec.Phi());
                                                    MC_validhits->Fill((*nvec_MC_validhits)[j]);
                                                    MC_trackschi2->Fill((*fvec_MC_trackschi2n)[j]);
                                                    ++dMC_trk;
                                                }
                                            }

                                            if (MC_vec.Pt() >= pt_cut)
                                            {
                                                if((abs(fMC_ptErr_pt) < ptErr_pt_cut && abs(fMC_d0_d0Err) < d0_d0Err_cut) && abs(fMC_dz_dzErr) < dz_dzErr_cut)
                                                {
                                                    MC_eta_histo->Fill(MC_vec.Eta());
                                                }
                                            }

                                            if (abs (MC_vec.Eta()) <= eta_cut)
                                            {
                                                if ((abs(fMC_ptErr_pt) < ptErr_pt_cut && abs(fMC_d0_d0Err) < d0_d0Err_cut) && abs(fMC_dz_dzErr) < dz_dzErr_cut)
                                                {
                                                    MC_pt_histo->Fill(MC_vec.Pt());
                                                }
                                            }

                                            ++nMC_HighPurity;
                                        }
                                    }
                                    MC_multiplicity->Fill(dMC_trk);
                                }


                            }
                        }
                    }
                }
            }
        }

        //-----------------------------output MC histo------------------------------

        TFile MC_histo("MC_histo.root", "Recreate");
        gPad->SetLogy();
        MC_multiplicity->Write();

        gPad->SetLogy();
        MC_dz_sigmadzcalc->Scale(1/dMC_trkdzErrcalc);
        MC_dz_sigmadzcalc->Write();
        MC_dz_sigmadzcalc->Draw();

        gPad->SetLogy();
        MC_d0_sigmad0calc->Scale(1/dMC_trkd0Errcalc);
        MC_d0_sigmad0calc->Write();
        MC_d0_sigmad0calc->Draw();

        TCanvas *cMC_validhits = new TCanvas ("MC validhits", "Normalized MC validhits");
        cMC_validhits->cd(1);
        gPad->SetLogy();
        MC_validhits->Scale(1/dMC_evt);
        MC_validhits->Draw();
        MC_validhits->Write();

        TCanvas *cMC_trackschi2n = new TCanvas ("MC trackschi2n", "Normalized MC tracks2n");
        cMC_trackschi2n->cd(1);
        gPad->SetLogy();
        MC_trackschi2->Scale(1/dMC_evt);
        MC_trackschi2->Draw();
        MC_trackschi2->Write();

        MC_canvas->Divide(3,2);

        MC_canvas->cd(1);
        gPad->SetLogy();
        MC_pt_histo->Scale(1/dMC_evt);
        MC_pt_histo->Draw();
        MC_pt_histo->Write();
        //MC_multiplicity->SetLineColor(kBlue+2);
        //MC_multiplicity->DrawNormalized("", 1);

        MC_canvas->cd(2);
        gPad->SetLogy();
        MC_eta_histo->SetLineColor(kBlue-4);
        //MC_eta_histo->Scale(1/dMC_evt, "SAME");
        MC_eta_histo->Scale(1/dMC_evt);
        MC_eta_histo->Draw();
        MC_eta_histo->Write();

        MC_canvas->cd(3);
        gPad->SetLogy();
        //MC_phi_histo->SetLineColor(kGreen);
        MC_phi_histo->Scale(1/dMC_evt);
        MC_phi_histo->Draw();
        MC_phi_histo->Write();

        MC_canvas->cd(4);
        gPad->SetLogy();
        MC_dz_sigmadz->Scale(1/dMC_trkdzErr);
        MC_dz_sigmadz->Draw();
        MC_dz_sigmadz->Write();

        MC_canvas->cd(5);
        gPad->SetLogy();
        MC_d0_sigmad0->Scale(1/dMC_trkd0Err);
        MC_d0_sigmad0->Draw();
        MC_d0_sigmad0->Write();

        MC_canvas->cd(6);
        gPad->SetLogy();
        MC_sigmapt_pt->Scale(1/dMC_trkptErr);
        MC_sigmapt_pt->Draw();
        MC_sigmapt_pt->Write();

        MC_histo.Write();

        MC_canvas->Write();
        MC_canvas->Update();
        //MC_canvas->SaveAs("MC Plots.pdf");
        cout << "End of MC loop" << endl;
    }

    else
    {
        cout << "You have entered an invalid number. Please try again." << endl;
    }

/*cout << "The number of tracks for data dz_sigmadz is " << dreco_trkdzErr << endl;
cout << "The number of tracks for data d0_sigmad0 is " << dreco_trkd0Err << endl;
cout << "The number of tracks for data sigmapt_pt is " << dreco_trkptErr << endl;
cout << "The number of tracks for MC dz_sigmadz is " << dMC_trkdzErr << endl;
cout << "The number of tracks for MC d0_sigmad0 is " << dMC_trkd0Err<< endl;
cout << "The number of tracks for MC sigmapt_pt is " << dMC_trkptErr << endl;
cout << "Total number of data events is " << ddata_evt << endl;
cout << "Total number of MC events is " << dMC_evt << endl;*/



//---------------------Efficiency-------------------------------------

    /*TH1F *pt_efficiency = new TH1F("pt_efficiency", "pt efficiency", 100, 0, 200);
    TH1F *eta_efficiency = new TH1F("eta_efficiency", "eta efficiency", 100, 0, 200);
    TH1F *phi_efficiency = new TH1F("phi_efficiency", "phi efficiency", 100, 0, 200);
    float fpt_eff, feta_eff, fphi_eff;*/

}
