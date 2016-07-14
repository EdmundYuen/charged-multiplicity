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

using namespace std;
using namespace ROOT::Math;

void trees1_CUETP8M1_66()
{
   //gROOT->ProcessLine(".L Loader.C+");

    // attach "UETree/data" of tree1.root as the main root file for this program
    TFile *myFile = TFile::Open("tree1.root", "READ");
    TTree *tree = (TTree*)myFile->Get("UETree/data");

    //for this dataset we want lumisection of 90 and above
    //data from branch 62
    int nlumi_section;
    tree->SetBranchAddress("lumi", &nlumi_section);

    int nZeroBias; //data from Branch 60
    tree->SetBranchAddress("trgZeroBias",&nZeroBias);

//---------------------Variables for reco tracks-------------------------------------------

    vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *reco_tracks = 0;
    tree->SetBranchAddress("recoTracksp4", &reco_tracks);
    //vector<TLorentzVector> *tracks = new vector<TLorentzVector>;

    //variables to check if there is only 1 vertex
    vector<float> *fvec_reco_vtxz = 0; //require initialisation to 0 to avoid crash
    tree->SetBranchAddress("vtxz",&fvec_reco_vtxz);
    float fvec_reco_vtxz_size;

    vector<float> *fvec_reco_vtxzBS = 0; //from branch 54
    tree->SetBranchAddress("vtxzBS", &fvec_reco_vtxzBS);
    float n_reco_BeamSize;

    vector<float> *fvec_dz = 0;
    tree->SetBranchAddress ("recoTracksdz", &fvec_dz);

    vector<float> *fvec_dzErr = 0;
    tree->SetBranchAddress ("recoTracksdzErr", &fvec_dzErr);

    vector<float> *fvec_d0 = 0;
    tree->SetBranchAddress ("recoTracksd0", &fvec_d0);

    vector<float> *fvec_d0Err = 0;
    tree->SetBranchAddress ("recoTracksd0Err", &fvec_d0Err);

    vector<float> *fvec_ptErr = 0;
    tree->SetBranchAddress("recoTracksptErr", &fvec_ptErr);

    vector<int> *nvec_HighPurity = 0;
    tree->SetBranchAddress ("recoTrackshighPurity", &nvec_HighPurity);

    vector<float> *fvec_reco_vtxy = 0; //require initialisation to 0 to avoid crash
    tree->SetBranchAddress("vtxy",&fvec_reco_vtxy);

	vector<float> *fvec_reco_vtxx = 0; //require initialisation to 0 to avoid crash
    tree->SetBranchAddress("vtxx",&fvec_reco_vtxx);

    vector<float> *fvec_reco_trackschi2n = 0;
    tree->SetBranchAddress("recoTrackschi2n", &fvec_reco_trackschi2n);
    float freco_chi2n;

    vector<int> *nvec_reco_validhits = 0;
    tree->SetBranchAddress("recoTracksnValidHits", &nvec_reco_validhits);
    int nreco_validhits;

    vector<int> *nvec_reco_ndof = 0;
    tree->SetBranchAddress("vtxndof", &nvec_reco_ndof);


    double dzleaf = 0;
	double dzcalc = 0;
	double dzminus= 0;

	double d0leaf = 0;
	double d0calc = 0;
	double d0minus = 0;

    int nreco_totaltrk;
    double dreco_trk;
    double dreco_trkpterr = 0;
    double dreco_trkd0err = 0;
    double dreco_trkdzerr = 0;
    double dreco_trkpt = 0;
    double dreco_trketa = 0;
    double dreco_trkphi = 0;
    double dreco_evt = 0;

    int nreco_trk_normalized = 0;
    int nreco_HighPurity = 0;
    float dz_dzErr, d0_d0Err, ptErr_pt;

//---------------------Histograms for reco Tracks-------------------------------------------

    TCanvas *reco_canvas = new TCanvas;
    //TCanvas *gen_canvas = new TCanvas;
    //TH1F *event_histo = new TH1F ("reco_evt", "reco_evt", 100, 0, 200);
    TH1F *reco_pt_histo = new TH1F ("reco_pt", "Normalized_reco_pT", 80, 0, 4);
    TH1F *reco_eta_histo = new TH1F ("reco_eta", "Normalized_reco_Eta", 60, -3, 3);
    TH1F *reco_phi_histo = new TH1F ("reco_phi", "Normalized_reco_Phi", 80, -4, 4);
    TH1F *lumi_histo = new TH1F ("lumi_section", "lumi_section", 160, 80, 230);
    TH1F *reco_multiplicity = new TH1F ("reco_multiplicity", "Normalized_reco_Multiplicity", 200, 0, 200);
    TH1F *reco_dz_sigmadz = new TH1F ("dz_sigmadz", "reco_dz_sigmadz", 100, -10, 10);
    TH1F *reco_d0_sigmad0 = new TH1F ("d0_sigmad0", "reco_d0_sigmad0", 200, -10, 10);
    TH1F *reco_sigmapt_pt = new TH1F ("sigmapt_pt", "reco_sigmapt_pt", 80, 0, 0.8);
    //TH1F *normalized_multiplicity_histo = new TH1F ("normalized_multiplicity", "normalized_multiplicity", 200, 0, 200);
    TH1F *reco_validhits = new TH1F ("reco_validhits", "reco_validhits", 100, 0, 50);
    TH1F *reco_trackschi2 = new TH1F ("reco_trackschi2", "reco_trackschi2", 80, 0, 40);

//---------------------Efficiency-------------------------------------

    /*TH1F *pt_efficiency = new TH1F("pt_efficiency", "pt efficiency", 100, 0, 200);
    TH1F *eta_efficiency = new TH1F("eta_efficiency", "eta efficiency", 100, 0, 200);
    TH1F *phi_efficiency = new TH1F("phi_efficiency", "phi efficiency", 100, 0, 200);
    float fpt_eff, feta_eff, fphi_eff;*/

//---------------------Cuts-------------------------------------------

    const int lumi_cut = 90;
    const double eta_cut = 2.4;
    const double pt_cut = 0.5;
    const double vtxz_number = 1.;
    const double vtxz_size = 10;
    const float d0_d0Err_cut = 3;
    const float dz_dzErr_cut = 3;
    const float ptErr_pt_cut = 0.05;

    //select events with only 1 vertex using information from Branch 46/47/48. Here we use Branch 48 (z-axis)
    Int_t nEvt = (Int_t)tree->GetEntries();
    cout << "Tree Entries " << nEvt << endl;


    //start loop over events
    for(Int_t i = 0; i != nEvt; ++i)
    {
        //if (i%10000 == 0)
            //cout<< "Processing event number: " << i << endl; //1 TTree entry = 1 event
        cout << "Entry " << i << endl;
        tree -> GetEntry(i);

        //selects for events from lumisection >= 90
        if (nlumi_section >= lumi_cut)
        {
            cout << "lumisection is " << nlumi_section << endl;


            //we select only events that are triggered by ZeroBias Trigger. "True" in ZeroBias evaluate to 1
            if (nZeroBias == 1)
            {
                fvec_reco_vtxz_size = fvec_reco_vtxz->size();
                //vtxz_plot->Fill(fvecVtxz_size);
                /*ntrk_MC = MC_tracks->size();
                nMC_Multi = 0;

//---------------------Loops for MC tracks-------------------------------------------

                //looping over generated tracks
                for(int q = 0; q != ntrk_MC; ++q)
                {
                    XYZTVector gen_vec = (*gen_tracks)[q];

                    if (abs (gen_vec.Eta()) <= eta_cut && gen_vec.Pt() >= pt_cut)
                    {
                        MC_phi_histo->Fill(gen_vec.Phi());
                        MC_pt_histo->Fill(gen_vec.Pt());
                        MC_eta_histo->Fill(gen_vec.Eta());
                        ++nMC_Multi;
                    }
                }
                gen_multiplicity->Fill(ngen_Multi);*/


//---------------------Loops for reco tracks-------------------------------------------

                if (fvec_reco_vtxz_size == vtxz_number)
                {
                    //looping over vertices
                    for (int k = 0; k != fvec_reco_vtxz_size; ++k)
                    {
                        n_reco_BeamSize = fabs((*fvec_reco_vtxz)[k] - (*fvec_reco_vtxzBS)[0]);
                        cout << "Beam Size is " << n_reco_BeamSize << endl;

                        if (n_reco_BeamSize <= vtxz_size)
                        {
                            nreco_totaltrk = reco_tracks->size();

                            //fill with the lumi sections that meet the above event-level cuts
                            lumi_histo->Fill(nlumi_section);
                            ++dreco_evt;

                            //looping over tracks
                            for (int j = 0; j != nreco_totaltrk; ++j)
                            {
                                XYZTVector reco_vec = (*reco_tracks)[j];

                                if (j == 9)
								{
									//dzleaf = (*fvec_dz)[j];
									dzcalc = ((reco_vec.Z())-(*fvec_reco_vtxz)[k])-((((reco_vec.X())-(*fvec_reco_vtxx)[k])*(reco_vec.px())+((reco_vec.Y())-(*fvec_reco_vtxy)[k])*(reco_vec.py()))/(reco_vec.Pt())*(reco_vec.pz()/reco_vec.Pt()));
									//dzminus = ((vec.Z())-(*fvec_vtxz)[k]);
                                    //tr_d0= (- (tr_x-vtx_x)  track.py() + (tr_y-vtx_y)  track.px() ) / track.pt()
									d0leaf = (*fvec_d0)[j];
									//d0calc = ((-(reco_vec.X() - (*fvec_reco_vtxx)[k])*reco_vec.Py()) + ((reco_vec.Y()-(*fvec_reco_vtxy)[k])*reco_vec.Px()))/reco_vec.Pt();
                                    d0calc = ((reco_vec.Y() - (*fvec_reco_vtxy)[k])*reco_vec.Px() - (reco_vec.X() - (*fvec_reco_vtxx)[k])*reco_vec.Py())/reco_vec.Pt();
								}

                                if ((*nvec_HighPurity)[j] == 1)
                                {
                                    if (abs (reco_vec.Eta()) <= eta_cut && reco_vec.Pt() >= pt_cut)
                                    {
                                        dz_dzErr = ((*fvec_dz)[j])/((*fvec_dzErr)[j]);
                                        d0_d0Err = ((*fvec_d0)[j])/((*fvec_d0Err)[j]);
                                        ptErr_pt = ((*fvec_ptErr)[j]/(reco_vec.Pt()));


                                        //fphi_eff = reco_vec.Phi()/gen_vec.Phi();
                                        //feta_eff = reco_vec.Eta()/gen_vec.Eta();
                                        //fpt_eff = reco_vec.Pt()/gen_vec.Pt();
                                        //reco_dz_sigmadz->Fill(dz_dzErr);
                                        //reco_d0_sigmad0->Fill(d0_d0Err);
                                        //reco_sigmapt_pt->Fill(ptErr_pt);

                                        if (abs(ptErr_pt) < ptErr_pt_cut && abs(d0_d0Err) < d0_d0Err_cut)
                                        {
                                            reco_dz_sigmadz->Fill(dz_dzErr);
                                            ++dreco_trkdzerr;
                                        }

                                        if (abs(ptErr_pt) < ptErr_pt_cut && abs(dz_dzErr) < dz_dzErr_cut)
                                        {
                                            reco_d0_sigmad0->Fill(d0_d0Err);
                                            ++dreco_trkd0err;
                                        }

                                        if (abs(d0_d0Err) < d0_d0Err_cut && abs(dz_dzErr) < dz_dzErr_cut)
                                        {
                                            reco_sigmapt_pt->Fill(ptErr_pt);
                                            ++dreco_trkpterr;
                                        }

                                        if((abs(ptErr_pt) < ptErr_pt_cut && abs(d0_d0Err) < d0_d0Err_cut) && abs(dz_dzErr) < dz_dzErr_cut)
                                        {
                                            reco_phi_histo->Fill(reco_vec.Phi());
                                            reco_validhits->Fill((*nvec_reco_validhits)[j]);
                                            reco_trackschi2->Fill((*fvec_reco_trackschi2n)[j]);
                                            ++dreco_trk;
                                        }
                                    }

                                    if(abs(reco_vec.Eta()) <= eta_cut)
                                    {
                                        if ((abs(ptErr_pt) < ptErr_pt_cut && abs(d0_d0Err) < d0_d0Err_cut) && abs(dz_dzErr) < dz_dzErr_cut)
                                        {
                                            reco_pt_histo->Fill(reco_vec.Pt());
                                            ++dreco_trkpt;
                                        }
                                    }

                                    if(abs(reco_vec.Pt()) >= pt_cut)
                                    {
                                        if ((abs(ptErr_pt) < ptErr_pt_cut && abs(d0_d0Err) < d0_d0Err_cut) && abs(dz_dzErr) < dz_dzErr_cut)
                                        {
                                            reco_eta_histo->Fill(reco_vec.Eta());
                                            ++dreco_trketa;
                                        }
                                    }
                                    ++nreco_HighPurity;
                                }



                                /*if (abs (vec.Eta()) <= eta_cut && vec.Pt() >= pt_cut)
                                {
                                    pt_histo->Fill(vec.Pt());
                                }

                                if (vec.Pt() >= pt_cut)
                                {
                                    eta_histo->Fill(vec.Eta());
                                }*/


                            }

                            reco_multiplicity->Fill(dreco_trk);
                            nreco_trk_normalized += dreco_trk;
                        }


                    }

                }

            }

        }

    }

    //---------------------loop for MC-------------------------------------------------------
    TFile *myMCFile = TFile::Open("treesCUETP8M1_66.root", "READ");
    TTree *MCtree = (TTree*)myMCFile->Get("UETree/data");

    //for this dataset we want lumisection of 90 and above
    //data from branch 62
    nlumi_section = 0;
    MCtree->SetBranchAddress("lumi", &nlumi_section);

    //data from Branch 60
    nZeroBias = 0;
    MCtree->SetBranchAddress("trgZeroBias",&nZeroBias);

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

    int nMC_trk = 0;
    int nMC_Multi = 0;
    int nMC_HighPurity = 0;
    nEvt = 0;
    double dMC_evt = 0;
    double dMC_trk = 0;
    double dMC_trkd0err = 0;
    double dMC_trkdzerr = 0;
    double dMC_trkpterr = 0;
    int nMC_totaltrk;


//---------------------Histograms for MC tracks-------------------------------------------

    TCanvas *MC_canvas = new TCanvas;
    TH1F *MC_pt_histo = new TH1F ("MC_pt", "MC_pT", 80, 0, 4);
    TH1F *MC_eta_histo = new TH1F ("MC_eta", "MC_Eta", 60, -3, 3);
    TH1F *MC_phi_histo = new TH1F ("MC_phi", "MC_Phi", 80, -4, 4);
    TH1F *MC_multiplicity = new TH1F ("MC_multiplicity", "MC_Multiplicity", 200, 0, 200);
    TH1F *MC_dz_sigmadz = new TH1F ("dz_sigmadz", "MC_dz_sigmadz", 100, -10, 10);
    TH1F *MC_d0_sigmad0 = new TH1F ("d0_sigmad0", "MC_d0_sigmad0", 200, -10, 10);
    TH1F *MC_sigmapt_pt = new TH1F ("sigmapt_pt", "MC_sigmapt_pt", 80, 0, 0.8);
    //TH1F *event_histo = new TH1F ("reco_evt", "reco_evt", 100, 0, 200);
    //TH1F *lumi_histo = new TH1F ("lumi_section", "lumi_section", 160, 80, 230);
    TH1F *MC_validhits = new TH1F ("MC_validhits", "MC_validhits", 100, 0, 50);
    TH1F *MC_trackschi2 = new TH1F ("MC_trackschi2", "MC_trackschi2", 80, 0, 40);


    nEvt = (Int_t)MCtree->GetEntries();
    cout << "MC Tree Entries " << nEvt << endl;

    for(Int_t i = 0; i != nEvt; ++i)
    {
        //if (i%10000 == 0)
            //cout<< "Processing event number: " << i << endl; //1 TTree entry = 1 event
        cout << "Entry " << i << endl;
        MCtree -> GetEntry(i);

        //selects for events from lumisection >= 90
        if (nlumi_section >= lumi_cut)
        {
            cout << "MC lumisection is " << nlumi_section << endl;

            //we select only events that are triggered by ZeroBias Trigger. "True" in ZeroBias evaluate to 1
            if (nZeroBias == 1)
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
                            lumi_histo->Fill(nlumi_section);
                            ++dMC_evt;

                            //looping over tracks
                            for (int j = 0; j != nMC_totaltrk; ++j)
                            {
                                XYZTVector MC_vec = (*MC_tracks)[j];
                                dMC_trk = 0;

                                 if ((*nvec_MC_HighPurity)[j] == 1)
                                 {
                                    if (abs (MC_vec.Eta()) <= eta_cut && MC_vec.Pt() >= pt_cut)
                                    {
                                        dz_dzErr = ((*fvec_MC_dz)[j])/((*fvec_MC_dzErr)[j]);
                                        d0_d0Err = ((*fvec_MC_d0)[j])/((*fvec_MC_d0Err)[j]);
                                        ptErr_pt = ((*fvec_MC_ptErr)[j]/(MC_vec.Pt()));

                                        //fphi_eff = reco_vec.Phi()/gen_vec.Phi();
                                        //feta_eff = reco_vec.Eta()/gen_vec.Eta();
                                        //fpt_eff = reco_vec.Pt()/gen_vec.Pt();

                                        //pt_efficiency->Fill(fpt_eff);
                                        //eta_efficiency->Fill(feta_eff);
                                        //phi_efficiency->Fill(fphi_eff);
                                        MC_phi_histo->Fill(MC_vec.Phi());

                                        if ((abs(ptErr_pt)) < ptErr_pt_cut && abs(d0_d0Err) < d0_d0Err_cut)
                                        {
                                            MC_dz_sigmadz->Fill(dz_dzErr);
                                            ++dMC_trkdzerr;
                                        }

                                        if (abs(ptErr_pt) < ptErr_pt_cut && abs(dz_dzErr) < dz_dzErr_cut)
                                        {
                                            MC_d0_sigmad0->Fill(d0_d0Err);
                                            ++dMC_trkd0err;
                                        }

                                        if (abs(d0_d0Err) < d0_d0Err_cut && abs(dz_dzErr) < dz_dzErr_cut)
                                        {
                                            MC_sigmapt_pt->Fill(ptErr_pt);
                                            ++dMC_trkpterr;
                                        }

                                        if((abs(ptErr_pt) < ptErr_pt_cut && abs(d0_d0Err) < d0_d0Err_cut) && abs(dz_dzErr) < dz_dzErr_cut)
                                        {
                                            MC_phi_histo->Fill(MC_vec.Phi());
                                            MC_validhits->Fill((*nvec_MC_validhits)[j]);
                                            MC_trackschi2->Fill((*fvec_MC_trackschi2n)[j]);
                                            ++dMC_trk;
                                        }
                                    }

                                    if (MC_vec.Pt() >= pt_cut)
                                    {
                                        if((abs(ptErr_pt) < ptErr_pt_cut && abs(d0_d0Err) < d0_d0Err_cut) && abs(dz_dzErr) < dz_dzErr_cut)
                                        {
                                            MC_eta_histo->Fill(MC_vec.Eta());
                                        }
                                    }

                                    if (abs (MC_vec.Eta()) <= eta_cut)
                                    {
                                        if ((abs(ptErr_pt) < ptErr_pt_cut && abs(d0_d0Err) < d0_d0Err_cut) && abs(dz_dzErr) < dz_dzErr_cut)
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

    //cout << "Total number of selected tracks is " << ntrk_normalized << endl;


    //vtxz_plot->Draw();
    //reco_canvas->Update();
    //reco_canvas->SaveAs("vtxz_number.png");
    //delete reco_canvas;
    //delete vtxz_plot;


    //after the loop over all events, draw the resulting plots
    cout << "The number of tracks for data dz_sigmadz is " << dreco_trkdzerr << endl;
    cout << "The number of tracks for data d0_sigmad0 is " << dreco_trkd0err << endl;
    cout << "The number of tracks for data sigmapt_pt is " << dreco_trkpterr << endl;
    cout << "The number of tracks for MC dz_sigmadz is " << dMC_trkdzerr << endl;
    cout << "The number of tracks for MC d0_sigmad0 is " << dMC_trkd0err<< endl;
    cout << "The number of tracks for MC sigmapt_pt is " << dMC_trkpterr << endl;
    cout << "Total number of data events is " << dreco_evt << endl;
    cout << "Total number of MC events is " << dMC_evt << endl;

    TFile histo("histo.root", "recreate");

    //-----------------------------output reco histo------------------------------

    gPad->SetLogy();
    reco_multiplicity->Write();

    TCanvas *creco_validhits = new TCanvas ("valid hits", "normalized valid hits");
    creco_validhits->cd(1);
    gPad->SetLogy();
    reco_validhits->Scale(1/dreco_evt);
    reco_validhits->Draw();
    reco_validhits->Write();

    TCanvas *creco_chi2n = new TCanvas ("chi2n", "normalized chi2n");
    creco_chi2n->cd(1);
    gPad->SetLogy();
    TLegend *reco_leg_validhits = new TLegend (0.1, 0.75, 0.38, 0.9);
    reco_leg_validhits->AddEntry(reco_trackschi2, "Reconstructed chi2", "L");
    reco_trackschi2->Scale(1/dreco_evt);
    reco_trackschi2->Draw();
    reco_trackschi2->Write();

    reco_canvas->Divide(3,2);

    reco_canvas->cd(1);
    gPad->SetLogy();
    reco_pt_histo->Scale(1/dreco_evt);
    reco_pt_histo->Draw();
    reco_pt_histo->Write();

    reco_canvas->cd(2);
    gPad->SetLogy();
    reco_eta_histo->Scale(1/dreco_evt);
    reco_eta_histo->Draw();
    reco_eta_histo->Write();

    reco_canvas->cd(3);
    gPad->SetLogy();
    reco_phi_histo->Scale(1/dreco_evt);
    reco_phi_histo->Draw();
    reco_phi_histo->Write();

    reco_canvas->cd(4);
    gPad->SetLogy();
    reco_dz_sigmadz->Scale(1/dreco_trkdzerr);
    reco_dz_sigmadz->Draw();
    reco_dz_sigmadz->Write();

    reco_canvas->cd(5);
    gPad->SetLogy();
    reco_d0_sigmad0->Scale(1/dreco_trkd0err);
    reco_d0_sigmad0->Draw();
    reco_d0_sigmad0->Write();

    reco_canvas->cd(6);
    gPad->SetLogy();
    reco_sigmapt_pt->Scale(1/dreco_trkpt);
    reco_sigmapt_pt->Draw();
    reco_sigmapt_pt->Write();

    reco_canvas->Write();
    reco_canvas->Update();
    //reco_canvas->SaveAs("data plots.pdf");

    //-----------------------------output MC histo------------------------------

    gPad->SetLogy();
    MC_multiplicity->Write();

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
    MC_dz_sigmadz->Scale(1/dMC_trkdzerr);
    MC_dz_sigmadz->Draw();
    MC_dz_sigmadz->Write();

    MC_canvas->cd(5);
    gPad->SetLogy();
    MC_d0_sigmad0->Scale(1/dMC_trkd0err);
    MC_d0_sigmad0->Draw();
    MC_d0_sigmad0->Write();

    MC_canvas->cd(6);
    gPad->SetLogy();
    MC_sigmapt_pt->Scale(1/dMC_trkpterr);
    MC_sigmapt_pt->Draw();
    MC_sigmapt_pt->Write();

    histo.Write();

    MC_canvas->Write();
    MC_canvas->Update();
    //MC_canvas->SaveAs("MC Plots.pdf");
}
