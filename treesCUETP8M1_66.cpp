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

using namespace std;
using namespace ROOT::Math;

void treesCUETP8M1_66()
{
   //gROOT->ProcessLine(".L Loader.C+");

    // attach "UETree/data" of tree1.root as the main root file for this program
    TFile *myFile = TFile::Open("treesCUETP8M1_66.root", "READ");
    TTree* tree = (TTree*)myFile->Get("UETree/data");

    //for this dataset we want lumisection of 90 and above
    //data from branch 62
    int nlumi_section;
    tree->SetBranchAddress("lumi", &nlumi_section);

    int iZeroBias; //data from Branch 60
    tree->SetBranchAddress("trgZeroBias",&iZeroBias);

//---------------------Variables for reco tracks-------------------------------------------

    //variables to check if there is only 1 vertex
    vector<float> *fVec_reco_Vtxz = 0; //require initialisation to 0 to avoid crash
    tree->SetBranchAddress("vtxz",&fVec_reco_Vtxz);
    float fVec_reco_Vtxz_size;

    vector<float> *fVec_reco_VtxzBS = 0; //from branch 54
    tree->SetBranchAddress("vtxzBS", &fVec_reco_VtxzBS);
    float n_reco_BeamSize;

    vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *reco_tracks = 0;
    tree->SetBranchAddress("recoTracksp4", &reco_tracks);
    //vector<TLorentzVector> *tracks = new vector<TLorentzVector>;

    vector<float> *fVec_dz = 0;
    tree->SetBranchAddress ("recoTracksdz", &fVec_dz);

    vector<float> *fVec_dzErr = 0;
    tree->SetBranchAddress ("recoTracksdzErr", &fVec_dzErr);

    vector<float> *fVec_d0 = 0;
    tree->SetBranchAddress ("recoTracksd0", &fVec_d0);
    vector<float> *fVec_d0Err = 0;
    tree->SetBranchAddress ("recoTracksd0Err", &fVec_d0Err);

    vector<float> *fVec_ptErr = 0;
    tree->SetBranchAddress("recoTracksptErr", &fVec_ptErr);

    vector<int> *nVec_HighPurity = 0;
    tree->SetBranchAddress ("recoTrackshighPurity", &nVec_HighPurity);

    vector<float> *fvec_reco_Vtxy = 0; //require initialisation to 0 to avoid crash
    tree->SetBranchAddress("vtxy",&fvec_reco_Vtxy);

	vector<float> *fvec_reco_Vtxx = 0; //require initialisation to 0 to avoid crash
    tree->SetBranchAddress("vtxx",&fvec_reco_Vtxx);

    double dzleaf = 0;
	double dzcalc = 0;
	double dzminus= 0;

	double d0leaf = 0;
	double d0calc = 0;
	double d0minus = 0;


    //declare variable to hold track number
    int ntrk_reco_total,
    double d_reco_trk;
    double d_reco_evt = 0;
    int ntrk_normalized = 0;
    int nHigh_Purity = 0;
    float dz_dzErr, d0_d0Err, ptErr_pt;

//---------------------Histograms for reco Tracks-------------------------------------------

    TCanvas *reco_canvas = new TCanvas;
    TCanvas *gen_canvas = new TCanvas;
    //TH1F *event_histo = new TH1F ("reco_evt", "reco_evt", 100, 0, 200);
    TH1F *reco_pt_histo = new TH1F ("reco_pt", "Normalized_reco_pT", 100, 0, 100);
    TH1F *reco_eta_histo = new TH1F ("reco_eta", "Normalized_reco_Eta", 25, -3, 3);
    TH1F *reco_phi_histo = new TH1F ("reco_phi", "Normalized_reco_Phi", 100, -4, 4);
    TH1F *lumi_histo = new TH1F ("lumi_section", "lumi_section", 160, 80, 230);
    TH1F *reco_multiplicity = new TH1F ("reco_multiplicity", "Normalized_reco_Multiplicity", 200, 0, 200);
    TH1F *dz_sigmadz = new TH1F ("dz_sigmadz", "Normalized_dz_sigmadz", 160, -20, 20);
    TH1F *d0_sigmad0 = new TH1F ("d0_sigmad0", "Normalized_d0_sigmad0", 160, -20, 20);
    TH1F *sigmapt_pt = new TH1F ("sigmapt_pt", "Normalized_sigmapt_pt", 20, 0, 0.2);
    //TH1F *normalized_multiplicity_histo = new TH1F ("normalized_multiplicity", "normalized_multiplicity", 200, 0, 200);

//---------------------Variables for gen tracks-------------------------------------------

    float *f_gen_Vtxz = 0;
    tree->SetBranchAddress("simvtxz",&f_gen_Vtxz);

    vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *gen_tracks = 0;
    tree->SetBranchAddress("genParticlesp4", &gen_tracks);
    //vector<TLorentzVector> *tracks = new vector<TLorentzVector>;

    int ntrk_gen, ngen_Multi;

//---------------------Histograms for gen tracks-------------------------------------------

    TH1F *gen_pt_histo = new TH1F ("gen_pt", "Normalized_gen_pT", 100, 0, 100);
    TH1F *gen_eta_histo = new TH1F ("gen_eta", "Normalized_gen_Eta", 25, -3, 3);
    TH1F *gen_phi_histo = new TH1F ("gen_phi", "Normalized_gen_Phi", 100, -4, 4);
    TH1F *gen_multiplicity = new TH1F ("gen_multiplicity", "Normalized_Multiplicity", 200, 0, 200);

//---------------------Efficiency-------------------------------------

    TH1F *pt_efficiency = new TH1F("pt_efficiency", "pt efficiency", 100, 0, 200);
    TH1F *eta_efficiency = new TH1F("eta_efficiency", "eta efficiency", 100, 0, 200);
    TH1F *phi_efficiency = new TH1F("phi_efficiency", "phi efficiency", 100, 0, 200);
    float fpt_eff, feta_eff, fphi_eff;

//---------------------Cuts-------------------------------------------

    const int lumi_cut = 90;
    const double eta_cut = 2.4;
    const double pt_cut = 0.5;
    const double vtxz_number = 1.;
    const double vtxz_size = 10;

    //select events with only 1 vertex using information from Branch 46/47/48. Here we use Branch 48 (z-axis)
    Int_t nEvt = (Int_t)tree->GetEntries();
    cout << "Tree Entries " << nEvt << endl;


    //start loop over events
    for(Int_t i = 0; i < nEvt; ++i)
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
            if (iZeroBias == 1)
            {
                cout << "track pass zero bias" << endl;
                fVec_reco_Vtxz_size = fVec_reco_Vtxz->size();
                //vtxz_plot->Fill(fvecVtxz_size);

//---------------------Loops for gen tracks-------------------------------------------

                ntrk_gen = gen_tracks->size();
                ngen_Multi = 0;

                //looping over generated tracks
                for(int q = 0; q != ntrk_gen; ++q)
                {
                    XYZTVector gen_vec = (*gen_tracks)[q];

                    if (abs (gen_vec.Eta()) <= eta_cut && gen_vec.Pt() >= pt_cut)
                    {
                        gen_phi_histo->Fill(gen_vec.Phi());
                        gen_pt_histo->Fill(gen_vec.Pt());
                        gen_eta_histo->Fill(gen_vec.Eta());
                        ++ngen_Multi;
                    }
                }
                gen_multiplicity->Fill(ngen_Multi);


//---------------------Loops for reco tracks-------------------------------------------

                if (fVec_reco_Vtxz_size == vtxz_number)
                {
                    //cout << "number of vertex for event " << i << " is " << fVec_reco_Vtxz_size << endl;

                    //looping over vertices
                    for (int k = 0; k != fVec_reco_Vtxz_size; ++k)
                    {
                        n_reco_BeamSize = fabs((*fVec_reco_Vtxz)[k] - (*fVec_reco_VtxzBS)[0]);
                        cout << "Beam Size is " << n_reco_BeamSize << endl;

                        if (n_reco_BeamSize <= vtxz_size)
                        {
                            ntrk_reco_total = reco_tracks->size();

                            //fill with the lumi sections that meet the above event-level cuts
                            lumi_histo->Fill(nlumi_section);

                            d_reco_trk = 0;
                            ++d_reco_evt;


                            //looping over tracks
                            for (int j = 0; j != ntrk_reco_total; ++j)
                            {
                                XYZTVector reco_vec = (*reco_tracks)[j];

                                if (j == 9)
								{
									//dzleaf = (*fVec_dz)[j];
									dzcalc = ((reco_vec.Z())-(*fVec_reco_Vtxz)[k])-((((reco_vec.X())-(*fvec_reco_Vtxx)[k])*(reco_vec.px())+((reco_vec.Y())-(*fvec_reco_Vtxy)[k])*(reco_vec.py()))/(reco_vec.Pt())*(reco_vec.pz()/reco_vec.Pt()));
									//dzminus = ((vec.Z())-(*fvecVtxz)[k]);
                                    //tr_d0= (- (tr_x-vtx_x)  track.py() + (tr_y-vtx_y)  track.px() ) / track.pt()
									d0leaf = (*fVec_d0)[j];
									//d0calc = ((-(reco_vec.X() - (*fvec_reco_Vtxx)[k])*reco_vec.Py()) + ((reco_vec.Y()-(*fvec_reco_Vtxy)[k])*reco_vec.Px()))/reco_vec.Pt();
                                    d0calc = ((reco_vec.Y() - (*fvec_reco_Vtxy)[k])*reco_vec.Px() - (reco_vec.X() - (*fvec_reco_Vtxx)[k])*reco_vec.Py())/reco_vec.Pt();
								}

                                if ((*nVec_HighPurity)[j] == 1)
                                {
                                    if (abs (reco_vec.Eta()) <= eta_cut && reco_vec.Pt() >= pt_cut)
                                    {
                                        dz_dzErr = ((*fVec_dz)[j])/((*fVec_dzErr)[j]);
                                        d0_d0Err = ((*fVec_d0)[j])/((*fVec_d0Err)[j]);
                                        ptErr_pt = ((*fVec_ptErr)[j]/(reco_vec.Pt()));

                                        //fphi_eff = reco_vec.Phi()/gen_vec.Phi();
                                        //feta_eff = reco_vec.Eta()/gen_vec.Eta();
                                        //fpt_eff = reco_vec.Pt()/gen_vec.Pt();

                                        pt_efficiency->Fill(fpt_eff);
                                        eta_efficiency->Fill(feta_eff);
                                        phi_efficiency->Fill(fphi_eff);
                                        reco_phi_histo->Fill(reco_vec.Phi());
                                        reco_pt_histo->Fill(reco_vec.Pt());
                                        reco_eta_histo->Fill(reco_vec.Eta());
                                        dz_sigmadz->Fill(dz_dzErr);
                                        d0_sigmad0->Fill(d0_d0Err);
                                        sigmapt_pt->Fill(ptErr_pt);
                                        ++d_reco_trk;
                                    }
                                    ++nHigh_Purity;
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

                            reco_multiplicity->Fill(d_reco_trk);
                            ntrk_normalized += d_reco_trk;
                        }


                    }

                }

            }

        }

    }

    //cout << "Total number of selected tracks is " << ntrk_normalized << endl;
    cout << "d0leaf is " << d0leaf << " d0calc is "<< d0calc << endl;
    cout << "Total number of selected events is " << d_reco_evt << endl;
    cout << "Total number of selected tracks is " << d_reco_trk << endl;

    //vtxz_plot->Draw();
    //reco_canvas->Update();
    //reco_canvas->SaveAs("vtxz_number.png");
    //delete reco_canvas;
    //delete vtxz_plot;


    //after the loop over all events, draw the resulting plots

    reco_canvas->Divide(3,2);

    reco_canvas->cd(1);
    gPad->SetLogy();
    gPad->SetLogx();
    reco_pt_histo->Scale(1/n_reco_evt);
    reco_pt_histo->Draw();
    gen_pt_histo->SetLineColor(kRed-2);
    gen_pt_histo->Scale(1/n_reco_evt, "SAME");
    gen_pt_histo->Draw();

    reco_canvas->cd(2);
    gPad->SetLogy();
    gPad->SetLogx();
    reco_eta_histo->Scale(1/n_reco_evt);
    reco_eta_histo->Draw();
    gen_eta_histo->SetLineColor(kBlue-4);
    gen_eta_histo->Scale(1/n_reco_evt, "SAME");
    gen_eta_histo->Draw();

    reco_canvas->cd(3);
    gPad->SetLogy();
    reco_phi_histo->Scale(1/n_reco_evt);
    reco_phi_histo->Draw();
    gen_phi_histo->SetLineColor(kGreen);
    gen_phi_histo->Scale(1/n_reco_evt, "SAME");
    gen_phi_histo->Draw();

    /*reco_canvas->cd(4);
    gPad->SetLogy();
    dz_sigmadz->DrawNormalized("", 1);

    reco_canvas->cd(5);
    gPad->SetLogy();
    d0_sigmad0->DrawNormalized("", 1);

    reco_canvas->cd(6);
    gPad->SetLogy();
    sigmapt_pt->DrawNormalized("", 1);*/
    reco_canvas->Update();
    reco_canvas->SaveAs("treesCUETP8M1_66_relative_errors_compare.pdf");

    /*gen_canvas->Divide(2,2);

    gen_canvas->cd(1);
    gPad->SetLogy();
    gen_multiplicity->SetLineColor(kBlue+2);
    gen_multiplicity->DrawNormalized("", 1);
    reco_multiplicity->SetLineColor(kBlack);
    reco_multiplicity->DrawNormalized("SAME", 1);

    gen_canvas->cd(2);
    pt_efficiency->Draw();

    gen_canvas->cd(3);
    eta_efficiency->Draw();

    gen_canvas->cd(4);
    phi_efficiency->Draw();

    gen_canvas->Update();
    gen_canvas->SaveAs("multiplicity_efficiency_plots_compare.pdf");*/


}
