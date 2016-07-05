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

    vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *gen_tracks = 0;
    tree->SetBranchAddress("genParticlesp4", &gen_tracks);

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

    //declare variable to hold track number
    int ntrk, nMulti;
    int ntrk_normalized = 0;
    int nHigh_Purity = 0;
    float dz_dzErr, d0_d0Err, pt_ptErr;

    const int lumi_cut = 90;
    const double eta_cut = 2.4;
    const double pt_cut = 0.5;
    const double vtxz_number = 1.;
    const double vtxz_size = 10;
    const int ntracks = 533084;


    TCanvas *canvas = new TCanvas;
    //TH1F *event_histo = new TH1F ("reco_evt", "reco_evt", 100, 0, 200);
    TH1F *pt_histo = new TH1F ("reco_pt", "Normalized_reco_pT", 100, 0, 100);
    TH1F *eta_histo = new TH1F ("reco_eta", "Normalized_reco_Eta", 25, -3, 3);
    TH1F *phi_histo = new TH1F ("reco_phi", "Normalized_reco_Phi", 100, -4, 4);
    TH1F *lumi_histo = new TH1F ("lumi_section", "lumi_section", 160, 80, 230);
    TH1F *multiplicity = new TH1F ("multiplicity", "Normalized_Multiplicity", 200, 0, 200);
    TH1F *dz_sigmadz = new TH1F ("dz_sigmadz", "Normalized_dz_sigmadz", 200, 0, 200);
    TH1F *d0_sigmad0 = new TH1F ("d0_sigmad0", "Normalized_d0_sigmad0", 200, 0, 200);
    TH1F *pt_sigmapt = new TH1F ("pt_sigmapt", "Normalized_pt_sigmapt", 200, 0, 200);
    //TH1F *normalized_multiplicity_histo = new TH1F ("normalized_multiplicity", "normalized_multiplicity", 200, 0, 200);

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
        cout << "set branch address for zero bias" << endl;

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

                if (fVec_reco_Vtxz_size == vtxz_number)
                {
                    cout << "number of vertex for event " << i << " is " << fVec_reco_Vtxz_size << endl;

                    //looping over vertices
                    for (int k = 0; k != fVec_reco_Vtxz_size; ++k)
                    {
                        n_reco_BeamSize = fabs((*fVec_reco_Vtxz)[k] - (*fVec_reco_VtxzBS)[0]);
                        cout << "Beam Size is " << n_reco_BeamSize << endl;

                        if (n_reco_BeamSize <= vtxz_size)
                        {
                            ntrk = reco_tracks->size();

                            //fill with the lumi sections that meet the above event-level cuts
                            lumi_histo->Fill(nlumi_section);

                            nMulti = 0;

                            //looping over tracks
                            for (int j = 0; j != ntrk; ++j)
                            {
                                XYZTVector vec = (*reco_tracks)[j];
                                dz_dzErr = ((*fVec_dz)[j])/((*fVec_dzErr)[j]);
                                d0_d0Err = ((*fVec_d0)[j])/((*fVec_d0Err)[j]);
                                pt_ptErr = ((vec.Pt())/(*fVec_ptErr)[j]);

                                if ((*nVec_HighPurity)[j] == 1)
                                {
                                    if (abs (vec.Eta()) <= eta_cut && vec.Pt() >= pt_cut)
                                    {
                                        phi_histo->Fill(vec.Phi());
                                        pt_histo->Fill(vec.Pt());
                                        eta_histo->Fill(vec.Eta());
                                        dz_sigmadz->Fill(dz_dzErr);
                                        d0_sigmad0->Fill(d0_d0Err);
                                        pt_sigmapt->Fill(pt_ptErr);
                                        ++nMulti;
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

                            multiplicity->Fill(nMulti);
                            ntrk_normalized += nMulti;
                        }


                    }

                }

            }

        }

    }

    cout << "Total number of selected tracks is " << ntrk_normalized << endl;

//vtxz_plot->Draw();
//canvas->Update();
//canvas->SaveAs("vtxz_number.png");
//delete canvas;
//delete vtxz_plot;


    //after the loop over all events, draw the resulting plots

    canvas->Divide(3,2);
    /*canvas->cd(1);
    gPad->SetLogy();
    dz_sigmadz->DrawNormalized("", 1);

    canvas->cd(2);
    gPad->SetLogy();
    d0_sigmad0->DrawNormalized("", 1);

    gPad->SetLogy();
    pt_sigmapt->DrawNormalized("", 1);
    canvas->Update();*/

    canvas->cd(1);
    gPad->SetLogy();
    pt_histo->DrawNormalized("", 1);

    canvas->cd(2);
    gPad->SetLogy();
    eta_histo->DrawNormalized("", 1);

    canvas->cd(3);
    gPad->SetLogy();
    phi_histo->DrawNormalized("", 1);

    canvas->cd(4);
    gPad->SetLogy();
    dz_sigmadz->DrawNormalized("", 1);

    canvas->cd(5);
    gPad->SetLogy();
    d0_sigmad0->DrawNormalized("", 1);

    canvas->cd(6);
    gPad->SetLogy();
    pt_sigmapt->DrawNormalized("", 1);
    canvas->Update();

    canvas->SaveAs("treesCUETP8M1_66_reco_relative_errors.pdf");
}
