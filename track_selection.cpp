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

void track_selection()
{
   //gROOT->ProcessLine(".L Loader.C+");

    // attach "UETree/data" of tree1.root as the main root file for this program
    TFile *myFile = TFile::Open("root://eoscms.cern.ch//eos/cms/store/user/wei/multiplicity/trees_10.root", "READ");
    TTree* tree = (TTree*)myFile->Get("UETree/data");

    //for this dataset we want lumisection of 90 and above
    //data from branch 62
    int nlumi_section;
    tree->SetBranchAddress("lumi", &nlumi_section);

    int iZeroBias; //data from Branch 60
    tree->SetBranchAddress("trgZeroBias",&iZeroBias);

    //variables to check if there is only 1 vertex
    vector<float> *fvecVtxz = 0; //require initialisation to 0 to avoid crash
    tree->SetBranchAddress("vtxz",&fvecVtxz);
    float fvecVtxz_size;

    //variables to check for within 15cm
    vector<float> *fvecVtxzBS = 0; //from branch 54
    tree->SetBranchAddress("vtxzBS", &fvecVtxzBS);
    float nBeamSize;

    vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *tracks = 0;
    tree->SetBranchAddress("recoTracksp4", &tracks);
    //vector<TLorentzVector> *tracks = new vector<TLorentzVector>;

    //declare variable to hold track number
    int ntrk;

    const int lumi_cut = 90;
    const double eta_cut = 2.4;
    const double pt_cut = 0.5;
    const double vtxz_number = 1.;
    const double vtxz_size = 10;


    TCanvas *canvas = new TCanvas;
    //TH1F *event_histo = new TH1F ("reco_evt", "reco_evt", 100, 0, 200);
    TH1F *pt_histo = new TH1F ("reco_pt", "reco_pT", 200, 0, 50);
    TH1F *eta_histo = new TH1F ("reco_eta", "reco_Eta", 100, -4, 4);
    TH1F *phi_histo = new TH1F ("reco_phi", "reco_Phi", 100, -4, 4);
    TH1F *lumi_histo = new TH1F ("lumi_section", "lumi_section", 160, 80, 230);
    TH1F *multiplicity = new TH1F ("multiplicity", "Multiplicity", 250, 0, 250);

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
                fvecVtxz_size = fvecVtxz->size();
                //vtxz_plot->Fill(fvecVtxz_size);

                if (fvecVtxz_size == vtxz_number)
                {
                    cout << "number of vertex for event " << i << " is " << fvecVtxz_size << endl;

                    //looping over vertices
                    for (int k = 0; k != fvecVtxz_size; ++k)
                    {
                        nBeamSize = fabs((*fvecVtxz)[k] - (*fvecVtxzBS)[0]);
                        cout << "Beam Size is " << nBeamSize << endl;

                        if (nBeamSize <= vtxz_size)
                        {
                            ntrk = tracks->size();
                            lumi_histo->Fill(nlumi_section);

                            //looping over tracks
                            for (int j = 0; j != ntrk; ++j)
                            {
                                XYZTVector vec = (*tracks)[j];

                                (*multiplicity) = (*eta_histo)/(*phi_histo);
                                if (abs (vec.Eta()) <= eta_cut && vec.Pt() >= pt_cut)
                                {
                                    phi_histo->Fill(vec.Phi());
                                }

                                if (abs (vec.Eta()) <= eta_cut)
                                {
                                    pt_histo->Fill(vec.Pt());
                                }

                                if (vec.Pt() >= pt_cut)
                                {
                                    eta_histo->Fill(vec.Eta());
                                }


                            }

                        }

                    }

                }

            }

        }

    }




//vtxz_plot->Draw();
//canvas->Update();
//canvas->SaveAs("vtxz_number.png");
//delete canvas;
//delete vtxz_plot;


    //after the loop over all events, draw the resulting plots

    canvas->Divide(2,2);

    canvas->cd(1);
    gPad->SetLogy();
    pt_histo->Draw();
    canvas->Update();

    canvas->cd(2);
    eta_histo->Draw();
    canvas->Update();

    canvas->cd(3);
    gPad->SetLogy();
    phi_histo->Draw();
    canvas->Update();

    canvas->cd(4);
    multiplicity->Draw();
    canvas->Update();

    /*canvas->cd();
    gPad->SetLogy();
    event_histo->Draw();
    canvas->Update();*/

    //canvas->SaveAs("track_selection.pdf");
    canvas->SaveAs("track_selection_eta X phi trial.png");
}
