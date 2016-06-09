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
    TFile *myFile = TFile::Open("tree1.root", "READ");
    TTree* tree = (TTree*)myFile->Get("UETree/data");

    vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *tracks = 0;
    tree->SetBranchAddress("recoTracksp4", &tracks);
    //vector<TLorentzVector> *tracks = new vector<TLorentzVector>;

    Int_t nentries = (Int_t)tree->GetEntries();
    cout << "Tree Entries " << nentries << endl;

    TCanvas *canvas = new TCanvas;
    TH1F *event_histo = new TH1F ("reco_evt", "reco_evt", 100, 0, 200);
    TH1F *pt_histo = new TH1F ("reco_pt", "reco_pT", 200, 0, 20);
    TH1F *eta_histo = new TH1F ("reco_eta", "reco_Eta", 100, -3, 3);
    TH1F *phi_histo = new TH1F ("reco_phi", "reco_Phi", 100, -4, 4);

     //start loop over events
     for(Int_t i = 0; i < nentries; ++i)
     {
        if (i%10000==0)
            cout<< "Processing event number: " << i << endl; //1 TTree entry = 1 event
        tree->GetEntry(i);

        Int_t *iZeroBias; //data from Branch 60
        tree->GetBranchAddress("trgZeroBias",&iZeroBias);

        //we select only events that are triggered by ZeroBias Trigger. "True" in ZeroBias evaluate to 1
        if (*iZeroBias == 1)
        {
            //select events with only 1 vertex using information from Branch 46/47/48. Here we use Branch 48
            vector<float> *fvecVtxz = 0;
            tree->SetBranchAddress("vtxx",&fvecVtxz);

            if (fvecVtxz->size() == 1)
            {
                //select vertices with vertex - beamspot < 15cm
                vector<float> *fvecVtxzBS;
                tree->GetBranchAddress("vtxzBS", &fvecVtxzBS);
                vector<float> fSmallBeam = abs(*fvecVtxz - *fvecVtxzBS);

                if (fSmallBeam < 0.15)
                {
                    int ntrk = tracks->size();
                    event_histo->Fill(tracks->size());

                    if (i%10000==0)
                        cout<< tracks->size() << endl;

                    for (Int_t j = 0; j != ntrk; ++j)
                    {
                        XYZTVector vec = (*tracks)[j];

                        if (abs (vec.Eta()) < 2.4)
                        {
                            if (vec.Pt() > 0.5)
                            {
                                pt_histo->Fill(vec.Pt());
                                eta_histo->Fill(vec.Eta());
                                phi_histo->Fill(vec.Phi());
                            }

                            else
                                continue;
                        }

                        else
                            continue;
                    }

                }

                else
                    continue;
            }

            else
                continue;
        }

        else
            continue;

        cout << "all is well" << endl;

        /*


    }
    //after the loop over all events, draw the resulting plots

    canvas->Divide(2,2);

    canvas->cd(1);
    gPad->SetLogy();
    event_histo->Draw();
    canvas->Update();

    canvas->cd(2);
    gPad->SetLogy();
    pt_histo->Draw();
    canvas->Update();

    canvas->cd(3);
    eta_histo->Draw();
    canvas->Update();

    canvas->cd(4);
    gPad->SetLogy();
    phi_histo->Draw();
    canvas->Update();

    canvas->SaveAs("track_selection.pdf");
    canvas->SaveAs("track_selection.png"); */
}
