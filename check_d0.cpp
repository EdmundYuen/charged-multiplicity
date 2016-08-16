#include "TROOT.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
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

void check_d0()
{
    //===========================Implement cuts================================
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

    //==============================================Variables==============================================================
    float fdata_vtxxysize, fdata_vtxzsize;
    float fdata_totaltrk, fdata_multiplicity;

    //===========================Retrieve ROOT file============================

    //TFile *datafile = TFile::Open("treesCUETP8M1_66.root", "READ");
    TFile *datafile = TFile::Open("tree1.root", "READ");
    TTree *datatree = (TTree*)datafile->Get("UETree/data");

    //===========================Define variables to read TTree================

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

    //==============================================Histograms==============================================================

    TH1F *vtxx_vs_evts = new TH1F ("vtxx_vs_evts", "vtxx_vs_evts", 100, -0.5, 0.5);
    TH2F *vtxx_vs_evts2D = new TH2F ("vtxx_vs_evts2D", "vtxx_vs_evts2D", 100, -0.5, 0.5, 100000, 0, 100000);
    TH1F *vtxy_vs_evts = new TH1F ("vtxy_vs_evts", "vtxy_vs_evts", 100, -0.5, 0.5);
    TH2F *vtxy_vs_evts2D = new TH2F ("vtxy_vs_evts2D", "vtxy_vs_evts2D", 100, -0.5, 0.5, 100000, 0, 100000);
    TH1F *vtxz_vs_evts = new TH1F ("vtxz_vs_evts", "vtxz_vs_evts", 400, -10, 10);
    TH2F *vtxz_vs_evts2D = new TH2F ("vtxz_vs_evts2D", "vtxz_vs_evts2D", 400, -2, -1, 100000, 0, 100000);
    TH2F *vtxz_vs_evts2Dnocut = new TH2F ("vtxz_vs_evts2Dnocut", "vtxz_vs_evts2Dnocut", 400, -2, -1, 100000, 0, 100000);
    TH2F *vtxz_vs_evts2D_aftvtxselection = new TH2F ("vtxz_vs_evts2D_with_vtxselection", "vtxz_vs_evts2D with vtx selection", 400, -2, -1, 100000, 0, 100000);
    TH3F *vtxz_vs_evts_vs_multiplicity = new TH3F ("vtxz_vs_evts3D", "vtxz_vs_evts3D", 900, -2, 1, 100000, 0, 100000, 300, 0, 300);

    TH1F *data_d0 = new TH1F("data_d0", "data d_{0}", 400, -1, 1);
    TH1F *data_d0nocut = new TH1F("data_d0nocut", "data d_{0} nocut", 200, -1, 1);
    TH1F *data_d0_sigmad0 = new TH1F ("d0_sigmad0", " #frac{d_{0}}{#sigma_{xy}}", 400, -20, 20);
    TH1F *data_dz = new TH1F("data_dz", "data d_{z}", 200, -10, 10);
    TH1F *data_dznocut = new TH1F("data_dznocut", "data d_{z} nocut", 400, -10, 10);
    TH1F *data_d0err = new TH1F ("data d0err", "data #sigma_{0}", 120, 0, 1.2);


    int nevt = (Int_t)datatree->GetEntries();
    cout << "There is a total of " << nevt << " entries." << endl;

    for (int evt = 0; evt != nevt; ++evt)
    {
        datatree->GetEntry(evt);

        if (ndata_lumi >= lumi_cut)
        {
            if (ndata_zerobias == 1)
            {
                for (int vtxx = 0; vtxx < fvecdata_vtxx->size(); ++vtxx)
                {
                    vtxx_vs_evts->Fill((*fvecdata_vtxx)[vtxx]);
                    vtxx_vs_evts2D->Fill((*fvecdata_vtxx)[vtxx], evt);
                }

                for (int vtxy = 0; vtxy < fvecdata_vtxy->size(); ++vtxy)
                {
                    vtxy_vs_evts->Fill((*fvecdata_vtxy)[vtxy]);
                    vtxy_vs_evts2D->Fill((*fvecdata_vtxy)[vtxy], evt);
                }

                fdata_totaltrk = data_tracks->size();
                //apply cuts to vtxz
                //loop through vertex
                for (int vtxz = 0; vtxz < fvecdata_vtxz->size(); ++vtxz)
                {
                    vtxz_vs_evts2Dnocut->Fill((*fvecdata_vtxz)[vtxz], evt);
                    vtxz_vs_evts->Fill((*fvecdata_vtxz)[vtxz]);

                    fdata_vtxxysize = sqrt(pow(((*fvecdata_vtxx)[vtxz]) - ((*fvecdata_vtxxBS)[vtxz]), 2) + pow(((*fvecdata_vtxy)[vtxz]) - ((*fvecdata_vtxyBS)[vtxz]), 2));
                    fdata_vtxzsize = fabs(((*fvecdata_vtxz)[vtxz]) - ((*fvecdata_vtxzBS)[vtxz]));
                    //fdata_multiplicity = 0;
                    //data_vtxz-BS->Fill((*fvecdata_vtxz)[vtxz]-(*fvecdata_vtxzBS)[vtxz]);

                    //if ((*fvecdata_vtxz)[vtxz] < -1.69 || (*fvecdata_vtxz)[vtxz] > -1.68)
                    //if ((*fvecdata_vtxz)[vtxz] < -1.69)
                    //if ((*fvecdata_vtxz)[vtxz] < -1.506 && (*fvecdata_vtxz)[vtxz] > -1.51)
                    //{
                        //if ((*fvecdata_vtxz)[vtxz] < -1.63 || (*fvecdata_vtxz)[vtxz] > -1.62)
                        //{
                            //if ((*fvecdata_vtxz)[vtxz] >= -1.505)
                            vtxz_vs_evts2D->Fill(evt, (*fvecdata_vtxz)[vtxz]);

                            //looping through tracks
                            for (int trk = 0; trk < data_tracks->size(); ++trk)
                            {
                                XYZTVector data_trk = (*data_tracks)[trk];
                                data_d0nocut->Fill((*fvecdata_d0)[trk]);
                                data_dznocut->Fill((*fvecdata_dz)[trk]);
                                data_d0err->Fill((*fvecdata_d0err)[trk]);
                                //++fdata_multiplicity;

                            }

                        //}
                        //cout << "Just before filling 3D-histo." << endl;
                        //vtxz_vs_evts_vs_multiplicity->Fill(evt, ((*fvecdata_vtxz)[vtxz]), fdata_multiplicity);
                        //cout << "Just after filling 3D-histo." << endl;
                    //}



                }
            }
        }
    }
    TFile vtxposition ("Histos/vtxposn.root", "recreate");

    vtxx_vs_evts->Draw();
    vtxx_vs_evts->Write();
    vtxx_vs_evts2D->Draw();
    vtxx_vs_evts2D->Write();

    vtxy_vs_evts->Draw();
    vtxy_vs_evts->Write();
    vtxy_vs_evts2D->Draw();
    vtxy_vs_evts2D->Write();

    vtxz_vs_evts_vs_multiplicity->Draw();
    vtxz_vs_evts_vs_multiplicity->Write();

    vtxz_vs_evts->Draw();
    vtxz_vs_evts->Write();
    vtxz_vs_evts2D->Draw();
    vtxz_vs_evts2D->Write();
    vtxz_vs_evts2Dnocut->Draw();
    vtxz_vs_evts2Dnocut->Write();

    data_d0nocut->Draw();
    data_d0nocut->Write();
    data_d0err->SetLineColor(kBlack);
    data_d0err->Scale(1/fdata_totaltrk);
    data_d0err->Write();


    data_dznocut->Draw();
    data_dznocut->Write();
}
