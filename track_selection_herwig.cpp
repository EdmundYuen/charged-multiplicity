#include "TROOT.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <algorithm>
#include "math.h"
#include <vector>
#include <iostream>
#include <fstream>
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
#include "TString.h"
#include "TChain.h"
#include "TChainElement.h"

using namespace std;
using namespace ROOT::Math;
vector<TString> *getListOfFiles(TString);
void track_selection_herwig()
{
    TH1::SetDefaultSumw2(1);

    //=========================== Implement Cuts ================================
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

    //TFile *recofile = TFile::Open("treesCUETP8M1_66.root", "READ");
    //TTree *herwigtree = (TTree*)recofile->Get("UETree/data");

//=========================== Define variables to read TTree ================

    vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *reco_tracks = 0;
    int nreco_lumi;
    int nreco_zerobias;
    vector<float> *fvecreco_vtxx = 0;
    vector<float> *fvecreco_vtxy = 0;
    vector<float> *fvecreco_vtxz = 0;
    vector<float> *fvecreco_vtxxBS = 0;
    vector<float> *fvecreco_vtxyBS = 0;
    vector<float> *fvecreco_vtxzBS = 0;
    vector<int> *nvecreco_highpurity = 0;
    vector<int> *nvecreco_vtxndof = 0;
    vector<float> *fvecreco_dz = 0;
    vector<float> *fvecreco_d0 = 0;
    vector<float> *fvecreco_dzerr = 0;
    vector<float> *fvecreco_d0err = 0;
    vector<float> *fvecreco_pterr = 0;
    vector<float> *fvecreco_vtxxerr = 0;
    vector<float> *fvecreco_vtxyerr = 0;
    vector<float> *fvecreco_vtxzerr = 0;
    vector<int> *nvecreco_validhits = 0;
    vector<float> *fvecreco_trackschi2n = 0;
    vector<int> *fvecreco_isfake = 0;
    vector<float> *fvecreco_trkx = 0;
    vector<float> *fvecreco_trky = 0;
    vector<float> *fvecreco_trkz = 0;

//============================================== Histos for pT, eta, phi ==============================================================

    TH1F *reco_pt_histo = new TH1F ("reco pT", "Working Plot #sqrt{s} = 13TeV", 200, 0, 10);
    TH1F *reco_eta_histo = new TH1F("reco eta", "Working Plot #sqrt{s} = 13TeV", 50, -2.5, 2.5);
    TH1F *reco_phi_histo = new TH1F ("reco phi", "Working Plot #sqrt{s} = 13TeV", 60, -3, 3);

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

    TH1F *reco_multiplicity = new TH1F("Normalized_Multiplicity", "Working Plot #sqrt{s} = 13TeV", 200, 0, 200);
    TH1F *reco_vtxzminusvtxz = new TH1F ("vtxzminusvtxz", "vtxzminusvtxz", 800, -4, 4);
    TH1F *reco_vtxzposn = new TH1F ("vtxzposn", "vtxzpos^{n}", 400, -20, 20);

    //=========================== Retrieve ROOT file============================

    vector<TString> *herwigfiles = new vector<TString>();
    cout << "Getting list of files..." << endl;
    //herwigfiles = getListOfFiles("filelistherwig.txt");
    herwigfiles = getListOfFiles("filelistherwigtest.txt");
    cout << "File list stored" << endl;

    TFile *herwigfile;
    TTree *herwigtree;

    //========================================================= Start of Herwig File Loop ================================================================

    //while (getline(filedata, datastr))
    for(vector<TString>::iterator itlistherwigfiles = herwigfiles->begin() ; itlistherwigfiles != herwigfiles->end(); ++itlistherwigfiles)
    {
        cout << "Opening new file " << *itlistherwigfiles << endl;
        //TString Tdatastr(datastr);
        //herwigfile = TFile::Open("treesCUETP8M1_66.root", "READ");
        //herwigtree = (TTree*)herwigfile->Get("UETree/data");
        //datafile = new TFile(*itlistdatafiles, "READ");
        //datafile = TFile::Open("root://eoscms.cern.ch//eos/cms/store/user/wei/multiplicity/data/ZeroBias1_trees_10.root", "READ");
        //datafile = new TFile(*itlistdatafiles, "READ"); //fail "no matching function for call to 'Open'"
        herwigfile = TFile::Open(*itlistherwigfiles, "READ");
        cout << "Opened " << *itlistherwigfiles << endl;
        herwigtree = (TTree*)herwigfile->Get("UETree/data");
        cout<< "Congratulations you have succeeded in looping over the damn Herwig files!\n";

        herwigtree->SetBranchAddress("recoTracksp4", &reco_tracks);
        herwigtree->SetBranchAddress("lumi", &nreco_lumi);
        herwigtree->SetBranchAddress("trgZeroBias", &nreco_zerobias);
        herwigtree->SetBranchAddress("vtxx", &fvecreco_vtxx);
        herwigtree->SetBranchAddress("vtxy", &fvecreco_vtxy);
        herwigtree->SetBranchAddress("vtxz", &fvecreco_vtxz);
        herwigtree->SetBranchAddress("vtxxBS", &fvecreco_vtxxBS);
        herwigtree->SetBranchAddress("vtxyBS", &fvecreco_vtxyBS);
        herwigtree->SetBranchAddress("vtxzBS", &fvecreco_vtxzBS);
        herwigtree->SetBranchAddress("recoTrackshighPurity", &nvecreco_highpurity);
        herwigtree->SetBranchAddress("vtxndof", &nvecreco_vtxndof);
        herwigtree->SetBranchAddress("recoTracksdz", &fvecreco_dz);
        herwigtree->SetBranchAddress("recoTracksd0", &fvecreco_d0);
        herwigtree->SetBranchAddress("recoTracksdzErr", &fvecreco_dzerr);
        herwigtree->SetBranchAddress("recoTracksd0Err", &fvecreco_d0err);
        herwigtree->SetBranchAddress("recoTracksptErr", &fvecreco_pterr);
        herwigtree->SetBranchAddress("vtxxErr", &fvecreco_vtxxerr);
        herwigtree->SetBranchAddress("vtxyErr", &fvecreco_vtxyerr);
        herwigtree->SetBranchAddress("vtxzErr", &fvecreco_vtxzerr);
        herwigtree->SetBranchAddress("recoTracksnValidHits", &nvecreco_validhits);
        herwigtree->SetBranchAddress("recoTrackschi2n", &fvecreco_trackschi2n);
        herwigtree->SetBranchAddress("vtxisFake", &fvecreco_isfake);
        herwigtree->SetBranchAddress("recoTracksvx", &fvecreco_trkx);
        herwigtree->SetBranchAddress("recoTracksvy", &fvecreco_trky);
        herwigtree->SetBranchAddress("recoTracksvz", &fvecreco_trkz);

        Int_t nreco_totalEvt = (Int_t)herwigtree->GetEntries();
        cout << "There is a total of " << nreco_totalEvt << " events." << endl;

        for (int ev = 0; ev < nreco_totalEvt; ++ev)
        {
            herwigtree->GetEntry(ev);

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

    //========================================================= Start of Herwig Evt Loop ================================================================

        for (Int_t i = 0; i < nreco_totalEvt; ++i)
        {
            herwigtree->GetEntry(i);
            //cout << "At entry " << i << endl;

            if (nreco_lumi >= lumi_cut)
            {
                if (nreco_zerobias == 1)
                {
                    nreco_totaltrk = reco_tracks->size();
                    ++freco_evt;
                    freco_multiplicity = 0;

                    int vtxdof = 0;

                    nreco_numberofvtxx = fvecreco_vtxx->size();
                    nreco_numberofvtxy = fvecreco_vtxy->size();
                    nreco_numberofvtxz = fvecreco_vtxz->size();


    //========================================================= Start of Vertex Loop ================================================================

                    for (int vtxnumber = 0; vtxnumber != nreco_numberofvtxx; ++vtxnumber)
                    {
                        if((*nvecreco_vtxndof)[vtxdof] > dof_cut)
                        {
                            freco_vtxxysize = sqrt(pow(fabs(((*fvecreco_vtxx)[vtxnumber]) - ((*fvecreco_vtxxBS)[vtxnumber])), 2) + pow(fabs((*fvecreco_vtxy)[vtxnumber] - (*fvecreco_vtxyBS)[vtxnumber]), 2));
                            freco_vtxzsize = fabs(((*fvecreco_vtxz)[vtxnumber]) - ((*fvecreco_vtxzBS)[vtxnumber]));

                            //if (freco_vtxxysize <= vtxxysize && freco_vtxzsize <= vtxzsize)
                            //{
                                if (nreco_numberofvtxz == vtx_number_cut && (*fvecreco_isfake)[0] == 0)
                                {
                                    reco_vtxzposn->Fill((*fvecreco_vtxz)[vtxnumber]);
                                    ++freco_numselectedvtxz;

    //========================================================= Start of Trk Loop ================================================================

                                    for (int t = 0; t != nreco_totaltrk; ++t)
                                    {
                                        XYZTVector reco_vec = (*reco_tracks)[t];
                                        //cout << "Within track " << endl;
                                        //using formula from paper of run 1 result, dz is leaf value

                            //======================================= dz ===================================================

                                        //freco_dz = (*fvecreco_dz)[t];
                                        freco_dz = ((*fvecreco_trkz)[t] - (*fvecreco_vtxz)[t]) - (((*fvecreco_trkx)[t] - (*fvecreco_vtxx)[t])*reco_vec.Px() + ((*fvecreco_trky)[t] - (*fvecreco_vtxy)[t])*reco_vec.Py())/reco_vec.Pt()*(reco_vec.Pz()/reco_vec.Pt());//WY
                                        reco_dzleaf->Fill((*fvecreco_dz)[t]);
                                        freco_sigmadz = sqrt(pow(((*fvecreco_dzerr)[t]), 2) + pow(freco_wz, 2));

                                        freco_dz_sigmadz = (freco_dz)/((*fvecreco_dzerr)[t]); //leaf
                                        freco_dz_sigmadzcalc = ((freco_dz)/sqrt(pow(((*fvecreco_dzerr)[t]),2)+pow((*fvecreco_vtxzerr)[t],2))); //WY
                                        freco_dz_sigmadzrun1 = ((freco_dz) / (freco_sigmadz)); //run 1 formula

                            //======================================= d0 ===================================================

                                        //freco_d0 = (*fvecreco_d0)[t];
                                        freco_d0 = (-((*fvecreco_trkx)[t] - (*fvecreco_vtxx)[t])*reco_vec.Py() + ((*fvecreco_trky)[t] - (*fvecreco_vtxy)[t])*reco_vec.Px())/reco_vec.Pt();//WY
                                        freco_sigmad0run1 = sqrt(pow(((*fvecreco_d0err)[t]), 2) + (freco_wx)*(freco_wy));
                                        freco_sigmad0calc = sqrt(pow(((*fvecreco_d0err)[t]),2)+pow((*fvecreco_vtxxerr)[t],2)+pow((*fvecreco_vtxyerr)[t],2)); //WY

                                        freco_d0_sigmad0 = (freco_d0)/((*fvecreco_d0err)[t]); //both leaf values
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
                    ++freco_evt;
                }
            }
        }
    //========================================================= End of Herwig Evt Loop ================================================================

    }

    //========================================================= End of Herwig File Loop ================================================================
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
    reco_eta_histo->GetYaxis()->SetTitle("#frac{1}{N_{ch}} #frac{dN_{ch}}{d#eta}");
    reco_eta_histo->SetLineColor(4);
    reco_eta_histo->Draw();
    reco_eta_histo->Write();

    //canvas->cd(3);
    reco_phi_histo->Scale(1/freco_trkphi);
    reco_phi_histo->SetMinimum(0.014);
    reco_phi_histo->SetMaximum(0.024);
    reco_phi_histo->GetXaxis()->SetTitle("#phi");
    reco_phi_histo->GetYaxis()->SetTitleOffset(1.3);
    reco_phi_histo->GetYaxis()->SetTitle("#frac{1}{N_{ch}} #frac{dN_{ch}}{d#phi}");
    reco_phi_histo->SetLineColor(4);
    reco_phi_histo->Write();

    //canvas->cd(1);
    gPad->SetLogy();
    reco_pt_histo->Scale(1/freco_trkpt);
    reco_pt_histo->GetXaxis()->SetTitle("p_{T}");
    reco_pt_histo->GetYaxis()->SetTitleOffset(1.3);
    reco_pt_histo->GetYaxis()->SetTitle("#frac{1}{N_{ch}} #frac{dN_{ch}}{d#p_{T}}");
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

    reco_d0_sigmad0->SetMinimum(1E-5);
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
    reco_multiplicity->Scale(1/freco_evt);
    reco_multiplicity->GetXaxis()->SetTitle("N_{ch}");
    reco_multiplicity->GetYaxis()->SetTitleOffset(1.1);
    reco_multiplicity->GetYaxis()->SetTitle("#frac{1}{N_{ev}} #frac{dN_{ch}}{dN_{ev}}");
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
    recoleg_multiplicity->SetFillStyle(0);
    recoleg_multiplicity->SetBorderSize(0);
    recoleg_multiplicity->SetTextSize(0.02);
    recoleg_multiplicity->AddEntry(reco_multiplicity, "reco", "lf");
    recoleg_multiplicity->Draw();
    recoleg_multiplicity->Write();

    reco_pad->Modified();
    reco_pad->Update();
    reco_plot.Write();
}

vector<TString> *getListOfFiles(TString strfiles)
{

    vector<TString>* vfiles = new vector<TString>;

    if(strfiles.Contains(".root"))
    {
        TChain chain("evt","");
        chain.Add(strfiles);
        TObjArray* fileElements=chain.GetListOfFiles();
        TIter next(fileElements);
        TChainElement *chEl=0;
        while (( chEl=(TChainElement*)next() ))
            vfiles->push_back(TString(chEl->GetTitle()));
    }
    else if(strfiles.Contains(".txt"))
    {
        ifstream txtfile;
        txtfile.open(strfiles);
        if(!txtfile)
        {
            cout<<"Unable to read the txt file where the rootfiles are." << endl ;
            cout << strfiles << " doesn't exist." << endl << "Aborting ...";
            exit(0);
        }

        string filename;

        while(txtfile>>filename && filename!="EOF")
        vfiles->push_back(TString(filename));
        txtfile.close();
    }
    else
    {
        cout << "Unknown type of input to get files. Must contain either .root or .txt extension." << endl << "Aborting ..." << endl;
        exit(0);
    }

    cout << "[getListOfFiles] Will run on " << vfiles->size() << " files" << endl;
    return vfiles;
}
