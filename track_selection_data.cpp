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


void track_selection_redo()
{
    TH1::SetDefaultSumw2(1);

//=========================== Filling vector<TString> to loop over files ====================================================

    //vector<TString> *vfiles = new vector<TString>();


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
    float fdata_multiplicity_norm = 0;
    float freco_multiplicity_norm = 0;
    float freco_totalevt, fdata_totalevt;

    //TH1F *data_multiplicity = new TH1F("Normalized_Data_Multiplicity", "Normalized Data Multiplicity", 200, 0, 200);
    //TH1F *reco_multiplicity = new TH1F("Normalized_Reco_Multiplicity", "Normalized Reco Multiplicity", 200, 0, 200);

    //int nselect = 0;
    //cout << "Enter 0 for data, 1 for Herwig or 2 for both: ";
    //cin >> nselect;

    //while (nselect > 2)
    //{
        //cout << "Please enter 0, 1 or 2 only: ";
        //cin >> nselect;
    //}
//==================================================== Data Loop ===================================================================

    //if ((nselect == 0) || (nselect == 2))
    //{
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
        int ndata_lumi;
        int ndata_zerobias;

    //=========================== Declaration of variables to read TTree ================

        vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *data_tracks = 0; //pointer declaration has to be assigned to 0 to qualify as mull pointer
        vector<float> *fvecdata_vtxx = 0;
        vector<float> *fvecdata_vtxy = 0;
        vector<float> *fvecdata_vtxz = 0;
        vector<float> *fvecdata_vtxxBS = 0;
        vector<float> *fvecdata_vtxyBS = 0;
        vector<float> *fvecdata_vtxzBS = 0;
        vector<int> *nvecdata_highpurity = 0;
        vector<int> *nvecdata_vtxndof = 0;
        vector<float> *fvecdata_dz = 0;
        vector<float> *fvecdata_d0 = 0;
        vector<float> *fvecdata_dzerr = 0;
        vector<float> *fvecdata_d0err = 0;
        vector<float> *fvecdata_pterr = 0;
        vector<float> *fvecdata_vtxxerr = 0;
        vector<float> *fvecdata_vtxyerr = 0;
        vector<float> *fvecdata_vtxzerr = 0;
        vector<int> *nvecdata_validhits = 0;
        vector<int> *nvecdata_isfake = 0;
        vector<float> *fvecdata_trackschi2n = 0;
        vector<float> *fvecdata_trkx = 0;
        vector<float> *fvecdata_trky = 0;
        vector<float> *fvecdata_trkz = 0;

    //============================================== Histos for pT, eta, phi ==============================================================

        TH1F *data_pt_histo = new TH1F ("data pT", "Working Plot #sqrt{s} = 13TeV", 120, 0, 6);
        TH1F *data_eta_histo = new TH1F("data eta", "Working Plot #sqrt{s} = 13TeV", 50, -2.5, 2.5);
        TH1F *data_phi_histo = new TH1F ("data phi", "Working Plot #sqrt{s} = 13TeV", 60, -3, 3);

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

        TH1F *data_multiplicity = new TH1F("Normalized_Multiplicity", "Working Plot #sqrt{s} = 13TeV", 200, 0, 200);
        TH1F *data_vtxzminusvtxz = new TH1F ("vtxzminusvtxz", "vtxzminusvtxz", 800, -4, 4);
        TH1F *data_vtxzposn = new TH1F ("vtxzposn", "vtxzpos^{n}", 400, -20, 20);

    //=========================== Retrieve ROOT file============================

        vector<TString> *datafiles = new vector<TString>();
        cout << "Getting list of files..." << endl;
        //datafiles = getListOfFiles("filelistdata.txt");
        datafiles = getListOfFiles("filelistdatatest.txt");
        cout << "File list stored" << endl;

        TFile *datafile;
        TTree *datatree;
    //while (getline(filedata, datastr))
    for(vector<TString>::iterator itlistdatafiles = datafiles->begin() ; itlistdatafiles != datafiles->end(); ++itlistdatafiles)
    {
        cout << "Opening new file " << *itlistdatafiles << endl;
        //TString Tdatastr(datastr);
        //datafile = TFile::Open("tree1.root", "READ");
        //datatree = (TTree*)datafile->Get("UETree/data");
        datafile = new TFile(*itlistdatafiles, "READ");
        //datafile = TFile::Open("root://eoscms.cern.ch//eos/cms/store/user/wei/multiplicity/data/ZeroBias1_trees_10.root", "READ");
        //datafile = new TFile(*itlistdatafiles, "READ"); //fail "no matching function for call to 'Open'"
        //datafile = TFile::Open(*itlistdatafiles, "READ");
        cout << "Opened " << *itlistdatafiles << endl;
        datatree = (TTree*)datafile->Get("UETree/data");
        cout << "Congratulations you have succeeded in looping over the damn data files!\n";
    //============================================== Assignment of TTree Branches ====================================================

        datatree->SetBranchAddress("recoTracksp4", &data_tracks);
        datatree->SetBranchAddress("lumi", &ndata_lumi);
        datatree->SetBranchAddress("trgZeroBias", &ndata_zerobias);
        datatree->SetBranchAddress("vtxx", &fvecdata_vtxx);
        datatree->SetBranchAddress("vtxy", &fvecdata_vtxy);
        datatree->SetBranchAddress("vtxz", &fvecdata_vtxz);
        datatree->SetBranchAddress("vtxxBS", &fvecdata_vtxxBS);
        datatree->SetBranchAddress("vtxyBS", &fvecdata_vtxyBS);
        datatree->SetBranchAddress("vtxzBS", &fvecdata_vtxzBS);
        datatree->SetBranchAddress("recoTrackshighPurity", &nvecdata_highpurity);
        datatree->SetBranchAddress("vtxndof", &nvecdata_vtxndof);
        datatree->SetBranchAddress("recoTracksdz", &fvecdata_dz);
        datatree->SetBranchAddress("recoTracksd0", &fvecdata_d0);
        datatree->SetBranchAddress("recoTracksdzErr", &fvecdata_dzerr);
        datatree->SetBranchAddress("recoTracksd0Err", &fvecdata_d0err);
        datatree->SetBranchAddress("recoTracksptErr", &fvecdata_pterr);
        datatree->SetBranchAddress("vtxxErr", &fvecdata_vtxxerr);
        datatree->SetBranchAddress("vtxyErr", &fvecdata_vtxyerr);
        datatree->SetBranchAddress("vtxzErr", &fvecdata_vtxzerr);
        datatree->SetBranchAddress("recoTracksnValidHits", &nvecdata_validhits);
        datatree->SetBranchAddress("recoTrackschi2n", &fvecdata_trackschi2n);
        datatree->SetBranchAddress("vtxisFake", &nvecdata_isfake);
        datatree->SetBranchAddress("recoTracksvx", &fvecdata_trkx);
        datatree->SetBranchAddress("recoTracksvy", &fvecdata_trky);
        datatree->SetBranchAddress("recoTracksvz", &fvecdata_trkz);

    //============================================== End of Assignment TTree Branches ====================================================

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

                    ndata_numberofvtxx = fvecdata_vtxx->size();
                    ndata_numberofvtxy = fvecdata_vtxy->size();
                    ndata_numberofvtxz = fvecdata_vtxz->size();
                    fdata_multiplicity = 0;
                    ++fdata_evt;

    //========================================================= Start of Vertex Loop ================================================================

                    for (int vtxnumber = 0; vtxnumber != ndata_numberofvtxx; ++vtxnumber)
                    {
                        if((*nvecdata_vtxndof)[vtxdof] > dof_cut)
                        {
                            fdata_vtxxysize = sqrt(pow(fabs(((*fvecdata_vtxx)[vtxnumber]) - ((*fvecdata_vtxxBS)[vtxnumber])), 2) + pow(fabs((*fvecdata_vtxy)[vtxnumber] - (*fvecdata_vtxyBS)[vtxnumber]), 2));
                            fdata_vtxzsize = fabs(((*fvecdata_vtxz)[vtxnumber]) - ((*fvecdata_vtxzBS)[vtxnumber]));

                            //if (fdata_vtxxysize <= vtxxysize && fdata_vtxzsize <= vtxzsize)
                            //{
                            if (ndata_numberofvtxz == vtx_number_cut && (*nvecdata_isfake)[0] == 0)
                            {
                                data_vtxzposn->Fill((*fvecdata_vtxz)[vtxnumber]);
                                ++fdata_numselectedvtxz;

    //========================================================= Start of Trk Loop ================================================================

                                for (int t = 0; t != ndata_totaltrk; ++t)
                                {
                                    XYZTVector data_vec = (*data_tracks)[t];
                                    //cout << "Within track " << endl;
                                    //using formula from paper of run 1 result, dz is leaf value

                        //======================================= dz ===================================================

                                    //fdata_dz = (*fvecdata_dz)[t];
                                    fdata_dz = ((*fvecdata_trkz)[t] - (*fvecdata_vtxz)[t]) - (((*fvecdata_trkx)[t] - (*fvecdata_vtxx)[t])*data_vec.Px() + ((*fvecdata_trky)[t] - (*fvecdata_vtxy)[t])*data_vec.Py())/data_vec.Pt()*(data_vec.Pz()/data_vec.Pt());//WY
                                    data_dzleaf->Fill((*fvecdata_dz)[t]);

                                    fdata_sigmadz = sqrt(pow(((*fvecdata_dzerr)[t]), 2) + pow(fdata_wz, 2));
                                    //fdata_dz_sigmadz = ((*fvecdata_dz)[t])/((*fvecdata_dzerr)[t]); //leaf
                                    fdata_dz_sigmadz = fdata_dz /((*fvecdata_dzerr)[t]);
                                    //fdata_dz_sigmadzcalc = ((fdata_dz)/sqrt(pow(((*fvecdata_dzerr)[t]),2)+pow((*fvecdata_vtxzerr)[t],2))); //WY
                                    fdata_dz_sigmadzrun1 = ((fdata_dz) / (fdata_sigmadz)); //run 1 formula

                        //======================================= d0 ===================================================

                                    //fdata_d0 = (*fvecdata_d0)[t];
                                    fdata_d0 = (-((*fvecdata_trkx)[t] - (*fvecdata_vtxx)[t])*data_vec.Py() + ((*fvecdata_trky)[t] - (*fvecdata_vtxy)[t])*data_vec.Px())/data_vec.Pt(); //WY
                                    fdata_sigmad0run1 = sqrt(pow(((*fvecdata_d0err)[t]), 2) + (fdata_wx)*(fdata_wy));
                                    fdata_sigmad0calc = sqrt(pow(((*fvecdata_d0err)[t]),2)+pow((*fvecdata_vtxxerr)[t],2)+pow((*fvecdata_vtxyerr)[t],2)); //WY

                                    //fdata_d0_sigmad0 = ((*fvecdata_d0)[t])/((*fvecdata_d0err)[t]); //both leaf values
                                    fdata_d0_sigmad0 = fdata_d0/((*fvecdata_d0err)[t]);
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
                                        if (data_vec.Pt() >= pt_cut && fabs(fdata_sigmapt_pt) < ptErr_pt_cut)
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

    }

    //========================================================= End of File Loop ================================================================
        cout << "Before plotting." << endl;
        cout << "wx is " << fdata_wx << endl;
        cout << "wy is " << fdata_wy << endl;
        cout << "wz is " << fdata_wz << endl;
        cout << "cut on d0 is " << sqrt(0.0025 - fdata_wx*fdata_wy) << endl;
        cout << "Largest z-coordinate of BS is " << fdata_vtxzBSupper << endl;

        TFile data_plot("Histos/data_trees.root", "recreate");
        //gStyle->SetOptLogy();

        //canvas->Divide (2,2);

        //TCanvas *eta_canvas = new TCanvas ("eta canvas", "Eta Preliminary Results", 2);
        //eta_canvas->cd();
        data_eta_histo->Scale(1/fdata_trketa);
        data_eta_histo->SetMaximum(0.04);
        data_eta_histo->GetXaxis()->SetTitle("#eta");
        //data_eta_histo->GetYaxis()->SetTitleOffset(1.3);
        //data_eta_histo->GetYaxis()->SetTitleSize(0.01);
        data_eta_histo->GetYaxis()->SetTitle("#frac{1}{N_{ch}} #frac{dN_{ch}}{d#eta}");
        //data_eta_histo->GetYaxis()->SetLabelSize(0.02);
        //data_eta_histo->GetYaxis()->SetLabelOffset(0.001);
        //data_eta_histo->SetLineStyle(0);
        //data_eta_histo->SetLineColorAlpha(0, 0);
        //data_eta_histo->SetMarkerStyle(33);
        data_eta_histo->Draw("P");
        data_eta_histo->Write();
        //eta_canvas->Update();
        //eta_canvas->SaveAs("eta_data.png");

        //canvas->cd(3);
        //TCanvas *phi_canvas = new TCanvas ("phi canvas", "Phi Preliminary Results", );
        //phi_canvas->cd();
        data_phi_histo->Scale(1/fdata_trkphi);
        data_phi_histo->SetMinimum(0.014);
        data_phi_histo->SetMaximum(0.024);
        //data_phi_histo->SetLineStyle(0);
        //data_phi_histo->SetLineColorAlpha(0,0);
        //data_phi_histo->SetMarkerStyle(33);
        data_phi_histo->GetXaxis()->SetTitle("#phi");
        //data_phi_histo->GetYaxis()->SetTitleOffset(1.2);
        //data_phi_histo->GetYaxis()->SetTitleSize(0.01);
        //data_phi_histo->GetYaxis()->SetLabelOffset(0.001);
        //data_phi_histo->GetYaxis()->SetLabelSize(0.02);
        data_phi_histo->GetYaxis()->SetTitle("#frac{1}{N_{ch}} #frac{dN_{ch}}{d#phi}");
        data_phi_histo->Draw("P");
        data_phi_histo->Write();
        //phi_canvas->Update();
        //phi_canvas->SaveAs("phi_data.png");

        //canvas->cd(1);
        //TCanvas *pt_canvas = new TCanvas ("pt canvas", "pt Preliminary Results", 4);
        //pt_canvas->cd();
        //gPad->SetLogy();
        data_pt_histo->Scale(1/fdata_trkpt);
        data_pt_histo->GetXaxis()->SetTitle("p_{T}");
        //data_pt_histo->GetYaxis()->SetTitleOffset(1.3);
        data_pt_histo->GetYaxis()->SetTitle("Fraction of Tracks");
        //data_pt_histo->SetLineColorAlpha(0, 0);
        //data_pt_histo->SetMarkerStyle(33);
        data_pt_histo->Draw("P");
        data_pt_histo->Write();
        //pt_canvas->Update();
        //pt_canvas->SaveAs("pt_data.png");


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

        TCanvas *multiplicity = new TCanvas ("multiplicity", "Data Multiplicity");
        TPad *data_pad = new TPad ("data_pad", "Data Pad", 0, 0.3, 1.0, 1.0);
        //TPad *ratio_pad = new TPad ("ratio_pad", "Ratio Pad", 0, 0, 1.0, 0.3);

        data_pad->Draw();
        //ratio_pad->Draw();

        data_pad->cd();
        data_multiplicity->Scale(1/fdata_evt);
        data_multiplicity->GetXaxis()->SetTitle("N_{ch}");
        data_multiplicity->GetYaxis()->SetTitleOffset(1.1);
        data_multiplicity->GetYaxis()->SetTitle("P");
        data_multiplicity->SetLineColorAlpha(0,0);
        data_multiplicity->SetMarkerStyle(8);
        data_multiplicity->SetMarkerColor(kViolet);
        data_multiplicity->SetMarkerSize(1.5);
        data_multiplicity->SetLineWidth(2);
        data_multiplicity->SetLineColor(kGreen);
        gPad->SetTickx();
        gPad->SetTicky();
        data_multiplicity->Draw("P");
        data_multiplicity->Write();

        TLegend *dataleg_multiplicity = new TLegend (0.6, 0.6, 0.8, 0.8);
        dataleg_multiplicity->SetFillColor(0);
        dataleg_multiplicity->SetFillStyle(0);
        dataleg_multiplicity->SetBorderSize(0);
        dataleg_multiplicity->SetTextSize(0.04);
        dataleg_multiplicity->AddEntry(data_multiplicity, "Data", "lf");
        dataleg_multiplicity->Draw();
        dataleg_multiplicity->Write();

        data_pad->Modified();
        data_pad->Update();
        data_plot.Write();
        fdata_totalevt = fdata_evt;
    }
