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
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"

using namespace std;
using namespace ROOT::Math;

void Track_selection_2()
{

	//---------------------Histograms for reco Tracks-------------------------------------------

    TCanvas *reco_canvas = new TCanvas;
    TCanvas *gen_canvas = new TCanvas;
    //TH1F *event_histo = new TH1F ("reco_evt", "reco_evt", 100, 0, 200);
    TH1F *reco_pt_histo = new TH1F ("reco_pt", "Normalized_reco_pT", 99, 0.1, 10);
    TH1F *reco_eta_histo = new TH1F ("reco_eta", "Normalized_reco_Eta", 48, -2.4, 2.4);
    TH1F *reco_phi_histo = new TH1F ("reco_phi", "Normalized_reco_Phi", 62, -3.14, 3.14);
    TH1F *lumi_histo = new TH1F ("lumi_section", "lumi_section", 160, 80, 230);
    TH1F *reco_multiplicity = new TH1F ("reco_multiplicity", "Normalized_reco_Multiplicity", 50, 0, 50);
    TH1F *dz_sigmadz = new TH1F ("dz_sigmadz", "Normalized_dz_sigmadz", 200, -20, 20);
    TH1F *d0_sigmad0 = new TH1F ("d0_sigmad0", "Normalized_d0_sigmad0", 200, -20, 20);
    TH1F *sigmapt_pt = new TH1F ("sigmapt_pt", "Normalized_sigmapt_pt", 20, 0, 0.2);
	TH1F *norm_dz_sigmadz = new TH1F ("dz_sigmadz", "Normalized_dz_sigmadz", 200, -20, 20);
    TH1F *norm_d0_sigmad0 = new TH1F ("d0_sigmad0", "Normalized_d0_sigmad0", 200, -20, 20);
    TH1F *norm_sigmapt_pt = new TH1F ("sigmapt_pt", "Normalized_sigmapt_pt", 20, 0, 0.2);
    //TH1F *normalized_multiplicity_histo = new TH1F ("normalized_multiplicity", "normalized_multiplicity", 200, 0, 200);
	
	TH1F *MCpileup = new TH1F ("MCpileup", "MCpileup", 20, 0, 20);
	
	//---------------------Histograms for data Tracks-------------------------------------------

    TCanvas *data_reco_canvas = new TCanvas;
    TCanvas *data_gen_canvas = new TCanvas;
    TH1F *data_reco_pt_histo = new TH1F ("data_reco_pt", "data_Normalized_reco_pT", 99, 0.1, 10);
    TH1F *data_reco_eta_histo = new TH1F ("data_reco_eta", "data_Normalized_reco_Eta", 48, -2.4, 2.4);
    TH1F *data_reco_phi_histo = new TH1F ("data_reco_phi", "data_Normalized_reco_Phi", 62, -3.14, 3.14);
    TH1F *data_lumi_histo = new TH1F ("data_lumi_section", "data_lumi_section", 160, 80, 230);
    TH1F *data_reco_multiplicity = new TH1F ("data_reco_multiplicity", "data_Normalized_reco_Multiplicity", 50, 0, 50);
    TH1F *data_dz_sigmadz = new TH1F ("data_dz_sigmadz", "data_Normalized_dz_sigmadz", 200, -20, 20);
    TH1F *data_d0_sigmad0 = new TH1F ("data_d0_sigmad0", "data_Normalized_d0_sigmad0", 200, -20, 20);
    TH1F *data_sigmapt_pt = new TH1F ("data_sigmapt_pt", "data_Normalized_sigmapt_pt", 20, 0, 0.2);
	TH1F *data_norm_dz_sigmadz = new TH1F ("data_dz_sigmadz", "data_Normalized_dz_sigmadz", 200, -20, 20);
    TH1F *data_norm_d0_sigmad0 = new TH1F ("data_d0_sigmad0", "data_Normalized_d0_sigmad0", 200, -20, 20);
    TH1F *data_norm_sigmapt_pt = new TH1F ("data_sigmapt_pt", "data_Normalized_sigmapt_pt", 20, 0, 0.2);
    //TH1F *normalized_multiplicity_histo = new TH1F ("normalized_multiplicity", "normalized_multiplicity", 200, 0, 200);
	
	TH1F *data_pileup = new TH1F ("data_pileup", "data_pileup", 20, 0, 20);
	
	Double_t norm = 1;
	
	//RooUnfold
	
	RooUnfoldResponse response(reco_multiplicity, data_reco_multiplicity);
	
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

    vector<float> *fVec_reco_Vtxy = 0; //require initialisation to 0 to avoid crash
    tree->SetBranchAddress("vtxy",&fVec_reco_Vtxy);

	vector<float> *fVec_reco_Vtxx = 0; //require initialisation to 0 to avoid crash
    tree->SetBranchAddress("vtxx",&fVec_reco_Vtxx);
	
	vector<int> *fVec_reco_Vtxndof = 0;
    tree->SetBranchAddress("vtxndof",&fVec_reco_Vtxndof);
	
	vector<int> *fVec_reco_Vtxex = 0;
    tree->SetBranchAddress("vtxxErr",&fVec_reco_Vtxex);
	
	vector<int> *fVec_reco_Vtxey = 0;
    tree->SetBranchAddress("vtxyErr",&fVec_reco_Vtxey);
	
	vector<int> *fVec_reco_Vtxez = 0;
    tree->SetBranchAddress("vtxzErr",&fVec_reco_Vtxez);

    double dzleaf = 0;
	double dzcalc = 0;
	double dzminus= 0;

	double d0leaf = 0;
	double d0calc = 0;
	double d0minus = 0;


    //declare variable to hold track number
    int ntrk_reco_total;
    double d_reco_trk;
    int nreco_evt = 0;
	int ntrkd0 = 0;
	int ntrkdz = 0;
	int ntrkdpt = 0;
	int ntrkptcut= 0;
	int ntrketa = 0;
	int naccepttrk = 0;
    int ntrk_normalized = 0;
    int nHigh_Purity = 0;
	int npileup = 0;
    float dz_dzErr, d0_d0Err, ptErr_pt;



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
/*
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
*/

//---------------------Loops for reco tracks-------------------------------------------
				npileup++;
				
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
                            
							if((*fVec_reco_Vtxndof)[k] > 4)
							{
								++nreco_evt;
								
								//looping over tracks
								for (int j = 0; j != ntrk_reco_total; ++j)
								{
									XYZTVector reco_vec = (*reco_tracks)[j];
								
									if((*nVec_HighPurity)[j] == 1)
									{

										//dzcalc = ((reco_vec.Z())-(*fVec_reco_Vtxz)[k])-((((reco_vec.X())-(*fVec_reco_Vtxx)[k])*(reco_vec.px())+((reco_vec.Y())-(*fVec_reco_Vtxy)[k])*(reco_vec.py()))/(reco_vec.Pt())*(reco_vec.pz()/reco_vec.Pt()));
										//d0calc = ((-(reco_vec.X() - (*fVec_reco_Vtxx)[k])*reco_vec.Py()) + ((reco_vec.Y()-(*fVec_reco_Vtxy)[k])*reco_vec.Px()))/reco_vec.Pt();
										
										//fabs(tr.vtxdxy.at(vtxnum) / sqrt(pow(tr.ed0, 2) + pow(goodVtx->ex, 2) + pow(goodVtx->ey, 2))) > IPsig
										//fabs(tr.vtxdz.at(vtxnum) / sqrt(pow(tr.edz, 2) + pow(goodVtx->ez, 2))) > IPsig
										dz_dzErr = (((*fVec_dz)[j])/sqrt(pow(((*fVec_dzErr)[j]),2)+pow((*fVec_reco_Vtxez)[k],2)));
										d0_d0Err = (((*fVec_d0)[j])/sqrt(pow(((*fVec_d0Err)[j]),2)+pow((*fVec_reco_Vtxex)[k],2)+pow((*fVec_reco_Vtxey)[k],2)));
										//dz_dzErr = (dzcalc)/((*fVec_dzErr)[j]);
										//d0_d0Err = (d0calc)/((*fVec_d0Err)[j]);
										//dz_dzErr = ((*fVec_dz)[j])/((*fVec_dzErr)[j]);
										//d0_d0Err = ((*fVec_d0)[j])/((*fVec_d0Err)[j]);
										ptErr_pt = ((*fVec_ptErr)[j]/(reco_vec.Pt()));
										pt_efficiency->Fill(fpt_eff);
										eta_efficiency->Fill(feta_eff);
										phi_efficiency->Fill(fphi_eff);
									
										if (abs(dz_dzErr) <= 3 && abs(ptErr_pt) <= 0.05 && abs(reco_vec.Eta()) <= 2.4 && reco_vec.Pt() >= 0.5)
										{
										d0_sigmad0->Fill(d0_d0Err);
										++ntrkd0;
										}
										if (abs(reco_vec.Eta()) <= 2.4 && reco_vec.Pt() >= 0.5 && abs(d0_d0Err) <= 3 && abs(ptErr_pt) <= 0.05 )
										{
										dz_sigmadz->Fill(dz_dzErr);
										++ntrkdz;
										}
										if (abs(reco_vec.Eta()) <= 2.4 && reco_vec.Pt() >= 0.5 && abs(d0_d0Err) <= 3 && abs(dz_dzErr) <= 3 )
										{
										sigmapt_pt->Fill(ptErr_pt);
										++ntrkdpt;
										}	
     
										//trackselection conditions to go here
										if (abs(reco_vec.Eta()) <= 2.4 && reco_vec.Pt() >= 0.5 && abs(d0_d0Err) <= 3 && abs(dz_dzErr) <= 3 && abs(ptErr_pt) <= 0.05 )
										{
											reco_phi_histo->Fill(reco_vec.Phi());
											++naccepttrk;
											++d_reco_trk;
										}
										if (abs(reco_vec.Eta()) <= 2.4 && abs(d0_d0Err) <= 3 && abs(dz_dzErr) <= 3 && abs(ptErr_pt) <= 0.05 )
										{	
											reco_pt_histo->Fill(reco_vec.Pt());
											++ntrkptcut;
										}
										//if (reco_vec.Pt() >= 0.5 && abs(d0_d0Err) <= 3 && abs(dz_dzErr) <= 3 && abs(ptErr_pt) <= 0.05 )
										if (reco_vec.Pt() >= 0.5 && abs(ptErr_pt) <= 0.05 )											
										{
											reco_eta_histo->Fill(reco_vec.Eta());
											++ntrketa;
										}
									}
								}

								reco_multiplicity->Fill(d_reco_trk);
								ntrk_normalized += d_reco_trk;
									
							}
							



                        }


                    }
					
                }
				
            }
			
			MCpileup->Fill(fVec_reco_Vtxz_size);
        }

    }

    //cout << "Total number of selected tracks is " << ntrk_normalized << endl;
    //cout << "d0leaf is " << d0leaf << " d0calc is "<< d0calc << endl;
    //cout << "Total number of selected events is " << d_reco_evt << endl;
    //cout << "Total number of selected tracks is " << d_reco_trk << endl;

    //vtxz_plot->Draw();
    //reco_canvas->Update();
    //reco_canvas->SaveAs("vtxz_number.png");
    //delete reco_canvas;
    //delete vtxz_plot;


    //after the loop over all events, draw the resulting plots

	double norm_evt = (1./nreco_evt);
	double norm_accepttrk = (1./naccepttrk);
	cout << "normevt is " << norm_evt << "nreco_evt" << nreco_evt<<endl;
	
	
    reco_canvas->Divide(3,2);
	
    reco_canvas->cd(1);
	//double norm_dz = (1./ntrkdz);
	//dz_sigmadz->Scale(norm_dz);
	dz_sigmadz->Scale(norm/dz_sigmadz->Integral(), "width");
	gPad->SetLogy();
	dz_sigmadz->Draw();
	
    //gPad->SetLogx();
    //reco_pt_histo->Scale(1/n_reco_evt);
    //reco_pt_histo->Draw();
    //gen_pt_histo->SetLineColor(kRed-2);
    //gen_pt_histo->Scale(1/n_reco_evt, "SAME");
    //gen_pt_histo->Draw();

    reco_canvas->cd(2);
	//double norm_d0 = (1./ntrkd0);
	//d0_sigmad0->Scale(norm_d0);
	d0_sigmad0->Scale(norm/d0_sigmad0->Integral(), "width");
	gPad->SetLogy();
	d0_sigmad0->Draw();
	
    //gPad->SetLogy();
    //gPad->SetLogx();
    //reco_eta_histo->Scale(1/n_reco_evt);
    //reco_eta_histo->Draw();
    //gen_eta_histo->SetLineColor(kBlue-4);
    //gen_eta_histo->Scale(1/n_reco_evt, "SAME");
    //gen_eta_histo->Draw();

    reco_canvas->cd(3);
    //double norm_dpt = (1./ntrkdpt);
	//sigmapt_pt->Scale(norm_dpt);
	sigmapt_pt->Scale(norm/sigmapt_pt->Integral(), "width");
	gPad->SetLogy();
	sigmapt_pt->Draw();
	
	reco_canvas->cd(4);
	//double norm_ptcut = (1./ntrkptcut);
    //reco_pt_histo->Scale(norm_ptcut);
	reco_pt_histo->Scale(norm/reco_pt_histo->Integral(), "width");
	gPad->SetLogy();
	reco_pt_histo->Draw();
	
	reco_canvas->cd(5);
	//double norm_eta = (1./ntrketa);
    //reco_eta_histo->Scale(norm_eta);
	reco_eta_histo->Scale(norm/reco_eta_histo->Integral(), "width");
	reco_eta_histo->Draw();
	
	reco_canvas->cd(6);
    //reco_phi_histo->Scale(norm_accepttrk);
	reco_phi_histo->Scale(norm/reco_phi_histo->Integral(), "width");
	reco_phi_histo->Draw();
	
	reco_multiplicity->Scale(norm/reco_multiplicity->Integral(), "width");
	//reco_multiplicity->Scale(norm_evt);
	
	
	//gPad->SetLogy();
    //reco_phi_histo->Scale(1/n_reco_evt);
    //reco_phi_histo->Draw();
    //gen_phi_histo->SetLineColor(kGreen);
    //gen_phi_histo->Scale(1/n_reco_evt, "SAME");
    //gen_phi_histo->Draw();

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
    reco_canvas->SaveAs("Track_selection.pdf");

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
	
	//running over data
	{   
   //gROOT->ProcessLine(".L Loader.C+");

    // attach "UETree/data" of tree1.root as the main root file for this program
    TFile *data_File = TFile::Open("tree1.root", "READ");
    TTree* data_tree = (TTree*)data_File->Get("UETree/data");

    //for this dataset we want lumisection of 90 and above
    //data from branch 62
    int data_nlumi_section;
    data_tree->SetBranchAddress("lumi", &data_nlumi_section);

    int data_iZeroBias; //data from Branch 60
    data_tree->SetBranchAddress("trgZeroBias",&data_iZeroBias);

//---------------------Variables for reco tracks-------------------------------------------

    //variables to check if there is only 1 vertex
    vector<float> *data_fVec_reco_Vtxz = 0; //require initialisation to 0 to avoid crash
    data_tree->SetBranchAddress("vtxz",&data_fVec_reco_Vtxz);
    float data_fVec_reco_Vtxz_size;

    vector<float> *data_fVec_reco_VtxzBS = 0; //from branch 54
    data_tree->SetBranchAddress("vtxzBS", &data_fVec_reco_VtxzBS);
    float data_n_reco_BeamSize;

    vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *data_reco_tracks = 0;
    data_tree->SetBranchAddress("recoTracksp4", &data_reco_tracks);
    //vector<TLorentzVector> *tracks = new vector<TLorentzVector>;

    vector<float> *data_fVec_dz = 0;
    data_tree->SetBranchAddress ("recoTracksdz", &data_fVec_dz);

    vector<float> *data_fVec_dzErr = 0;
    data_tree->SetBranchAddress ("recoTracksdzErr", &data_fVec_dzErr);

    vector<float> *data_fVec_d0 = 0;
    data_tree->SetBranchAddress ("recoTracksd0", &data_fVec_d0);
    vector<float> *data_fVec_d0Err = 0;
    data_tree->SetBranchAddress ("recoTracksd0Err", &data_fVec_d0Err);

    vector<float> *data_fVec_ptErr = 0;
    data_tree->SetBranchAddress("recoTracksptErr", &data_fVec_ptErr);

    vector<int> *data_nVec_HighPurity = 0;
    data_tree->SetBranchAddress ("recoTrackshighPurity", &data_nVec_HighPurity);

    vector<float> *data_fVec_reco_Vtxy = 0; //require initialisation to 0 to avoid crash
    data_tree->SetBranchAddress("vtxy",&data_fVec_reco_Vtxy);

	vector<float> *data_fVec_reco_Vtxx = 0; //require initialisation to 0 to avoid crash
    data_tree->SetBranchAddress("vtxx",&data_fVec_reco_Vtxx);
	
	vector<int> *data_fVec_reco_Vtxndof = 0;
    data_tree->SetBranchAddress("vtxndof",&data_fVec_reco_Vtxndof);

    double data_dzleaf = 0;
	double data_dzcalc = 0;
	double data_dzminus= 0;

	double data_d0leaf = 0;
	double data_d0calc = 0;
	double data_d0minus = 0;


    //declare variable to hold track number
    int data_ntrk_reco_total;
    double data_d_reco_trk;
    int data_nreco_evt = 0;
	int data_naccepttrk = 0;
	int data_ntrkd0 = 0;
	int data_ntrkdz = 0;
	int data_ntrkdpt = 0;
	int data_ntrkptcut= 0;
	int data_ntrketa = 0;
    int data_ntrk_normalized = 0;
    int data_nHigh_Purity = 0;
	int data_npileup = 0;
    float data_dz_dzErr, data_d0_d0Err, data_ptErr_pt;



//---------------------Cuts-------------------------------------------

    const int data_lumi_cut = 90;
    const double data_eta_cut = 2.4;
    const double data_pt_cut = 0.5;
    const double data_vtxz_number = 1.;
    const double data_vtxz_size = 10;

    //select events with only 1 vertex using information from Branch 46/47/48. Here we use Branch 48 (z-axis)
    Int_t data_nEvt = (Int_t)data_tree->GetEntries();
    cout << "data_Tree Entries " << data_nEvt << endl;


    //start loop over events
    for(Int_t i = 0; i < data_nEvt; ++i)
    {
        //if (i%10000 == 0)
            //cout<< "Processing event number: " << i << endl; //1 TTree entry = 1 event
        cout << "Entry " << i << endl;
        data_tree -> GetEntry(i);

        //selects for events from lumisection >= 90
        if (data_nlumi_section >= data_lumi_cut)
        {
            cout << "lumisection is " << data_nlumi_section << endl;


            //we select only events that are triggered by ZeroBias Trigger. "True" in ZeroBias evaluate to 1
            if (data_iZeroBias == 1)
            {
                cout << "track pass zero bias" << endl;
                data_fVec_reco_Vtxz_size = data_fVec_reco_Vtxz->size();
                //vtxz_plot->Fill(fvecVtxz_size);

//---------------------Loops for reco tracks-------------------------------------------
				data_npileup++;
                if (data_fVec_reco_Vtxz_size == data_vtxz_number)
                {
                    //cout << "number of vertex for event " << i << " is " << fVec_reco_Vtxz_size << endl;

                    //looping over vertices
                    for (int k = 0; k != data_fVec_reco_Vtxz_size; ++k)
                    {
                        data_n_reco_BeamSize = fabs((*data_fVec_reco_Vtxz)[k] - (*data_fVec_reco_VtxzBS)[0]);
                        cout << "Beam Size is " << data_n_reco_BeamSize << endl;

                        if (data_n_reco_BeamSize <= data_vtxz_size)
                        {
                            data_ntrk_reco_total = data_reco_tracks->size();

                            //fill with the lumi sections that meet the above event-level cuts
                            data_lumi_histo->Fill(data_nlumi_section);

                            data_d_reco_trk = 0;
                            
							if((*data_fVec_reco_Vtxndof)[k] > 4)
							{
							
								++data_nreco_evt;


								//looping over tracks
								for (int j = 0; j != data_ntrk_reco_total; ++j)
								{

									XYZTVector data_reco_vec = (*data_reco_tracks)[j];

									if((*data_nVec_HighPurity)[j] == 1)
									{
										data_dzcalc = ((data_reco_vec.Z())-(*data_fVec_reco_Vtxz)[k])-((((data_reco_vec.X())-(*data_fVec_reco_Vtxx)[k])*(data_reco_vec.px())+((data_reco_vec.Y())-(*data_fVec_reco_Vtxy)[k])*(data_reco_vec.py()))/(data_reco_vec.Pt())*(data_reco_vec.pz()/data_reco_vec.Pt()));
										data_d0calc = ((-(data_reco_vec.X() - (*data_fVec_reco_Vtxx)[k])*data_reco_vec.Py()) + ((data_reco_vec.Y()-(*data_fVec_reco_Vtxy)[k])*data_reco_vec.Px()))/data_reco_vec.Pt();
										data_dz_dzErr = (data_dzcalc)/((*data_fVec_dzErr)[j]);
										data_d0_d0Err = (data_d0calc)/((*data_fVec_d0Err)[j]);
										//data_dz_dzErr = ((*data_fVec_dz)[j])/((*data_fVec_dzErr)[j]);
										//data_d0_d0Err = ((*data_fVec_d0)[j])/((*data_fVec_d0Err)[j]);
										data_ptErr_pt = ((*data_fVec_ptErr)[j]/(data_reco_vec.Pt()));

										if (abs(data_reco_vec.Eta()) <= 2.4 && data_reco_vec.Pt() >= 0.5 && abs(data_dz_dzErr) <= 3 && abs(data_ptErr_pt) <= 0.05)
										{
										data_d0_sigmad0->Fill(data_d0_d0Err);
										++data_ntrkd0;
										}
										if (abs(data_reco_vec.Eta()) <= 2.4 && data_reco_vec.Pt() >= 0.5 && abs(data_d0_d0Err) <= 3 && abs(data_ptErr_pt) <= 0.05)
										{
										data_dz_sigmadz->Fill(data_dz_dzErr);
										++data_ntrkdz;
										}
										if (abs(data_reco_vec.Eta()) <= 2.4 && data_reco_vec.Pt() >= 0.5 && abs(data_d0_d0Err) <= 3 && abs(data_dz_dzErr) <= 3 )
										{
										data_sigmapt_pt->Fill(data_ptErr_pt);
										++data_ntrkdpt;
										}	

										//trackselection conditions to go here
										if (abs(data_reco_vec.Eta()) <= 2.4 && data_reco_vec.Pt() >= 0.5 && abs(data_d0_d0Err) <= 3 && abs(data_dz_dzErr) <= 3 && abs(data_ptErr_pt) <= 0.05)
										{
											data_reco_phi_histo->Fill(data_reco_vec.Phi());
											++data_naccepttrk;
											++data_d_reco_trk;
										}
									
										if (abs(data_reco_vec.Eta()) <= 2.4 && abs(data_d0_d0Err) <= 3 && abs(data_dz_dzErr) <= 3 && abs(data_ptErr_pt) <= 0.05)
										{
											data_reco_pt_histo->Fill(data_reco_vec.Pt());
											++data_ntrkptcut;
										}										
									
										if (data_reco_vec.Pt() >= 0.5 && abs(data_d0_d0Err) <= 3 && abs(data_dz_dzErr) <= 3 && abs(data_ptErr_pt) <= 0.05)
										{
											data_reco_eta_histo->Fill(data_reco_vec.Eta());
											++data_ntrketa;
										}										
									}
								}

								data_reco_multiplicity->Fill(data_d_reco_trk);
								data_ntrk_normalized += data_d_reco_trk;
							}
                        }


                    }

                }

            }
			data_pileup->Fill(data_fVec_reco_Vtxz_size);
        }

    }

    //cout << "Total number of selected tracks is " << ntrk_normalized << endl;
    //cout << "d0leaf is " << d0leaf << " d0calc is "<< d0calc << endl;
    //cout << "Total number of selected events is " << d_reco_evt << endl;
    //cout << "Total number of selected tracks is " << d_reco_trk << endl;

    //vtxz_plot->Draw();
    //reco_canvas->Update();
    //reco_canvas->SaveAs("vtxz_number.png");
    //delete reco_canvas;
    //delete vtxz_plot;


    //after the loop over all events, draw the resulting plots

	double data_norm_evt = (1./data_nreco_evt);
	double data_norm_accepttrk = (1./data_naccepttrk);
	cout << "data_nreco_evt" << data_nreco_evt << "data_ntrkdz    " << data_ntrkdz <<  "data_ntrkd0   " << data_ntrkd0 << "data_ntrkdpt   " << data_ntrkdpt << "data_naccepttrk   " << data_naccepttrk << "data_ntrketa   " << data_ntrketa << "data_ntrkptcut   " << data_ntrkptcut << endl;

    data_reco_canvas->Divide(3,2);
	
    data_reco_canvas->cd(1);
	//double data_norm_dz = (1./data_ntrkdz);
	//data_dz_sigmadz->Scale(data_norm_dz);
	data_dz_sigmadz->Scale(norm/data_dz_sigmadz->Integral(), "width");
	gPad->SetLogy();
	data_dz_sigmadz->Draw();
	
    //gPad->SetLogx();
    //reco_pt_histo->Scale(1/n_reco_evt);
    //reco_pt_histo->Draw();
    //gen_pt_histo->SetLineColor(kRed-2);
    //gen_pt_histo->Scale(1/n_reco_evt, "SAME");
    //gen_pt_histo->Draw();

    data_reco_canvas->cd(2);
	//double data_norm_d0 = (1./data_ntrkd0);
	//data_d0_sigmad0->Scale(data_norm_d0);
	data_d0_sigmad0->Scale(norm/data_d0_sigmad0->Integral(), "width");
	gPad->SetLogy();
	data_d0_sigmad0->Draw();
	
    //gPad->SetLogy();
    //gPad->SetLogx();
    //reco_eta_histo->Scale(1/n_reco_evt);
    //reco_eta_histo->Draw();
    //gen_eta_histo->SetLineColor(kBlue-4);
    //gen_eta_histo->Scale(1/n_reco_evt, "SAME");
    //gen_eta_histo->Draw();

    data_reco_canvas->cd(3);
	//double data_norm_dpt = (1./data_ntrkdpt);
    //data_sigmapt_pt->Scale(data_norm_dpt);
	data_sigmapt_pt->Scale(norm/data_sigmapt_pt->Integral(), "width");
	gPad->SetLogy();
	data_sigmapt_pt->Draw();
	
	data_reco_canvas->cd(4);
	//double data_norm_ptcut = (1./data_ntrkptcut);
    //data_reco_pt_histo->Scale(data_norm_ptcut);
	data_reco_pt_histo->Scale(norm/data_reco_pt_histo->Integral(), "width");
	gPad->SetLogy();
	data_reco_pt_histo->Draw();
	
	data_reco_canvas->cd(5);
	//double data_norm_eta = (1./data_ntrketa);	
    //data_reco_eta_histo->Scale(data_norm_eta);
	data_reco_eta_histo->Scale(norm/data_reco_eta_histo->Integral(), "width");
	data_reco_eta_histo->Draw();
	
	data_reco_canvas->cd(6);
    //data_reco_phi_histo->Scale(data_norm_accepttrk);
	data_reco_phi_histo->Scale(norm/data_reco_phi_histo->Integral(), "width");
	data_reco_phi_histo->Draw();
	
	data_reco_multiplicity->Scale(norm/data_reco_multiplicity->Integral(), "width");
	//data_reco_multiplicity->Scale(data_norm_evt);
	
	//gPad->SetLogy();
    //reco_phi_histo->Scale(1/n_reco_evt);
    //reco_phi_histo->Draw();
    //gen_phi_histo->SetLineColor(kGreen);
    //gen_phi_histo->Scale(1/n_reco_evt, "SAME");
    //gen_phi_histo->Draw();

    /*reco_canvas->cd(4);
    gPad->SetLogy();
    dz_sigmadz->DrawNormalized("", 1);

    reco_canvas->cd(5);
    gPad->SetLogy();
    d0_sigmad0->DrawNormalized("", 1);

    reco_canvas->cd(6);
    gPad->SetLogy();
    sigmapt_pt->DrawNormalized("", 1);*/
    data_reco_canvas->Update();
    data_reco_canvas->SaveAs("data_Track_selection.pdf");

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
	
	
	
	//Writing output
	TFile f("Track_selection_2_nocut.root", "recreate");
	d0_sigmad0->Write();
	dz_sigmadz->Write();
	sigmapt_pt->Write();
	reco_pt_histo->Write();
	reco_eta_histo->Write();
	reco_phi_histo->Write();
	reco_multiplicity->Write();
	MCpileup->Write();
	
	data_d0_sigmad0->Write();
	data_dz_sigmadz->Write();
	data_sigmapt_pt->Write();
	data_reco_pt_histo->Write();
	data_reco_eta_histo->Write();
	data_reco_phi_histo->Write();
	data_reco_multiplicity->Write();
	data_pileup->Write();
	
}