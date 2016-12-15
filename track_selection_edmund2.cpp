#include "TROOT.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
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

#ifdef __MAKECINT__
// // #pragma link C++ class vector<TLorentzVector>++
// // #pragma link C++ class vector<vector<TLorentzVector>>++
// // #pragma link C++ class ROOT::Math::PxPyPzE4D<double>++
// #pragma link C++ class ROOT::Math::PxPyPzE4D<double>++;
// #pragma link C++ class ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >++;
#pragma link C++ class vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >++;
#endif


void edmund2()
{
    // ===========================implement cuts================================

    // const float dz_dzErr_cut = 3;
	const int vtx_number_cut = 1;
    const int lumi_cut = 90;

    // ==============================================Variables==============================================================

    int ndata_numberofvtxx = 0;
	int	ndata_totalvtxx = 0;
	int fdata_numselectedvtxx = 0;
    int ndata_numberoftrk = 0;
    int ndata_totaltrk = 0;
	int ndata_numselectedtrk = 0;
	int ndata_numtrkupper = 0;
	float d0_calc = 0;
	float dz_calc = 0;
	float d0err = 0;
	float dzerr = 0;
	float d0_d0err = 0;
	float dz_dzerr = 0;
	float fdata_sigmapt_pt = 0;
	float dz_track = 0;
	float d0_qx=0;
	


    // ===========================retrieve ROOT file============================

    // TFile *datafile = TFile::Open("EdmundHerwigtrees_99.root", "READ");
    // TFile *datafile = TFile::Open("EdmundZBtrees_1.root", "READ");
    TFile *datafile = TFile::Open("treesCUETP8M1_66.root", "READ");
    // TFile *datafile = TFile::Open("EdmundZB7trees_8.root", "READ");
    TTree *datatree = (TTree*)datafile->Get("UETree/data");

    // ===========================define variables to read TTree================

    vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *data_tracks = 0;
    datatree->SetBranchAddress("recoTracksp4", &data_tracks);
	
    vector<float> *fvecdata_vtxx = 0;
    datatree->SetBranchAddress("vtxx", &fvecdata_vtxx);

    int ndata_lumi;
    datatree->SetBranchAddress("lumi", &ndata_lumi);

    int ndata_zerobias;
    datatree->SetBranchAddress("trgZeroBias", &ndata_zerobias);

    vector<float> *fvecdata_vtxy = 0;
    datatree->SetBranchAddress("vtxy", &fvecdata_vtxy);

    vector<float> *fvecdata_vtxz = 0;
    datatree->SetBranchAddress("vtxz", &fvecdata_vtxz);
	
	//new 


		// vector<float> *fvecdata_vtxxBS = 0;
		// datatree->SetBranchAddress("vtxxBS", &fvecdata_vtxxBS);

		// vector<float> *fvecdata_vtxyBS = 0;
		// datatree->SetBranchAddress("vtxyBS", &fvecdata_vtxyBS);

		// vector<float> *fvecdata_vtxzBS = 0;
		// datatree->SetBranchAddress("vtxzBS", &fvecdata_vtxzBS);

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

		vector<int> *nvecdata_validhits = 0;
		datatree->SetBranchAddress("recoTracksnValidHits", &nvecdata_validhits);

		vector<float> *fvecdata_trackschi2n = 0;
		datatree->SetBranchAddress("recoTrackschi2n", &fvecdata_trackschi2n);


// ============================================== Histos for vtx ==============================================================

    TH1F *data_vtxx = new TH1F ("vtxx", "vtxx", 400, -0.2, 0.4);
    TH1F *data_vtxy = new TH1F ("vtxy", "vtxy", 400, -0.2, 0.4);
    TH1F *data_vtxz = new TH1F ("vtxz", "vtxz", 400, -20, 20);
    TH1F *data_vtxperevent = new TH1F ("vtxperevent", "vtx per event", 50, 0, 50);

    TH1F *data_norm_vtxx = new TH1F ("norm_vtxx", "normalised vtxx", 400, -0.2, 0.4);
    TH1F *data_norm_vtxy = new TH1F ("norm_vtxy", "normalised vtxy", 400, -0.2, 0.4);
    TH1F *data_norm_vtxz = new TH1F ("norm_vtxz", "normalised vtxz", 400, -20, 20);
    TH1F *data_norm_vtxperevent = new TH1F ("norm_vtxperevent", "normalised vtx per event", 50, 0, 50);
	

	

// ============================================== Histos for trk ==============================================================

    TH1F *data_trkperevent = new TH1F ("trkperevent", "normalised trk per event", 400, 0, 400);
    TH1F *data_pt_histo = new TH1F ("data pT", "Normalized data p_{T}", 200, 0, 10);
    TH1F *data_eta_histo = new TH1F("data eta", "Normalised data #eta", 50, -2.6, 2.6);
    TH1F *data_phi_histo = new TH1F ("data phi", "Normalized_data #phi", 60, -4, 4);
	TH1F *data_d0_d0err = new TH1F ("d0_d0err", "d0_d0err", 200, -20, 20);
	TH1F *data_dz_dzerr = new TH1F ("dz_dzerr", "dz_dzerr", 200, -20, 20);
	TH1F *data_d0 = new TH1F ("data_d0", "data_d0", 200, -5, 5);
	TH1F *data_d0_qx = new TH1F ("data_d0_qx", "data_d0_qx", 200, -5, 5);	
	TH1F *data_dz = new TH1F ("data_dz", "data_dz", 200, -20, 20);
	TH1F *data_d0_err = new TH1F ("data_d0_err", "data_d0_err", 200, 0, 5);
	TH1F *data_dz_err = new TH1F ("data_dz_err", "data_dz_err", 200, 0, 5);
	TH1F *data_dz_track = new TH1F ("data_dz_track", "data_dz_track", 200, -20, 20);
	
	
	TH1F *data_sigmapt_pt = new TH1F ("data_sigmapt_pt", "data_sigmapt_pt", 20, 0, 0.2);
	
	
	
// ============================================== Histos for events with 1 vtx ==============================================================
	
	TH1F *data_trk_singlevtx = new TH1F ("trk_singlevtx", "normalised trk per vtx", 400, 0, 400);
	
    Int_t ndata_totalEvt = (Int_t)datatree->GetEntries();
    cout << "There is a total of " << ndata_totalEvt << " events." << endl;


// ========================================================= Start of Evt Loop ================================================================

    for (Int_t i = 0; i < ndata_totalEvt; ++i)
    {
        datatree->GetEntry(i);
        // cout << "At entry " << i << endl;
		
		if (ndata_lumi >= lumi_cut)
        {
            if (ndata_zerobias == 1)
            {
				
				ndata_numberofvtxx = fvecdata_vtxx->size();
				ndata_numberoftrk = data_tracks->size();
				
				ndata_totalvtxx += ndata_numberofvtxx ;
				ndata_totaltrk += ndata_numberoftrk ;
				
				if (ndata_numtrkupper < ndata_numberoftrk)
				{
					ndata_numtrkupper = ndata_numberoftrk;
				}
				
				data_vtxperevent->Fill(ndata_numberofvtxx);
				data_norm_vtxperevent->Fill(ndata_numberofvtxx);
				data_trkperevent->Fill(ndata_numberoftrk);
		

// ========================================================= Start of Vertex Loop ================================================================

				for (int vtxnumber = 0; vtxnumber < ndata_numberofvtxx; ++vtxnumber)
				{

					data_vtxx->Fill((*fvecdata_vtxx)[vtxnumber]);
					data_vtxy->Fill((*fvecdata_vtxy)[vtxnumber]);
					data_vtxz->Fill((*fvecdata_vtxz)[vtxnumber]);
					
					data_norm_vtxx->Fill((*fvecdata_vtxx)[vtxnumber]);
					data_norm_vtxy->Fill((*fvecdata_vtxy)[vtxnumber]);
					data_norm_vtxz->Fill((*fvecdata_vtxz)[vtxnumber]);
					
				// if (ndata_numberofvtxx == vtx_number_cut && (*nvecdata_vtxndof)[0] > 4 && (*fvecdata_vtxzerr)[0] < 4)
				if (ndata_numberofvtxx == vtx_number_cut && (*nvecdata_vtxndof)[0] > 4)
					{
						data_trk_singlevtx->Fill(ndata_numberoftrk);				
						++fdata_numselectedvtxx;
								
// ========================================================= Start of Trk Loop ================================================================

						for (int t = 0; t < ndata_numberoftrk; ++t)
						{
							XYZTVector data_vec = (*data_tracks)[t];
							
							data_pt_histo->Fill(data_vec.Pt());
							data_eta_histo->Fill(data_vec.Eta());
							data_phi_histo->Fill(data_vec.Phi());
							
							fdata_sigmapt_pt = (((*fvecdata_pterr)[t])/(data_vec.Pt())); 
							data_sigmapt_pt->Fill(fdata_sigmapt_pt);
							
							if(fabs(data_vec.Eta())< 2.4)
							// if(fabs(data_vec.Eta())< 2.4 && (*nvecdata_highpurity)[t] == 1 && (data_vec.Pt())>0.5 && fdata_sigmapt_pt < 0.05)
							{
								d0_qx = sqrt(pow((data_vec.x()-(*fvecdata_vtxx)[vtxnumber]),2)+pow((data_vec.y()-(*fvecdata_vtxy)[vtxnumber]),2));
								d0_calc = (-(data_vec.x()-(*fvecdata_vtxx)[vtxnumber])*(data_vec.py())+(data_vec.y()-(*fvecdata_vtxy)[vtxnumber])*(data_vec.px()))/(data_vec.pt());
								dz_calc = (data_vec.z()-(*fvecdata_vtxz)[vtxnumber])-((data_vec.x()-(*fvecdata_vtxx)[vtxnumber])*(data_vec.px())+(data_vec.y()-(*fvecdata_vtxy)[vtxnumber])*(data_vec.py()))*(data_vec.pz())/pow((data_vec.pt()),2);  
								d0err = sqrt(pow(((*fvecdata_d0err)[t]),2)+pow(((*fvecdata_vtxxerr)[vtxnumber]),2)+pow(((*fvecdata_vtxyerr)[vtxnumber]),2));
								dzerr = sqrt(pow(((*fvecdata_dzerr)[t]),2)+pow(((*fvecdata_vtxzerr)[vtxnumber]),2));
								dz_track = (data_vec.z()-(*fvecdata_vtxz)[vtxnumber]);
								data_dz_track->Fill(dz_track);
								d0_d0err = d0_qx/d0err;
								dz_dzerr = dz_calc/dzerr;
								data_d0_qx->Fill(d0_qx);
								data_d0->Fill(d0_calc);
								data_dz->Fill(dz_calc);
								data_d0_err->Fill(d0err);
								data_dz_err->Fill(dzerr);
								// if(fabs(dz_dzerr)<3)
								// {
									data_d0_d0err->Fill(d0_d0err);
								// }
								
								// if(fabs(d0_d0err)<3)
								// {
									data_dz_dzerr->Fill(dz_dzerr);
								// }
								
								++ndata_numselectedtrk;
							}
							
						}

					}
						
// ========================================================= End of Trk Loop ================================================================

				}
// ========================================================= End of Vertex Loop ================================================================
 
			}
		}
    }
// ========================================================= End of Evt Loop ================================================================

    TFile data_plot("HistosEdmundHerwigtrees_99.root", "recreate");
    // TFile data_plot("HistosEdmundZBtrees_1.root", "recreate");
    // TFile data_plot("HistostreesCUETP8M1_66.root", "recreate");
    // TFile data_plot("HistosEdmundZB7trees_8.root", "recreate");
    //TCanvas *canvas = new TCanvas ("foo bar", "foo #bar");

    cout << "ndata_totalvtxx is " << ndata_totalvtxx << endl;
    cout << "ndata_numselectedtrk is " << ndata_numselectedtrk << endl;
    cout << "fdata_numselectedvtxx is " << fdata_numselectedvtxx << endl;
    cout << "ndata_numtrkupper is " << ndata_numtrkupper << endl;
	
	data_norm_vtxx->Scale(1./ndata_totalvtxx);
    data_norm_vtxy->Scale(1./ndata_totalvtxx);
    data_norm_vtxz->Scale(1./ndata_totalvtxx);
    data_norm_vtxperevent->Scale(1./ndata_totalEvt);
    data_trkperevent->Scale(1./ndata_totalEvt);
    data_pt_histo->Scale(1./ndata_numselectedtrk);
    data_eta_histo->Scale(1./ndata_numselectedtrk);
    data_phi_histo->Scale(1./ndata_numselectedtrk);
    data_trk_singlevtx->Scale(1./fdata_numselectedvtxx);

    data_vtxx->Write();
    data_vtxy->Write();
    data_vtxz->Write();	
	data_vtxperevent->Write();
    data_norm_vtxx->Write();
    data_norm_vtxy->Write();
    data_norm_vtxz->Write();	
	data_norm_vtxperevent->Write();
	data_trkperevent->Write();
	data_pt_histo->Write();
    data_eta_histo->Write();
    data_phi_histo->Write();
	data_trk_singlevtx->Write();
	data_d0_d0err->Write();
	data_dz_dzerr->Write();
	data_d0->Write();
	data_dz->Write();
	data_d0_err->Write();
	data_dz_err->Write();
	data_dz_track->Write();
	data_d0_qx->Write();
	
	data_sigmapt_pt->Write();

    data_plot.Write();

}

void track_selection_edmund2()
{
	// Write histograms to file.
	edmund2();
	
	// TFile *f1 = new TFile("HistosEdmundZB7trees_8.root");
	// TH1F *h1; f1->GetObject("trk_singlevtx;1", h1);
	// TFile *f2 = new TFile("HistosEdmundHerwigtrees_99.root");
	// TH1F *h2; f2->GetObject("trk_singlevtx;1", h2);
	
	// TCanvas *c1 = new TCanvas ("c1", "Superimpose two histograms");
	
	// c1->SetLogy();
	// h1->SetLineColor(4);
	// h1->Draw();
	// h2->SetLineColor(2);
	// h2->Draw("SAME");
	
	// TLegend *leg = new TLegend(0.1,0.1,0.48,0.3);
	// leg->AddEntry(h1,"CMS data","l");
	// leg->AddEntry(h2,"Monte Carlo","l");
	// leg->Draw();
	// c1->SaveAs("HistosEdmundTest.jpg");
	
}