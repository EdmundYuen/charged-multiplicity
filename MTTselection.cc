#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TDirectory.h>
#include <TH1.h>
#include <algorithm>
#include "math.h"
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

using namespace ROOT::Math;

#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

#ifdef __MAKECINT__
    #pragma link C++ class vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >++;
#endif


vector<TString>* getListOfFiles(TString strfiles){

  vector<TString>* vfiles = new vector<TString>;

  if(strfiles.Contains(".root")){
    TChain chain("tree/tree","");
    chain.Add(strfiles);
    TObjArray* fileElements=chain.GetListOfFiles();
    TIter next(fileElements);
    TChainElement *chEl=0;
    while (( chEl=(TChainElement*)next() )) {
      vfiles->push_back(TString(chEl->GetTitle()));
    }
  }
  else if(strfiles.Contains(".txt")){
    ifstream txtfile;
    txtfile.open(strfiles);
    if(!txtfile) {
      cout<<"Unable to read the txt file where the rootfiles are." << endl ;
      cout << strfiles << " doesn't exist." << endl << "Aborting ...";
      exit(0);
    }
    string filename;
    while(txtfile>>filename && filename!="EOF")
      vfiles->push_back(TString(filename));
    txtfile.close();
  }
  else {
    cout << "Unknown type of input to get files. Must contain either .root or .txt extension." << endl << "Aborting ..." << endl;
    exit(0);
  }
  cout << "[getListOfFiles] Will run on " << vfiles->size() << " files" << endl;
  return vfiles;
}


void MTTselection(bool isMC=true)
{
    //===============Getting filelist of trees for training=====================
    
	vector<TString>* vfiles = new vector<TString>();
	cout << "Getting list of files..." << endl;
    
    if (isMC)
    {
        vfiles = getListOfFiles("FileListSelectionMC.txt");
    }
    else
    {
        vfiles = getListOfFiles("FileListSelection.txt");
    }
    
    cout<< "File list stored." << endl;
    
    // ===========================Event selections==============================
    
    const int vtx_number_cut = 1;
    const float vtxz_cut = 10.;
    const float vtxxy_cut = 2.; //0.2?
    
    // ===========================Track cuts====================================
    
    const double eta_cut = 2.4;
    const double pt_cut = 0.5;
    const int highpurity_cut = 1;
    const float d0_d0Err_cut = 3.;
    const float dz_dzErr_cut = 3.;
    const float ptErr_pt_cut = 0.05;
    const int dof_cut = 4;
    
    // ===========================Variables=====================================
    
    
    //=========================== Histos for pT, eta ===========================

    TH1F *hd0d0err = new TH1F ("d0d0err", "d0 d0err", 400, -20, 20);
    // TH1F *hgenPt = new TH1F ("gen Pt", "Normalized data p_{T}", 200, 0, 10);
    // TH1F *hrecoPt = new TH1F ("reco Pt", "Normalized data p_{T}", 200, 0, 10);
    // TH1F *hgenEta = new TH1F("gen eta", "Normalised data #eta", 50, -2.5, 2.5);
    // TH1F *hrecoEta = new TH1F("reco eta", "Normalised data #eta", 50, -2.5, 2.5);

    //============================Starting Loop over files====================== 
    for(vector<TString>::iterator itfiles = vfiles->begin() ; itfiles != vfiles->end(); ++itfiles)
	{
        
        // ===========================retrieve ROOT file========================

        TFile *oldfile = TFile::Open(*itfiles, "READ");
        TTree *oldtree;
        
        if (isMC)
        {
            oldtree = (TTree*)oldfile->Get("mytree/tree");
        }
        else
        {
            oldtree = (TTree*)oldfile->Get("tree/tree");
        }   
        
        // ===========================define variables for new TTree============

        int gennch, nch;
        bool isSelected, Charged;
        
        // ===========================define new ROOT file======================
        
        TFile *newfile = new TFile("selected_"+*itfiles,"recreate");
        TTree *newtree = new TTree("newtree","Selected");
        newtree->Branch("nch",&nch,"nch/I");
        newtree->Branch("isSelected",&isSelected,"isSelected/O");
        
        if(isMC)
        {
            newtree->Branch("gennch",&gennch,"gennch/I");
            newtree->Branch("Charged",&Charged,"Charged/O");
        }

        // ===========================define variables to read TTree============

        int ndata_lumi = 0;
        oldtree->SetBranchAddress("Lumi", &ndata_lumi);
        
        int ndata_run = 0;
        oldtree->SetBranchAddress("Run", &ndata_run);

        vector<double> *dvecdata_vtxx = 0;
        oldtree->SetBranchAddress("vtxx", &dvecdata_vtxx);

        vector<double> *dvecdata_vtxy = 0;
        oldtree->SetBranchAddress("vtxy", &dvecdata_vtxy);

        vector<double> *dvecdata_vtxz = 0;
        oldtree->SetBranchAddress("vtxz", &dvecdata_vtxz);
        
        vector<double> *dvecdata_vtxzBS = 0;
        oldtree->SetBranchAddress("vtxzBS", &dvecdata_vtxzBS);

        Double_t ddata_BSz = 0;
        oldtree->SetBranchAddress("BSz", &ddata_BSz);

        vector<int> *nvecdata_vtxndof = 0;
        oldtree->SetBranchAddress("vtxndof", &nvecdata_vtxndof);
        
        vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *reco_tracks = 0;
        oldtree->SetBranchAddress("trP4", &reco_tracks);
        
        vector<double> *dvecdata_treta = 0;
        oldtree->SetBranchAddress("trEta", &dvecdata_treta);
        
        vector<double> *dvecdata_trpt = 0;
        oldtree->SetBranchAddress("trPt", &dvecdata_trpt);
        
        vector<double> *dvecdata_pterr = 0;
        oldtree->SetBranchAddress("ptErr", &dvecdata_pterr);

        vector<double> *dvecdata_dxyvtxBS = 0;
		oldtree->SetBranchAddress("dxyvtxBS", &dvecdata_dxyvtxBS);
        
		vector<double> *dvecdata_d0err = 0;
		oldtree->SetBranchAddress("d0err", &dvecdata_d0err);
        
        vector<double> *dvecdata_dzvtxBS = 0;
		oldtree->SetBranchAddress("dzvtxBS", &dvecdata_dzvtxBS);
        
		vector<double> *dvecdata_dzerr = 0;
		oldtree->SetBranchAddress("dzerr", &dvecdata_dzerr);
        
        vector<int> *nvecdata_highpurity = 0;
        oldtree->SetBranchAddress("highPurity", &nvecdata_highpurity);
        
        int ndata_vtx = 0;
        
        vector<int> *dvecdata_partCharge = 0;
        
        vector<int> *nvecdata_partStatus = 0;
        
        vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *gen_tracks = 0;
        
        if (isMC)
        {     
            oldtree->SetBranchAddress("vtx", &ndata_vtx);

            oldtree->SetBranchAddress("partCharge", &dvecdata_partCharge);
            
            oldtree->SetBranchAddress("partStatus", &nvecdata_partStatus);

            oldtree->SetBranchAddress("partP4", &gen_tracks);
        }
        
        Int_t ndata_totalEvt = (Int_t)oldtree->GetEntries();
        cout << "There is a total of " << ndata_totalEvt << " events." << endl;

        // =========================== Start of Evt Loop =======================
        for (Int_t i = 0; i < ndata_totalEvt; ++i)
        {
            
            gennch = 0;
            nch = 0;
            isSelected = false;
            Charged = false;
            
            oldtree->GetEntry(i);
            // cout << "At entry " << i << endl;
            
            //good lumisections for runs
            if (
            // ((ndata_lumi >= 90) && (ndata_run == 251721)) ||
            // (isMC) ||
            ((ndata_lumi >= 1 && ndata_lumi <=103) && ndata_run == 254987) ||
            (((ndata_lumi >= 1 && ndata_lumi <=120) || (ndata_lumi >= 153 && ndata_lumi <=214) ) && ndata_run == 254989) ||
            ((ndata_lumi >= 1 && ndata_lumi <=28) && ndata_run == 254993) ||
            ((ndata_lumi >= 4 && ndata_lumi <=414) && ndata_run == 255019) ||
            (((ndata_lumi >= 1 && ndata_lumi <=34) || (ndata_lumi >= 36 && ndata_lumi <=209) || (ndata_lumi >= 306 && ndata_lumi <=325) || (ndata_lumi >= 327 && ndata_lumi <=343) ) && ndata_run == 255029) ||
            (((ndata_lumi >= 31 && ndata_lumi <=101) || (ndata_lumi >= 125 && ndata_lumi <=231) || (ndata_lumi >= 233 && ndata_lumi <=493) || (ndata_lumi >= 586 && ndata_lumi <=1054) || (ndata_lumi >= 1096 && ndata_lumi <=1199)) && ndata_run == 255031)
            )
            {
                if (isMC)
                {
            // ======================= Start of Vertex Loop ====================
                    int ndata_numberofgtrk = gen_tracks->size();
                                
                    // =================== Start of Trk Loop ===================
                    for (int gt = 0; gt < ndata_numberofgtrk; ++gt)
                    {
                        XYZTVector gen_vec = (*gen_tracks)[gt];

                        if (((*dvecdata_partCharge)[gt] != 0) && ((*nvecdata_partStatus)[gt] == 1))
                        {
                            //Indicate that event generated charged particle
                            Charged = true;
                            
                            if ((gen_vec.Pt() >= pt_cut) && (fabs(gen_vec.Eta()) <= eta_cut))
                            {
                                //Fill generated charged particle multiplicity
                                ++gennch;
                            }
                        } 
                    }
                    // =================== End of Trk Loop =====================
                }
            // ======================= End of Vertex Loop ======================

            
            // ======================= Start of Vertex Loop ====================
                if(!isMC) ndata_vtx = dvecdata_vtxz->size();
                
                for (int vtxnumber = 0; vtxnumber < ndata_vtx; ++vtxnumber)
                {
                    if (ndata_vtx == vtx_number_cut)
                    {
                        
                        // if(((*nvecdata_vtxndof)[vtxnumber] > dof_cut) && (fabs((*dvecdata_vtxzBS)[vtxnumber]-ddata_BSz) <= vtxz_cut))
                        if((*nvecdata_vtxndof)[vtxnumber] > dof_cut)
                        {
                            isSelected = true;
                            int ndata_numberofrtrk = reco_tracks->size();
                            
                // =================== Start of Trk Loop =======================
                            for (int rt = 0; rt < ndata_numberofrtrk; ++rt)
                            {
                                // XYZTVector reco_vec = (*reco_tracks)[rt];

                                if ((*nvecdata_highpurity)[rt] == highpurity_cut)
                                {
                                    if (((*dvecdata_trpt)[rt] >= pt_cut) && (fabs((*dvecdata_treta)[rt]) <= eta_cut))
                                    {
                                        if ((*dvecdata_pterr)[rt]/(*dvecdata_trpt)[rt]<ptErr_pt_cut)
                                        {
                                            if ((*dvecdata_dzvtxBS)[rt]/(*dvecdata_dzerr)[rt]<dz_dzErr_cut)
                                            {
                                                hd0d0err->Fill((*dvecdata_dxyvtxBS)[rt]/(*dvecdata_d0err)[rt]);
                                            
                                                if ((*dvecdata_dxyvtxBS)[rt]/(*dvecdata_d0err)[rt]<d0_d0Err_cut)
                                                {
                                                    ++nch;
                                                }
                                            }
                                        }
                                    }
                                }   
                            }
                // =================== End of Trk Loop =========================

                        }
                    }   
                }
            // ======================= End of Vertex Loop ======================

            }
            
            //Ignore events with no charged particles generated
            newtree->Fill();
        }
        // =========================== End of Evt Loop =========================

        
        newtree->Write();
        hd0d0err->Write();

        newfile->Close();
        oldfile->Close();
    }
    // =========================== End of File Loop ============================
    
}

#ifndef __CINT__
int main () { MTTselection(); return 0; }  // Main program when run stand-alone
#endif