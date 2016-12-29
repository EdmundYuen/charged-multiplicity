#if !defined(__CINT__) || defined(__MAKECINT__)

//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TDirectory.h>
#include <algorithm>

//Copied from track_selection_edmund2
#include "math.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
using namespace ROOT::Math;

//STANDARD C++ INCLUDES
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

//Unfold macros
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
//#include "RooUnfoldSvd.h"
//#include "RooUnfoldBinByBin.h"
#endif

#define pi 3.14159265358979323846


double deg2rad(double deg)
{
  return (deg * pi / 180);
}


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


void MTTunfold()
{
    #ifdef __CINT__
       gSystem->Load("RooUnfold-1.1.1/libRooUnfold");
       #pragma link C++ class ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >++;
    #endif

    
	TH1::SetDefaultSumw2(1);        
	TH1::AddDirectory(0);

    
	//Histograms for checking training and unfolding (MC)
	TH1F* hGenNch = new TH1F("hGenNch", "Monte Carlo Generated;nch particles;# Events", 300, -0.5, 299.5);
	TH1F* hRecoNch = new TH1F("hRecoNch", "Monte Carlo Reconstructed;nch particles;# Events", 300, -0.5, 299.5);
    TH1F* hUnfoldedNch = new TH1F("hUnfoldedNch", "Unfolded Monte Carlo Reconstructed;nch particles;# Events", 300, -0.5, 299.5);
    TH1F* hNoCh = new TH1F("hNoCh", "No Charged Particles Generated;nch particles;# Events", 300, -0.5, 299.5);

    
	//RooUnfoldResponse variable
	RooUnfoldResponse response_nch(hRecoNch, hGenNch);
    
    
    //Declaration of tree and its branches variables
    TTree* tree = new TTree("EventTree","");
    
    bool isMiss = true;
    bool noCharged = true;
    int gennch = 0;
    int nch = 0;
   
   
    //Getting filelist of trees for training
	vector<TString>* vfiles = new vector<TString>();
	cout<< "Getting list of files..." << endl;
    vfiles = getListOfFiles("FileListTraining.txt");
    cout<< "File list stored." << endl;
    
    
    int i_tot = 0;
   
	//****************************************************TRAINING**************************************************************
    cout<< "=======================Training=====================" <<endl;
    
	//============================================Starting Loop over files====================================================== 
    for(vector<TString>::iterator itfiles = vfiles->begin() ; itfiles != vfiles->end(); ++itfiles)
	{
        cout<< "Opening new file." << endl;
        TFile* file = TFile::Open(*itfiles,"READ");

        //getting the tree from the current file
        cout<< "Getting tree from file." << endl;
        tree = (TTree*) file->Get("newtree");
        
        //adding branches to the tree ----------------------------------------------------------------------
        //!Branch addresses are not updated to match new skims
        tree->SetBranchStatus("*", 1);
        tree->SetBranchAddress("isMiss", &isMiss);
        tree->SetBranchAddress("noCharged", &noCharged);
        tree->SetBranchAddress("gennch", &gennch);
        tree->SetBranchAddress("nch", &nch);
        cout<< "All branches set." << endl;

        
        //Getting number of events
        int nev = int(tree->GetEntriesFast());
        cout <<"The current file has " << nev << " entries." << endl;

        //---------------------------------------Starting loop over events------------------------------------------
        //----------------------------------------------------------------------------------------------------------
        //(Stops when reached end of file or nevt_max)
        //----------------------------------------------------------------------------------------------------------
        for(int i = 0; i < (int) nev; ++i , ++i_tot)
        {
            //cout<< "Event Number: " << i_tot << endl;
            //printing the number of events done every hundred events
            //if( ((i_tot+1) % 100) == 0) cout <<int(double(i_tot+1))<<" events done"<<endl;
            //printing the % of events done every 100k evts
            if( ((i_tot+1) % 100000) == 0) cout <<int(double(i_tot+1)/1000)<<"k done"<<endl;

            //Filling the variables defined setting branches
            gennch = -1, nch = -1;
            isMiss = true, noCharged = true;
            
            tree->GetEntry(i);
            
            //Filling Response
            {

            }

            int nch_training, gennch_training;
            
            {
                gennch_training = gennch;
                nch_training = nch;
            }
            
            if (isMiss==false && noCharged==false)
            {
                response_nch.Fill(nch_training, gennch_training);
                hGenNch->Fill(gennch_training);
                hRecoNch->Fill(nch_training);
            }
            else if (isMiss==true && noCharged==false)
            {
                response_nch.Miss(gennch_training);

            }
            else
            {
                // response_nch.Fake(nch_training);
                hNoCh->Fill(nch_training);
            }
            
        }
        //---------------------------------------------End of loop over events----------------------------------------

        //Closing current files
        file->Close();
    }
	//============================================End of loop over files==============================================
    
	cout<< "End Training." << endl;
	//*************************************************END TRAINING***************************************************	

	// double Jet20bias = 1.;    //weight factor to test bias of trigger matching

	// TString transtext = "";

	// TProfile* profile_n_mult_trans_unfolded = new TProfile("profile_nch_trans"+transtext+"_unfolded", "profile_nch_trans"+transtext+"_unfolded;<# ch>", 0, 10000);


	// TH1F* sample_nch = new TH1F("sample_nch"+transtext, "sample_nch"+transtext+";nch;# Events", 50, -0.5, 49.5);
	// TH1F* unfolded_nch = new TH1F("unfolded_nch"+transtext, "unfolded_nch"+transtext+";nch;# Events", 50, -0.5, 49.5);


	// //Ratio Plot
	// TH1F* ratio_nch_true_unfold;
    
	// //Declaration of tree and its branches variables
	// TTree* tree2 = new TTree("UEevt","");
	// delete mcbranches;
	// delete databranches;
	// mcbranches = new MCTreeBranches;
	// databranches = new DataTreeBranches;

	// //Getting filelist of trees for unfolding
	// delete vfiles;
	// vfiles = new vector<TString>();
	// cout<< "Getting list of files..." << endl;
   // vfiles = getListOfFiles("../filelistsmalltreesdataleftright_akt.txt");
   // cout<< "File list stored." << endl;

	// i_tot = 0;

	//****************************************************UNFOLDING*************************************************************
	// cout<< "=======================Testing=====================" <<endl;
	//============================================Starting Loop over files====================================================== 
	//==========================================================================================================================
	//(Stops at end of list of files or when reached nevt_max)
	//==========================================================================================================================
   // for(vector<TString>::iterator itfiles = vfiles->begin() ; itfiles != vfiles->end(); ++itfiles)
	// {
      // cout<< "Opening new file." << endl;
		// TFile* fileunfold = TFile::Open(*itfiles,"READ");

		// //getting the tree from the current file
		// cout<< "Getting tree from file." << endl;
		// if (isMC) tree2 = (TTree*) fileunfold->Get("MCTree");
		// else tree2 = (TTree*) fileunfold->Get("DataTree");

		// //adding branches to the tree ----------------------------------------------------------------------
		// tree2->SetBranchStatus("*", 1);
		// if (isMC) tree2->SetBranchAddress("UEevent", &(mcbranches->nch_trans));
		// else tree2->SetBranchAddress("UEevent", &(databranches->nch_trans));
	
		// cout<< "All branches set." << endl;

		// //Getting number of events
		// int nev = int(tree2->GetEntriesFast());
		// cout <<"The current file has " << nev << " entries." << endl;
		       
		// //---------------------------------------Starting loop over events------------------------------------------
		// //----------------------------------------------------------------------------------------------------------
		// //(Stops when reached end of file or nevt_max)
		// //----------------------------------------------------------------------------------------------------------
		// for(int i = (int) 0; i < (int) nev; ++i , ++i_tot)
		// {
			// //cout<< "Event Number: " << i_tot << endl;
			// //printing the number of events done every hundred events
			// //if( ((i_tot+1) % 100) == 0) cout <<int(double(i_tot+1))<<" events done"<<endl;
			// //printing the % of events done every 100k evts
			// if( ((i_tot+1) % 100000) == 0) cout <<int(double(i_tot+1)/1000)<<"k done"<<endl;

			// //Filling the variables defined setting branches
			// tree2->GetEntry(i);
			
			// if (isMC)
			// {
				// double gennchtransall = mcbranches->gen_nch_trans;

				// double nchtransall = mcbranches->nch_trans;

				// double nchtransallcorr = mcbranches->nch_trans_corr;

				// // {
					// // profile_n_mult_trans_true->Fill(genleadjetpt, (gennchtransall)/(2*2.0*deg2rad(120)), pthatweight);
				// // }

				// if (mcbranches->isSelected)
				// {


					// {
						// // profile_n_mult_trans_meas->Fill(leadjetpt, nchtransall/(2*2.0*deg2rad(120)), pthatweight);

						// {
							// sample_nch->Fill(nchtransall);
						// }
					// }		
				// }
			// }
			// else
			// {
				// double nchtransall = databranches->nch_trans, nchtransdiff = fabs(databranches->nch_trans_left-databranches->nch_trans_right);
				// double nchtransmax, nchtransmin;

				// double nchtransallcorr = databranches->nch_trans_corr, nchtransdiffcorr = fabs(databranches->nch_trans_left_corr-databranches->nch_trans_right_corr);
				// double nchtransmaxcorr, nchtransmincorr;
				
				// double leadjetpt = databranches->leadjetpt;

				// if (databranches->isSelected)
				// {

					// if ((leadjetpt <= 25) && databranches->trg_minbias)
					// {
						// {
							// {
								// sample_nch->Fill(leadjetpt, nchtransdiff);
							// }
						// }
					// }	
					// if ((leadjetpt > 25) && (leadjetpt <= 50) && databranches->trg_jet20)
					// {
							// {
								// sample_nch->Fill(leadjetpt, nchtransdiff, Jet20bias*(databranches->prsc_jet20)/(databranches->prsc_minbias));
							// }
						// }
					// }	
					// if ((leadjetpt > 50) && databranches->trg_jet40)
					// {
						// {
							// {
								// sample_nch->Fill(leadjetpt, nchtransdiff, (databranches->prsc_jet40)/(databranches->prsc_minbias));
							// }
						// }
					// }
				// }
			// }
		// }
		// //---------------------------------------------End of loop over events----------------------------------------

		// //Closing current files
		// fileunfold->Close();
	// }
	//============================================End of loop over files==============================================
	
	cout<< "=======================Unfolding=====================" <<endl;
	cout<< "Unfolding charged particle multiplicity..." << endl;
	RooUnfoldBayes nch_unfold(&response_nch, hRecoNch, 3);
	hUnfoldedNch = (TH1F*) nch_unfold.Hreco();
    
	cout<< "===================Unfold Complete!===================" << endl;
    
	//*************************************************END UNFOLDING**************************************************	
	

	//=========================================Start writing to output file===========================================
	TFile* outputunfold = new TFile("outputUnfold.root","RECREATE");
	// outputunfold->cd();

    hGenNch->SetLineColor(3);
	hGenNch->Write();

    hRecoNch->SetLineColor(2);
	hRecoNch->Write();

    hUnfoldedNch->SetLineColor(1);
    hUnfoldedNch->Write();
    
    hNoCh->Write();
    
    outputunfold->Close();

	// gDirectory->mkdir("UnfoldResult");
	// gDirectory->cd("UnfoldResult");

	// //===========================================End writing output files==========================================
}// end of main loop



#ifndef __CINT__
int main () { MTTunfold(); return 0; }  // Main program when run stand-alone
#endif