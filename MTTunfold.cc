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
	TH1F* hGenNch = new TH1F("hGenNch", "Normalised Monte Carlo Generated (Training);nch particles;# Events", 200, 0, 200);
	TH1F* hRecoNch = new TH1F("hRecoNch", "Normalised Monte Carlo Reconstructed (Training);nch particles;# Events", 200, 0, 200);
    TH1F* hUnfoldedNch = new TH1F("hUnfoldedNch", "Unfolded Monte Carlo Reconstructed (Training);nch particles;# Events", 200, 0, 200);
    // TH1F* hNoCh = new TH1F("hNoCh", "No Charged Particles Generated;nch particles;# Events", 300, -0.5, 299.5);

    // Normalisation variables
    int nGenNch, nRecoNch, nUnfoldedNch;
    
	//RooUnfoldResponse variable
	RooUnfoldResponse response_nch(hRecoNch, hGenNch);
    
    
    //Declaration of tree and its branches variables
    TTree* tree = new TTree("EventTree","");
    
    bool isSelected = false;
    bool Charged = false;
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
        tree->SetBranchAddress("isSelected", &isSelected);
        tree->SetBranchAddress("Charged", &Charged);
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
            isSelected = false, Charged = false;
            
            tree->GetEntry(i);
            
            //Filling Response
            {

            }

            int nch_training, gennch_training;
            
            {
                gennch_training = gennch;
                nch_training = nch;
            }
            
            if (isSelected==true && Charged==true)
            {
                response_nch.Fill(nch_training, gennch_training);
                
                hGenNch->Fill(gennch_training);
                nGenNch++;
                
                hRecoNch->Fill(nch_training);
                nRecoNch++;
            }
            else if (isSelected==false && Charged==true)
            {
                response_nch.Miss(gennch_training);
                
                hGenNch->Fill(gennch_training);
                nGenNch++;
            }
            else
            {
                // response_nch.Fake(nch_training);
                // hNoCh->Fill(nch_training);
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
    
    //Histograms for unfolding (data)
	TH1F* hdataRecoNch = new TH1F("hdataNch", "Normalised Data Reconstructed (Unfolding);nch particles;# Events", 200, 0, 200);
    TH1F* hdataUnfoldedNch = new TH1F("hdataUnfoldedNch", "Unfolded Data Reconstructed (Unfolding);nch particles;# Events", 200, 0, 200);
    
    //Normalisation variables
    int ndataRecoNch, ndataUnfoldedNch;
    
    
	//Declaration of tree and its branches variables
	TTree* tree2 = new TTree("EventTree","");

	//Getting filelist of data trees for unfolding
	delete vfiles;
	vfiles = new vector<TString>();
	cout<< "Getting list of files..." << endl;
    vfiles = getListOfFiles("FileListUnfolding.txt");
    cout<< "File list stored." << endl;

	// i_tot = 0;

	//****************************************************UNFOLDING*************************************************************
	cout<< "=======================Testing=====================" <<endl;
    
	//============================================Starting Loop over files====================================================== 
	//==========================================================================================================================
	//(Stops at end of list of files or when reached nevt_max)
	//==========================================================================================================================
   for(vector<TString>::iterator itfiles = vfiles->begin() ; itfiles != vfiles->end(); ++itfiles)
	{
        cout<< "Opening new file." << endl;
		TFile* fileunfold = TFile::Open(*itfiles,"READ");

		//getting the tree from the current file
		cout<< "Getting tree from file." << endl;
		tree2 = (TTree*) fileunfold->Get("newtree");
        
        //adding branches to the tree ----------------------------------------------------------------------
        //!Branch addresses are not updated to match new skims
        tree2->SetBranchStatus("*", 1);
        tree2->SetBranchAddress("isSelected", &isSelected);
        tree2->SetBranchAddress("nch", &nch);
        cout<< "All branches set." << endl;

        //Getting number of events
        int nev = int(tree2->GetEntriesFast());
        cout <<"The current file has " << nev << " entries." << endl;
        
		// //---------------------------------------Starting loop over events------------------------------------------
		// //----------------------------------------------------------------------------------------------------------
		// //(Stops when reached end of file or nevt_max)
		// //----------------------------------------------------------------------------------------------------------
        for(int i = 0; i < (int) nev; ++i , ++i_tot)
        {
            //cout<< "Event Number: " << i_tot << endl;
            //printing the number of events done every hundred events
            //if( ((i_tot+1) % 100) == 0) cout <<int(double(i_tot+1))<<" events done"<<endl;
            //printing the % of events done every 100k evts
            if( ((i_tot+1) % 100000) == 0) cout <<int(double(i_tot+1)/1000)<<"k done"<<endl;

            //Filling the variables defined setting branches
            tree2->GetEntry(i);
            
            //
            {

            }

            int nch_unfolding;
            
            {
                nch_unfolding = nch;
            }
            
            if (isSelected==true)
            {
                hdataRecoNch->Fill(nch_unfolding);
                ndataRecoNch++;
            }
            else if (isSelected==false)
            {

            }
            else
            {

            }
            
        }
		// //---------------------------------------------End of loop over events----------------------------------------

		//Closing current files
		fileunfold->Close();
	}
	//============================================End of loop over files==============================================
	
	cout<< "=======================Unfolding=====================" <<endl;
	cout<< "Unfolding charged particle multiplicity..." << endl;
    
	RooUnfoldBayes nch_unfold(&response_nch, hRecoNch, 3);
	// RooUnfoldBayes data_nch_unfold(&response_nch, hdataRecoNch, 3);
    
    //Temp code to unfold from histogram inistead of tree
    TFile *tempfile = new TFile("ZB_Multiplicity_2.root", "READ");
    TH1F * hdatatempRecoNch = (TH1F*)tempfile->Get("Multiplicity");
    tempfile->Close();    
	RooUnfoldBayes data_nch_unfold(&response_nch, hdatatempRecoNch, 3);
    
	hUnfoldedNch = (TH1F*) nch_unfold.Hreco();
	hdataUnfoldedNch = (TH1F*) data_nch_unfold.Hreco();
    
	cout<< "===================Unfold Complete!===================" << endl;
    
	//*************************************************END UNFOLDING**************************************************	
	

	//=========================================Start writing to output file===========================================
	TFile* outputunfold = new TFile("outputUnfold.root","RECREATE");
	// outputunfold->cd();
    
    hGenNch->Scale(1./nGenNch);
    hRecoNch->Scale(1./nRecoNch);
    // hdataRecoNch->Scale(1./ndataRecoNch);
    
    //Temp code to handle histogram input
    hdatatempRecoNch->Scale(1./hdatatempRecoNch->Integral());
    
    hUnfoldedNch->Scale(1./hUnfoldedNch->Integral());
    hdataUnfoldedNch->Scale(1./hdataUnfoldedNch->Integral());

    hGenNch->SetLineColor(3);
	hGenNch->Write();

    hRecoNch->SetLineColor(2);
	hRecoNch->Write();

    hUnfoldedNch->SetLineColor(1);
    hUnfoldedNch->Write();

    // hdataRecoNch->SetLineColor(2);
	// hdataRecoNch->Write();
    
    //Temp code to handle histogram input
    hdatatempRecoNch->SetLineColor(2);
	hdatatempRecoNch->Write();
    
    hdataUnfoldedNch->SetLineColor(1);
    hdataUnfoldedNch->SetNameTitle("Unfolded Multiplicity", "Unfolded Multiplicity (Normalised)");
    hdataUnfoldedNch->Write();
    
    // hNoCh->Write();
    
    outputunfold->Close();

	// gDirectory->mkdir("UnfoldResult");
	// gDirectory->cd("UnfoldResult");

	// //===========================================End writing output files==========================================
}// end of main loop



#ifndef __CINT__
int main () { MTTunfold(); return 0; }  // Main program when run stand-alone
#endif