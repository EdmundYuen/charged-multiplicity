#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include "TList.h"

using namespace ROOT::Math;

#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

#ifdef __MAKECINT__
    #pragma link C++ class vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >++;
#endif

void hist_operations()
{
	TFile *f1 = new TFile("MB1_Multiplicity.root", "READ");
	TH1D * h_multiplicity_1 = (TH1D*)f1->Get("Multiplicity");
	TH1D * h_multiplicity_zb_1 = (TH1D*)f1->Get("Multiplicity_ZB");
	TH1D * h_multiplicity_hm85_1 = (TH1D*)f1->Get("Multiplicity_HM85");

	TFile *f2 = new TFile("MB3_Multiplicity.root", "READ");
	TH1D * h_multiplicity_2 = (TH1D*)f1->Get("Multiplicity");
	TH1D * h_multiplicity_zb_2 = (TH1D*)f1->Get("Multiplicity_ZB");
	TH1D * h_multiplicity_hm85_2 = (TH1D*)f1->Get("Multiplicity_HM85");

	TFile *f3 = new TFile("MB7_Multiplicity.root", "READ");
	TH1D * h_multiplicity_3 = (TH1D*)f1->Get("Multiplicity");
	TH1D * h_multiplicity_zb_3 = (TH1D*)f1->Get("Multiplicity_ZB");
	TH1D * h_multiplicity_hm85_3 = (TH1D*)f1->Get("Multiplicity_HM85");	

	TFile *f4 = new TFile("MB8_Multiplicity.root", "READ");
	TH1D * h_multiplicity_4 = (TH1D*)f1->Get("Multiplicity");
	TH1D * h_multiplicity_zb_4 = (TH1D*)f1->Get("Multiplicity_ZB");
	TH1D * h_multiplicity_hm85_4 = (TH1D*)f1->Get("Multiplicity_HM85");	
	
	
	TFile *outfile = new TFile("histos.root", "RECREATE");
			
	TList *list = new TList;
	list->Add(h_multiplicity_1);
	list->Add(h_multiplicity_2);
	list->Add(h_multiplicity_3);
	list->Add(h_multiplicity_4);

	TH1D *h_multiplicity = (TH1D*)h_multiplicity_1->Clone();
	h_multiplicity->Reset();
	h_multiplicity->Merge(list);
	
	
	TList *list_zb = new TList;
	list_zb->Add(h_multiplicity_zb_1);
	list_zb->Add(h_multiplicity_zb_2);
	list_zb->Add(h_multiplicity_zb_3);
	list_zb->Add(h_multiplicity_zb_4);

	TH1D *h_multiplicity_zb = (TH1D*)h_multiplicity_zb_1->Clone();
	h_multiplicity_zb->Reset();
	h_multiplicity_zb->Merge(list_zb);
	
	
	TList *list_hm85 = new TList;
	list_hm85->Add(h_multiplicity_hm85_1);
	list_hm85->Add(h_multiplicity_hm85_2);
	list_hm85->Add(h_multiplicity_hm85_3);
	list_hm85->Add(h_multiplicity_hm85_4);

	TH1D *h_multiplicity_hm85 = (TH1D*)h_multiplicity_hm85_1->Clone();
	h_multiplicity_hm85->Reset();
	h_multiplicity_hm85->Merge(list_hm85);
	
	TH1D *h_efficiency_zb = (TH1D*)h_multiplicity_zb_1->Clone();
	h_efficiency_zb->Reset();
	h_efficiency_zb->Divide(h_multiplicity_zb, h_multiplicity);
	h_efficiency_zb->SetTitle("h_efficiency_zb");
	h_efficiency_zb->Write();

	TH1D *h_efficiency_hm85 = (TH1D*)h_multiplicity_hm85_1->Clone();
	h_efficiency_hm85->Reset();
	h_efficiency_hm85->Divide(h_multiplicity_hm85, h_multiplicity);
	h_efficiency_hm85->SetTitle("h_efficiency_hm85");
	h_efficiency_hm85->Write();
	
	outfile->Write();
	
	// f1->Close();
	// f2->Close();
	// f3->Close();
	// f4->Close();
	// outfile->Close();
}