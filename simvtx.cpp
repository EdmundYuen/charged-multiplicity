#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TFile.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <TMath.h>
#include "TSystem.h"

using namespace std;

void simvtx()
{
    TFile *myFile = new TFile("tree1.root","READ");
    TTree *tree = (TTree*)myFile->Get("UETree/data");
    TCanvas *myCanvas = new TCanvas("myCanvas","Plots", 400, 800);
    myCanvas->Divide(2,2);

    TH1F *myvtxxplot = new TH1F("vtxxplot", "myvtxx", 10, -0.2, 0.2);
    TH1F *simvtxyplot = new TH1F("simvtxyplot", "mysimvtxy", 100, 0, 5);
    TH1F *simvtxzplot = new TH1F("simvtxzplot", "mysimvtxz", 200, -10, 10);
    //TH3F *simvtx3D = new TH3F("simvtx3D","mysimvtx3D", 100, 0, 5, 100, 0, 5, 200, -10, 10);

    TTreeReader myReader("UETree/data", myFile);
    TTreeReaderValue <Float_t> Myvtxx (myReader,"vtxx");
    //TTreeReaderValue <Float_t> Mysimvtxy (myReader,"simvtxy");
    //TTreeReaderValue <Float_t> Mysimvtxz (myReader,"simvtxz");

    while(myReader.Next())
    {
        myvtxxplot->Fill(*Myvtxx);
        //simvtxyplot->Fill(*Mysimvtxy);
        //simvtxzplot->Fill(*Mysimvtxz);
        //simvtx3D->Fill(*Mysimvtxx, *Mysimvtxy, *Mysimvtxz);
    }

    myvtxxplot->Draw();
/*
    myCanvas->cd(1);
    simvtx3D->Draw();

    myCanvas->cd(2);
    simvtxxplot->Draw("COLZ");

    myCanvas->cd(3);
    simvtxyplot->Draw("LEGO");

    myCanvas->cd(4);
    simvtxzplot->Draw("SURF1");
*/
}
