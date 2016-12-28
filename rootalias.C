void RooUnfoldMagic(TString path="")
{
    //Adds the include path and loads the library for RooUnfold
    //Ensure that RooUnfold-1.1.1 and rootalias.C are in the working directory
    //root [] RooUnfoldMagic("MTTunfold.cc") to compile and run MTTunfold.cc
    
    gROOT->ProcessLine(".include RooUnfold-1.1.1/src");
    printf("\nRooUnfold headers added to include path with rootalias\n\n");
    gSystem->Load("RooUnfold-1.1.1/libRooUnfold");
    printf("\nlibRooUnfold loaded with rootalias\n\n");
    
    if (path!="") gROOT->ProcessLine(".x "+path+"++");
}