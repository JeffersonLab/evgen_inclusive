#include "eicOutput.h"

eicOutput::eicOutput(eicInput *inp, char *fileout){
    char defaultname[255] = "output.root";

    if( fileout ){
	strcpy(fOutName, fileout);
    } else {
	strcpy(fOutName, defaultname);
    }
    printf("output file set as %s\n", fOutName);

    fTree = new TTree("T", "EIC Monte Carlo");
    fRunTime = inp->GetRunTime();
    fNevt    = ((double) inp->GetNevt());

    InitTree();

    return;
}

eicOutput::~eicOutput(){
    delete fTree;

    return;
}

void eicOutput::InitTree(){
    fTree->Branch("weight", &fData.weight, "weight/D");
    fTree->Branch("rate", &fData.rate, "rate/D");
    fTree->Branch("theta", &fData.theta, "theta/D");
    fTree->Branch("phi", &fData.phi, "phi/D");
    fTree->Branch("x", &fData.x, "x/D");
    fTree->Branch("Ef", &fData.ef, "Ef/D");
    fTree->Branch("crs", &fData.crs, "crs/D");
    fTree->Branch("F1", &fData.f1, "F1/D");
    fTree->Branch("F2", &fData.f2, "F2/D");

    fTree->Branch("F1gZ", &fData.f1gz, "F1gZ/D");
    fTree->Branch("F3gZ", &fData.f3gz, "F3gZ/D");

    fTree->Branch("eta_gZ", &fData.eta_gZ, "eta_gZ/D");

    fTree->Branch("Q2", &fData.Q2, "Q2/D");
    fTree->Branch("pdf", fData.pdf, "pdf[6]/D");
    fTree->Branch("Dpdf", fData.Dpdf, "pdf[5]/D");


    fTree->Branch("Abeam", &fData.Abeam, "Abeam/D");

    fTree->Branch("A_L", &fData.A_L, "A_L/D");
    fTree->Branch("A_T", &fData.A_T, "A_T/D");

    fTree->Branch("Y_L", &fData.Y_L, "Y_L/D");
    fTree->Branch("Y_T", &fData.Y_T, "Y_T/D");

    fTree->Branch("eta_L", &fData.eta_L, "eta_L/D");
    fTree->Branch("eta_T", &fData.eta_T, "eta_T/D");

    fTree->Branch("g1gZ", &fData.g1gz, "g1gZ/D");
    fTree->Branch("g5gZ", &fData.g5gz, "g5gZ/D");

    fTree->Branch("W", &fData.W, "W/D");
    fTree->Branch("y", &fData.y, "y/D");

    fTree->Branch("Y", &fData.Y, "Y/D");
    

    fTree->Branch("Z_ion", &fData.Z_ion, "Z_ion/I");
    fTree->Branch("N_ion", &fData.N_ion, "N_ion/I");
    
    fTree->Branch("particle_id", &fData.particle_id, "particle_id/I");
    fTree->Branch("charge", &fData.charge, "charge/I");
    fTree->Branch("mass", &fData.mass, "mass/D");
    fTree->Branch("mom", &fData.mom, "mom/D");
    

    return;
}

void eicOutput::Write( eicEvent *ev ){
    fData = ev->GetEventData();

    fData.rate    = fData.weight/fNevt;
    fData.weight *= fRunTime/fNevt;

    fTree->Fill();

    return;
}

void eicOutput::Close(){
    fOutfile = new TFile(fOutName, "RECREATE");
    fTree->Write();
    fOutfile->Close();
    delete fOutfile;

    return;
}

