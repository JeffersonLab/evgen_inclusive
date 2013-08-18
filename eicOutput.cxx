#include "eicOutput.h"

eicOutput::eicOutput(eicInput *inp, char *fileout){
    char defaultname[255] = "output.root";

    if( fileout ){
	strcpy(fOutName, fileout);
    } else {
	strcpy(fOutName, defaultname);
    }
    printf("output file set as %s\n", fOutName);
    fOutfile = new TFile(fOutName, "RECREATE");
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
    fTree->Branch("jetth", &fData.jetth, "jetth/D");
    fTree->Branch("phi", &fData.phi, "phi/D");
    fTree->Branch("x", &fData.x, "x/D");
    fTree->Branch("Ef", &fData.ef, "Ef/D");
    fTree->Branch("crs", &fData.crs, "crs/D");
    fTree->Branch("F1", &fData.f1, "F1/D");
    fTree->Branch("F2", &fData.f2, "F2/D");

    fTree->Branch("F1gZ", &fData.f1gz, "F1gZ/D");
    fTree->Branch("F3gZ", &fData.f3gz, "F3gZ/D");

    fTree->Branch("F1Wp", &fData.f1wp, "F1Wp/D");
    fTree->Branch("F3Wp", &fData.f3wp, "F3Wp/D");
    fTree->Branch("F1Wm", &fData.f1wm, "F1Wm/D");
    fTree->Branch("F3Wm", &fData.f3wm, "F3Wm/D");

    fTree->Branch("eta_gZ", &fData.eta_gZ, "eta_gZ/D");

    fTree->Branch("Q2", &fData.Q2, "Q2/D");
    fTree->Branch("pdf", fData.pdf, "pdf[6]/D");
    fTree->Branch("Dpdf", fData.Dpdf, "pdf[5]/D");


    fTree->Branch("Abeam", &fData.Abeam, "Abeam/D");

    fTree->Branch("A_L", &fData.A_L, "A_L/D");
    fTree->Branch("A_L_g1", &fData.A_L_g1, "A_L_g1/D");
    fTree->Branch("A_L_g5", &fData.A_L_g5, "A_L_g5/D");

    fTree->Branch("A_T", &fData.A_T, "A_T/D");

    fTree->Branch("Y_L", &fData.Y_L, "Y_L/D");
    fTree->Branch("Y_T", &fData.Y_T, "Y_T/D");

    fTree->Branch("eta_L", &fData.eta_L, "eta_L/D");
    fTree->Branch("eta_T", &fData.eta_T, "eta_T/D");

    fTree->Branch("g1gZ", &fData.g1gz, "g1gZ/D");
    fTree->Branch("g5gZ", &fData.g5gz, "g5gZ/D");

    fTree->Branch("g1Wp", &fData.g1wp, "g1Wp/D");
    fTree->Branch("g5Wp", &fData.g5wp, "g5Wp/D");
    fTree->Branch("g1Wm", &fData.g1wm, "g1Wm/D");
    fTree->Branch("g5Wm", &fData.g5wm, "g5Wm/D");

    fTree->Branch("W", &fData.W, "W/D");
    fTree->Branch("y", &fData.y, "y/D");

    fTree->Branch("Y", &fData.Y, "Y/D");

    fTree->Branch("AWp", &fData.AWp, "AWp/D");
    fTree->Branch("AWp_g1", &fData.AWp_g1, "AWp_g1/D");
    fTree->Branch("AWp_g5", &fData.AWp_g5, "AWp_g5/D");

    fTree->Branch("AWm", &fData.AWm, "AWm/D");
    fTree->Branch("AWm_g1", &fData.AWm_g1, "AWm_g1/D");
    fTree->Branch("AWm_g5", &fData.AWm_g5, "AWm_g5/D");

    fTree->Branch("Wpweight", &fData.Wpweight, "Wpweight/D");
    fTree->Branch("Wprate",   &fData.Wprate, "Wprate/D");
    fTree->Branch("Wmweight", &fData.Wmweight, "Wmweight/D");
    fTree->Branch("Wmrate",   &fData.Wmrate, "Wmrate/D");

    fTree->Branch("jetp", &fData.jetp, "jetp/D");
    fTree->Branch("jetE", &fData.jetE, "jetE/D");

    fTree->Branch("Z_ion", &fData.Z_ion, "Z_ion/I");
    fTree->Branch("N_ion", &fData.N_ion, "N_ion/I");
    
    fTree->Branch("particle_id", &fData.particle_id, "particle_id/I");
    fTree->Branch("charge", &fData.charge, "charge/I");
    fTree->Branch("mass", &fData.mass, "mass/D");
    fTree->Branch("pf", &fData.pf, "pf/D");

    fTree->Branch("vx", &fData.vx, "vx/D");
    fTree->Branch("vy", &fData.vy, "vy/D");
    fTree->Branch("vz", &fData.vz, "vz/D");

    fTree->Branch("g1_theta", &fData.g1_theta, "g1_theta/D");
    fTree->Branch("g1_phi", &fData.g1_phi, "g1_phi/D");
    fTree->Branch("g1_p", &fData.g1_p, "g1_p/D");

    fTree->Branch("g2_theta", &fData.g2_theta, "g2_theta/D");
    fTree->Branch("g2_phi", &fData.g2_phi, "g2_phi/D");
    fTree->Branch("g2_p", &fData.g2_p, "g2_p/D");

    return;
}

void eicOutput::Write( eicEvent *ev ){
    fData = ev->GetEventData();

    fData.rate    = fData.weight/fNevt;
    fData.weight *= fRunTime/fNevt;

    fData.Wprate    = fData.Wpweight/fNevt;
    fData.Wpweight *= fRunTime/fNevt;

    fData.Wmrate    = fData.Wmweight/fNevt;
    fData.Wmweight *= fRunTime/fNevt;

    fTree->Fill();

    return;
}

void eicOutput::Close(){
    fTree->Write();
    fOutfile->Close();
    delete fOutfile;

    return;
}

void  eicOutput::MakeFileLUND(){

  double weight;
  double theta;
  double phi;
  double Ef;
  double x;
  double W;
  double y;
  double crs;
  double f1;
  double f2;
  double f1gz;
  double f3gz;
  double g1gz;
  double g5gz;
  double Q2;
  double eta_gZ;
  double rate;
  double pdf[6];
  double Dpdf[5];

  double Abeam;
  double Y;

  double A_L;
  double Y_L;
  double eta_L;
	
  double A_T;
  double Y_T;
  double eta_T;

  int particle_id;
  int charge;
  int Z_ion;
  int N_ion;
  double pf;
  double mass;

  TVector3 p_vertex;
  TLorentzVector pi0_g1;
  TLorentzVector pi0_g2;

  double vx, vy, vz;
  double g1_theta, g1_phi, g1_p;
  double g2_theta, g2_phi, g2_p;

 
  
  TString file(fOutName);
  file.ReplaceAll("root","lund"); 
  TChain input_chain("T");
  input_chain.Add(fOutName);

  input_chain.SetBranchAddress("weight",&weight);
  input_chain.SetBranchAddress("theta",&theta);
  input_chain.SetBranchAddress("phi",&phi);
  input_chain.SetBranchAddress("Ef",&Ef);
  input_chain.SetBranchAddress("x",&x);
  input_chain.SetBranchAddress("W",&W);
  input_chain.SetBranchAddress("y",&y);
  input_chain.SetBranchAddress("crs",&crs);
  input_chain.SetBranchAddress("F1",&f1);
  input_chain.SetBranchAddress("F2",&f2);
  input_chain.SetBranchAddress("F1gZ",&f1gz);
  input_chain.SetBranchAddress("F3gZ",&f3gz);
  input_chain.SetBranchAddress("g1gZ",&g1gz);
  input_chain.SetBranchAddress("g5gZ",&g5gz);
  input_chain.SetBranchAddress("Q2",&Q2);
  input_chain.SetBranchAddress("eta_gZ",&eta_gZ);
  input_chain.SetBranchAddress("rate",&rate);
  input_chain.SetBranchAddress("pdf",&pdf);
  input_chain.SetBranchAddress("Dpdf",&Dpdf);
  
  input_chain.SetBranchAddress("Abeam",&Abeam);
  input_chain.SetBranchAddress("Y",&Y);
  
  input_chain.SetBranchAddress("A_L",&A_L);
  input_chain.SetBranchAddress("Y_L",&Y_L);
  input_chain.SetBranchAddress("eta_L",&eta_L);
  
  input_chain.SetBranchAddress("A_T",&A_T);
  input_chain.SetBranchAddress("Y_T",&Y_T);
  input_chain.SetBranchAddress("eta_T",&eta_T);

  input_chain.SetBranchAddress("particle_id",&particle_id);
  input_chain.SetBranchAddress("charge",&charge);
  input_chain.SetBranchAddress("Z_ion",&Z_ion);
  input_chain.SetBranchAddress("N_ion",&N_ion);
  input_chain.SetBranchAddress("pf",&pf);
  input_chain.SetBranchAddress("mass",&mass);

  input_chain.SetBranchAddress("vx",&vx);
  input_chain.SetBranchAddress("vy",&vy);
  input_chain.SetBranchAddress("vz",&vz);

  input_chain.SetBranchAddress("g1_theta",&g1_theta);
  input_chain.SetBranchAddress("g1_phi",&g1_phi);
  input_chain.SetBranchAddress("g1_p",&g1_p);

  input_chain.SetBranchAddress("g2_theta",&g2_theta);
  input_chain.SetBranchAddress("g2_phi",&g2_phi);
  input_chain.SetBranchAddress("g2_p",&g2_p);

  
  Int_t nentries = (Int_t)input_chain.GetEntries();
  double pmod,px,py,pz,nu;
  double MASS_p = 0.938;
  TVector3 vtemp;

  ofstream OUT (file.Data());
  for (int i=0; i<nentries ; i++) {
    input_chain.GetEntry(i);
    if(i % 1000000 == 0 ){
      printf("Analyzed %09d events of total %09d \n",i,nentries);
    }

    if (x>0)   nu = Q2 / (2 * MASS_p * x) ;
    else nu = 0;
    pmod = pow(Ef,2) - pow(mass,2) ;
    if (pmod > 0) pmod = sqrt(pmod);
    else pmod =0;
    px = pmod * sin(theta) * cos(phi);
    py = pmod * sin(theta)* sin(phi);
    pz = pmod * cos(theta);
    // Vertex needs to be in units of cm
    p_vertex.SetXYZ(vx*100.0, vy*100.0, vz*100.0);
    vtemp.SetMagThetaPhi(g1_p, g1_theta, g1_phi);
    pi0_g1.SetVectM(vtemp, 0.0);
    vtemp.SetMagThetaPhi(g2_p, g2_theta, g2_phi);
    pi0_g2.SetVectM(vtemp, 0.0);


    // NB at now, rather than the mass of the particle, will be stored the weight. This is done because the mass is not used by GEMC. For full variables and mass value, one should use the SOLLUND format

//     if (particle_id != 111) {
//       OUT << "1" << " \t " << (Z_ion + N_ion)  << " \t " << Z_ion  << " \t " << "0"  << " \t " << "0" << " \t "  << x << " \t " << y  << " \t " << W  << " \t " << Q2  << " \t " << nu << endl;
//       OUT << " \t " << "1" << " \t " << charge << " \t " << "1" << " \t " << particle_id << " \t " << "0" << " \t " << "0" << " \t " << px << " \t " << py << " \t " << pz << " \t " << Ef << " \t " << weight << " \t " << p_vertex.X()  << " \t " << p_vertex.Y() << " \t " << p_vertex.Z() << endl;
//     }
//     else {
//       OUT << "3" << " \t " << (Z_ion + N_ion)  << " \t " << Z_ion  << " \t " << "0"  << " \t " << "0" << " \t "  << x << " \t " << y  << " \t " << W  << " \t " << Q2  << " \t " << nu << endl;
//       OUT << " \t " << "1" << " \t " << charge << " \t " << "1" << " \t " << particle_id << " \t " << "0" << " \t " << "0" << " \t " << px << " \t " << py << " \t " << pz << " \t " << Ef << " \t " << weight << " \t " << p_vertex.X()  << " \t " << p_vertex.Y() << " \t " << p_vertex.Z() << endl;
//       OUT << " \t " << "2" << " \t " << "0.0" << " \t " << "1" << " \t " << "22" << " \t " << "0" << " \t " << "0" << " \t " << pi0_g1.Px() << " \t " << pi0_g1.Py() << " \t " << pi0_g1.Pz() << " \t " << pi0_g1.E() << " \t " << "0.0" << " \t " << p_vertex.X()  << " \t " << p_vertex.Y() << " \t " << p_vertex.Z() << endl;
//       OUT << " \t " << "3" << " \t " << "0.0" << " \t " << "1" << " \t " << "22" << " \t " << "0" << " \t " << "0" << " \t " << pi0_g2.Px() << " \t " << pi0_g2.Py() << " \t " << pi0_g2.Pz() << " \t " << pi0_g2.E() << " \t " << "0.0" << " \t " << p_vertex.X()  << " \t " << p_vertex.Y() << " \t " << p_vertex.Z() << endl;
//     }
    
      OUT << "1" << " \t " << (Z_ion + N_ion)  << " \t " << Z_ion  << " \t " << "0"  << " \t " << "0" << " \t "  << x << " \t " << y  << " \t " << W  << " \t " << Q2  << " \t " << nu << endl;
      OUT << " \t " << "1" << " \t " << charge << " \t " << "1" << " \t " << particle_id << " \t " << "0" << " \t " << "0" << " \t " << px << " \t " << py << " \t " << pz << " \t " << Ef << " \t " << weight << " \t " << p_vertex.X()  << " \t " << p_vertex.Y() << " \t " << p_vertex.Z() << endl;
  }

  OUT.close();

}

void  eicOutput::MakeFileSOLLUND(){

  double weight;
  double theta;
  double phi;
  double Ef;
  double x;
  double W;
  double y;
  double crs;
  double f1;
  double f2;
  double f1gz;
  double f3gz;
  double g1gz;
  double g5gz;
  double Q2;
  double eta_gZ;
  double rate;
  double pdf[6];
  double Dpdf[5];

  double Abeam;
  double Y;

  double A_L;
  double Y_L;
  double eta_L;
	
  double A_T;
  double Y_T;
  double eta_T;

  int particle_id;
  int charge;
  int Z_ion;
  int N_ion;
  double mass;
  TVector3 p_vertex;
  TLorentzVector pi0_g1;
  TLorentzVector pi0_g2;

  double vx, vy, vz, pf;
  double g1_theta, g1_phi, g1_p;
  double g2_theta, g2_phi, g2_p;

  TString file(fOutName);
  file.ReplaceAll("root","sollund"); 
    TChain input_chain("T");
  input_chain.Add(fOutName);

  input_chain.SetBranchAddress("weight",&weight);
  input_chain.SetBranchAddress("theta",&theta);
  input_chain.SetBranchAddress("phi",&phi);
  input_chain.SetBranchAddress("Ef",&Ef);
  input_chain.SetBranchAddress("x",&x);
  input_chain.SetBranchAddress("W",&W);
  input_chain.SetBranchAddress("y",&y);
  input_chain.SetBranchAddress("crs",&crs);
  input_chain.SetBranchAddress("F1",&f1);
  input_chain.SetBranchAddress("F2",&f2);
  input_chain.SetBranchAddress("F1gZ",&f1gz);
  input_chain.SetBranchAddress("F3gZ",&f3gz);
  input_chain.SetBranchAddress("g1gZ",&g1gz);
  input_chain.SetBranchAddress("g5gZ",&g5gz);
  input_chain.SetBranchAddress("Q2",&Q2);
  input_chain.SetBranchAddress("eta_gZ",&eta_gZ);
  input_chain.SetBranchAddress("rate",&rate);
  input_chain.SetBranchAddress("pdf",&pdf);
  input_chain.SetBranchAddress("Dpdf",&Dpdf);
  
  input_chain.SetBranchAddress("Abeam",&Abeam);
  input_chain.SetBranchAddress("Y",&Y);
  
  input_chain.SetBranchAddress("A_L",&A_L);
  input_chain.SetBranchAddress("Y_L",&Y_L);
  input_chain.SetBranchAddress("eta_L",&eta_L);
  
  input_chain.SetBranchAddress("A_T",&A_T);
  input_chain.SetBranchAddress("Y_T",&Y_T);
  input_chain.SetBranchAddress("eta_T",&eta_T);

  input_chain.SetBranchAddress("particle_id",&particle_id);
  input_chain.SetBranchAddress("charge",&charge);
  input_chain.SetBranchAddress("Z_ion",&Z_ion);
  input_chain.SetBranchAddress("N_ion",&N_ion);
  input_chain.SetBranchAddress("pf",&pf);
  input_chain.SetBranchAddress("mass",&mass);

  input_chain.SetBranchAddress("vx",&vx);
  input_chain.SetBranchAddress("vy",&vy);
  input_chain.SetBranchAddress("vz",&vz);

  input_chain.SetBranchAddress("g1_theta",&g1_theta);
  input_chain.SetBranchAddress("g1_phi",&g1_phi);
  input_chain.SetBranchAddress("g1_p",&g1_p);

  input_chain.SetBranchAddress("g2_theta",&g2_theta);
  input_chain.SetBranchAddress("g2_phi",&g2_phi);
  input_chain.SetBranchAddress("g2_p",&g2_p);


  Int_t nentries = (Int_t)input_chain.GetEntries();
  double pmod,px,py,pz,nu;
  double MASS_p = 0.938;
  TVector3 vtemp;

  ofstream OUT (file.Data());
  for (int i=0; i<nentries ; i++) {
    input_chain.GetEntry(i);
    if(i % 1000000 == 0 ){
      printf("Translated %09d events of total %09d \n",i,nentries);
    }

    // Vertex needs to be in units of cm
    p_vertex.SetXYZ(vx*100.0, vy*100.0, vz*100.0);


    vtemp.SetMagThetaPhi(g1_p, g1_theta, g1_phi);
    pi0_g1.SetVectM(vtemp, 0.0);
    vtemp.SetMagThetaPhi(g2_p, g2_theta, g2_phi);
    pi0_g2.SetVectM(vtemp, 0.0);

    if (x>0)   nu = Q2 / (2 * MASS_p * x) ;
    else nu = 0;
    pmod = pow(Ef,2) - pow(mass,2) ;
    if (pmod > 0) pmod = sqrt(pmod);
    else pmod =0;
    px = pmod * sin(theta) * cos(phi);
    py = pmod * sin(theta)* sin(phi);
    pz = pmod * cos(theta);


    if (particle_id != 111) {
      OUT << "1" << " \t "  << weight << " \t " << (Z_ion + N_ion)  << " \t " << Z_ion  << " \t " << "0"  << " \t " << "0" << " \t "  << x << " \t " << y  << " \t " << W  << " \t " << Q2  << " \t " << nu << endl;
      OUT << " \t " << "1" << " \t " << charge << " \t " << "1" << " \t " << particle_id << " \t " << "0" << " \t " << "0" << " \t " << px << " \t " << py << " \t " << pz << " \t " << Ef << " \t " << mass << " \t " << p_vertex.X()  << " \t " << p_vertex.Y() << " \t " << p_vertex.Z() << endl;
    }
    else {
      OUT << "3" << " \t "  << weight << " \t " << (Z_ion + N_ion)  << " \t " << Z_ion  << " \t " << "0"  << " \t " << "0" << " \t "  << x << " \t " << y  << " \t " << W  << " \t " << Q2  << " \t " << nu << endl;
      OUT << " \t " << "1" << " \t " << charge << " \t " << "1" << " \t " << particle_id << " \t " << "0" << " \t " << "0" << " \t " << px << " \t " << py << " \t " << pz << " \t " << Ef << " \t " << mass << " \t " << p_vertex.X()  << " \t " << p_vertex.Y() << " \t " << p_vertex.Z() << endl;
      OUT << " \t " << "2" << " \t " << "0.0" << " \t " << "1" << " \t " << "22" << " \t " << "0" << " \t " << "0" << " \t " << pi0_g1.Px() << " \t " << pi0_g1.Py() << " \t " << pi0_g1.Pz() << " \t " << pi0_g1.E() << " \t " << "0.0" << " \t " << p_vertex.X()  << " \t " << p_vertex.Y() << " \t " << p_vertex.Z() << endl;
      OUT << " \t " << "3" << " \t " << "0.0" << " \t " << "1" << " \t " << "22" << " \t " << "0" << " \t " << "0" << " \t " << pi0_g2.Px() << " \t " << pi0_g2.Py() << " \t " << pi0_g2.Pz() << " \t " << pi0_g2.E() << " \t " << "0.0" << " \t " << p_vertex.X()  << " \t " << p_vertex.Y() << " \t " << p_vertex.Z() << endl;
    }    
  }

  OUT.close();


}
