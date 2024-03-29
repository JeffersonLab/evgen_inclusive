#include "eicPhysics.h"

#include "cteqpdf.h"

#include "wiser_pion.h"
//#include " "

#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "TRandom.h"

//#include "../G2PPhys/G2PPhysBase.h"
//#include "../G2PPhys/G2PPhysPB.h"


extern"C" {
void wiser_all_sig_(float *E0,float *P,float *THETA_DEG,float *RAD_LEN,int *TYPE,float *SIGMA);
}

eicPhysics::eicPhysics(){
    int mode = 400;
    fPDF = cteq_pdf_alloc_id(mode);
    
    fRandom = new TRandom2(0);
    printf("Seed number %d\n",fRandom->GetSeed());
    ReadPolTable();

    fHaveTotalXs = false;
    fTotalXs_pip     = 0.0;
    fTotalXs_pim     = 0.0;
    fTotalXs_Kp     = 0.0;
    fTotalXs_Km     = 0.0;
    fTotalXs_p     = 0.0;
    fTotalXs_pbar     = 0.0;

    return;
}

eicPhysics::~eicPhysics(){
    cteq_pdf_free(fPDF);
    return;
}

//eDIS generator
void eicPhysics::MakeEvent(eicBeam *beam, eicIon *ion, eicEvent *ev , eicModel *model){
  // target info
    double tglx = model->GetLx();
    double tgly = model->GetLy();
    double tglength = model->GetLength();

    TVector3 tgtoff = model->GetTgtOffset();

    // We're going to work in the lab frame
    
    // Get beta of ion
    double ionp  = sqrt( pow(ion->GetEnergy(),2.0) - pow(ion->GetMass(),2.0) );
    double beta  = ionp/ion->GetEnergy();
    double gamma = 1.0/sqrt(1.0-beta*beta);

    double e_lab = beam->GetEnergy()*gamma*(1.0+beta);

    TLorentzVector e_rest( 0.0, 0.0, e_lab, e_lab);
    TVector3 blab(0.0, 0.0, beta);

    // Sample evenly in theta (in EIC frame) and phi (lab frame)

    //double th_eic   = GetTheta();
    //double th       = acos((cos(th_eic) + beta)/(beta*cos(th_eic)+1.0));

    // Calculate J = dcos/dcos' so we can fix our cross section
    // which is dx/dcos
    //double J        = (1.0-beta*beta)/pow( 1.0 + beta*cos(th_eic),2.0);
    double J        = 1.0;

    double ph       = GetPhi();
    //double ef       = GetEf(e_rest.Energy(), th);

    //double V        = GetEfmax(e_rest.Energy(), th)*4.0*3.14159;

    double x     = fRandom->Uniform();
    double y     = fRandom->Uniform(e_lab/(MASS_P*x/2.0+e_lab));
    double V     = 2.0*3.14159*e_lab/(MASS_P*x/2.0+e_lab);

    double numax  = e_lab;
    double Q2     = 2.0*x*y*(MASS_P*e_lab);
    double ef     = e_lab*(1.0-y);

    double th     = acos(1.0-MASS_P*x*y/(numax*(1.0-y)));

    TLorentzVector ki = TLorentzVector(0.0, 0.0, e_rest.Energy(), e_rest.Energy());

    TLorentzVector kf = TLorentzVector(ef*sin(th), 0.0, ef*cos(th), ef);
    kf.SetPhi(ph);

    double M     = MASS_P;

    TLorentzVector pi;
    pi.SetXYZM( 0.0, 0.0, 0.0, M);

    TLorentzVector q = ki-kf;

//    double Q2    = -q*q;


    double W     = (q+pi).Mag();
//    double s     = (ki+pi).Mag2();

    nucl n;
    double A = ((double) (ion->GetZ()+ion->GetN()));
    double prot_prob = ((double) ion->GetZ())/A;
    // Determine which type of nucleon we hit
    if( fRandom->Uniform() < prot_prob ){
	n = kProton;
    } else {
	n = kNeutron;
    }

    double f1 = A*F1(x,Q2,n);
    double f2 = A*F2(x,Q2,n);

    double f1gz = A*F1gZ(x,Q2,n);
    double f3gz = A*F3gZ(x,Q2,n);

    double f1z = A*F1Z(x,Q2,n);
    double f3z = A*F3Z(x,Q2,n);
    ////////////////////////////////////////////////////
    //  Parity violation
    double G_fermi = 1.16637e-5; // GeV^-2
    double MZ = 91.1876; // GeV
    double MW = 80.399; // GeV
    double sin2thW = 0.23119;

    double Ynum = 2.0*x*(1.0-y/2.0);
    double Yden = 2.0*x*y + 4.0*x*(1.0-y-x*y*M/(2.0*numax))/y;

    double Y = Ynum/Yden;

    double eta_gZ = G_fermi*MZ*MZ*Q2*137.0/(2.0*sqrt(2.0)*3.14159)/(Q2+MZ*MZ);
    double eta_Z  = eta_gZ*eta_gZ;

    double gA = -1.0/2.0;
    double gV = -1.0/2.0 + 2.0*sin2thW;

    double f2nc = f2 - gV*eta_gZ*f1gz*2.0*x+(gV*gV+gA*gA)*eta_Z*f1z*2.0*x;
    double f1nc = f2nc/2.0/x;
    double f3nc = -gA*eta_gZ*f3gz+2.0*gV*gA*eta_Z*f3z;

    // From Halzen and Marten (Agrees with PDG)
//    double ds_dxdy 
//	= 2.0*3.14159*(x*y*y*f1 + ((1.0-y)-M*x*y/(2.0*numax))*f2)/(x*x*y*y*M*numax*137.0*137.0);

    //  From PDG
    double ds_dxdy 
//	= 4.0*3.14159*((1.0-y-pow(x*y*MASS_P,2.0)/Q2)*f2+y*y*x*f1)
//	  /(x*y*Q2*137.0*137.0);
	= 4.0*3.14159*((1.0-y-pow(x*y*MASS_P,2.0)/Q2)*f2nc+y*y*x*f1nc+(y-y*y/2.0)*x*f3nc)
	  /(x*y*Q2*137.0*137.0);
    
    
    
//***********************Added By Ziheng**************************
    double ds_dOmegadE = ds_dxdy*ef*J/(2.0*3.14159*M*e_lab*y);    // in GeV-2 / GeV
    ds_dOmegadE = ds_dOmegadE* (0.197*0.197*1e7); // in nb/GeV-sr
    // convert to nb/GeV-sr,
    //GeV-2 to nb:  1GeV-2=(0.197*0.197*1e-30)m2 = (0.197*0.197*1e7) nb

//    double ds_dOmegadE = ds_dxdy*ef*J/(2.0*M*e_rest.Energy()*3.14159*y);
//*****************************************************************

    
    

//    double Abeam =-eta_gZ*(gA*f1gz/f1 + Y*gV*f3gz/f1);
    double Abeam = (eta_gZ*gA*f1gz   - 2.0*eta_Z*gV*gA*f1z  
	          + eta_gZ*Y*gV*f3gz - Y*(gV*gV+gA*gA)*eta_Z*f3z)/(f1nc+Y*f3nc);

    double g1gz = A*g1gZ(x,Q2,n);
    double g5gz = A*g5gZ(x,Q2,n);

//    double g1z = A*g1Z(x,Q2,n);
//    double g5z = A*g5Z(x,Q2,n);

    /* Simplify.  eta_L here is high by a factor of 2...
    double eta_L = eta_gZ*2.0*(x*y*y + x*(1.0-y)*(2.0-x*y*M/numax))
	           /(x*y*y + 2.0*x*(1.0-y-x*y*M/(2.0*numax)));
    double Y_L   = y*(2.0-y-x*y*M/numax)/(y*y + (1.0-y)*(2.0-x*y*M/numax));
    */

    double eta_L = eta_gZ;
    double Y_L   = y*(2.0-y)/(y*y + (1.0-y)*2.0);

    double eta_T = 0.0;
    double Y_T   = 0.0;

    double A_L = eta_L*(gV*g5gz/f1 + Y_L*gA*g1gz/f1);
    double A_T = 0.0;

    double A_L_g1 = eta_L*(Y_L*gA*g1gz/f1);
    double A_L_g5 = eta_L*(gV*g5gz/f1);
    
    if (std::isnan(Abeam) ||  std::isinf(Abeam)) {cout << "Abeam " << Abeam << endl; Abeam=0;}
    if (std::isnan(A_L) ||  std::isinf(A_L)) {cout << "A_L " << A_L << endl; A_L=0;}
    if (std::isnan(A_L_g1) ||  std::isinf(A_L_g1)) {cout << "A_L_g1 " << A_L_g1 << endl; A_L_g1=0;}
    if (std::isnan(A_L_g5) ||  std::isinf(A_L_g5)) {cout << "A_L_g5 " << A_L_g5 << endl; A_L_g5=0;}    

    double Yp = 1.0 + pow(1.0-y,2.0);
    double Ym = 1.0 - pow(1.0-y,2.0);
    double etaW = 0.5*pow(G_fermi*MW*MW*137.0*Q2/(Q2+MW*MW)/(4.0*3.14159),2.0);

    double f2wp = A*F2W(x,Q2,n,1);
    double f3wp = A*F3W(x,Q2,n,1);
    double g1wp = A*g1W(x,Q2,n,1);
    double g4wp = A*g4W(x,Q2,n,1);

    double f2wm = A*F2W(x,Q2,n,-1);
    double f3wm = A*F3W(x,Q2,n,-1);
    double g1wm = A*g1W(x,Q2,n,-1);
    double g4wm = A*g4W(x,Q2,n,-1);

    double f1wp = A*F1W(x,Q2,n, 1);
    double f1wm = A*F1W(x,Q2,n,-1);
//    double g5wp = A*g5W(x,Q2,n, 1);
//    double g5wm = A*g5W(x,Q2,n,-1);

    // For now we are assuming that the electrons/positrons
    // completely polarized (as per Vogelsang)
    double dsWmdxy = (etaW*8.0*3.14159/(137.0*137.0*x*y*Q2))*
	(Yp*f2wm + Ym*x*f3wm);
    double dsWpdxy = (etaW*8.0*3.14159/(137.0*137.0*x*y*Q2))*
	(Yp*f2wp - Ym*x*f3wp);
	
    // W- exchange needs extra minus sign
    
    // Sigma+ and Sigma- each get half the luminosity,
    // So dsWmdxy needs a factor of 2 reduction to get 
    // the asymmetry 
    double AWm = (etaW*16.0*3.14159/(137.0*137.0*x*y*Q2))*
	(-Yp*g4wm + 2.0*Ym*x*g1wm)/(2.0*dsWmdxy);

    double AWp = (etaW*16.0*3.14159/(137.0*137.0*x*y*Q2))*
	(-Yp*g4wp - 2.0*Ym*x*g1wp)/(2.0*dsWpdxy);

    double AWm_g1 = (etaW*16.0*3.14159/(137.0*137.0*x*y*Q2))*
	(2.0*Ym*x*g1wm)/(2.0*dsWmdxy);
    double AWm_g5 = (etaW*16.0*3.14159/(137.0*137.0*x*y*Q2))*
	(-Yp*g4wm)/(2.0*dsWmdxy);

    double AWp_g1 = (etaW*16.0*3.14159/(137.0*137.0*x*y*Q2))*
	(-2.0*Ym*x*g1wp)/(2.0*dsWpdxy);
    double AWp_g5 = (etaW*16.0*3.14159/(137.0*137.0*x*y*Q2))*
	(-Yp*g4wp)/(2.0*dsWpdxy);
	
      if(ds_dxdy==0) cout << "ds_dxdy " << ds_dxdy << endl;
      if(dsWpdxy==0) {cout << "dsWpdxy " << dsWpdxy << endl; AWp=0;AWp_g1=0;AWp_g5=0;}
      if(dsWmdxy==0) {cout << "dsWmdxy " << dsWmdxy << endl; AWm=0;AWm_g1=0;AWm_g5=0;}

    ////////////////////////////////////////////////////

    // Boost q, kf back to EIC frame
    kf.Boost(-blab);
    TLorentzVector qlab = q;
    qlab.Boost(-blab);

//    double th_eic = kf.Theta();

    // Determine the vertex
    // Generating the vertex randomly in the target
    TVector3 vert;

    double vert_x, vert_y,vert_z;

    vert_x = fRandom->Uniform((-tglx/2),(tglx/2)) + tgtoff.X();
    vert_y = fRandom->Uniform((-tgly/2),(tgly/2)) + tgtoff.Y();
    vert_z = fRandom->Uniform((-tglength/2),(tglength/2))+ tgtoff.Z();;

    vert.SetXYZ(vert_x,vert_y,vert_z);
   
    eventdata data;

    /*
    if( 0.0 < ds_dOmegadE && ds_dOmegadE < 1e9 ){
	data.weight  = ds_dOmegadE*(0.197*0.197*1e-30)*V; // GeV^-2 to m^2
	data.weight *= beam->GetLumin();
    } else {
	// Unphysical for some reason
	data.weight = 0.0;
    }
    */

    double ds_dxdydphi = ds_dxdy/(2.0*3.14159);
    double dsWp_dxdydphi = dsWpdxy/(2.0*3.14159);
    double dsWm_dxdydphi = dsWmdxy/(2.0*3.14159);

    if( 0.0 < ds_dxdy && ds_dxdy < 1e9 ){
	data.weight  = ds_dxdydphi*(0.197*0.197*1e-30)*V*J;
    // GeV^-2 to m^2: 1GeV-1=0.197*1e-15 m, 1GeV-2=(0.197*0.197*1e-30)m2 = (0.197*0.197*1e7) nb
	data.weight *= beam->GetLumin();
    } else {
	// Unphysical for some reason
	data.weight = 0.0;
    }

    if( 0.0 < dsWp_dxdydphi && dsWp_dxdydphi < 1e9 ){
	data.Wpweight  = dsWp_dxdydphi*(0.197*0.197*1e-30)*V*J; // GeV^-2 to m^2
	data.Wpweight *= beam->GetLumin();
    } else {
	// Unphysical for some reason
	data.Wpweight = 0.0;
    }

    if( 0.0 < dsWm_dxdydphi && dsWm_dxdydphi < 1e9 ){
	data.Wmweight  = dsWm_dxdydphi*(0.197*0.197*1e-30)*V*J; // GeV^-2 to m^2
	data.Wmweight *= beam->GetLumin();
    } else {
	// Unphysical for some reason
	data.Wmweight = 0.0;
    }


    data.rate   = 0.0;
    data.Wprate   = 0.0;
    data.Wmrate   = 0.0;

    data.crs    = ds_dOmegadE;
//    data.crs    = ds_dxdy;

    data.ef     = kf.Energy();
    data.theta  = kf.Theta();
    data.jetth  = qlab.Theta();

    data.jetE   = 0.0;
    data.jetp   = qlab.Vect().Mag();

    data.phi    = kf.Phi();
    data.x      = x;
    data.y      = y;

    data.f1     = f1;
    data.f2     = f2;

    data.g1gz   = g1gz;
    data.g5gz   = g5gz;

    data.eta_gZ = eta_gZ;

    data.f1gz     = f1gz;
    data.f3gz     = f3gz;

    data.Q2     = Q2;
    data.W      = W;

    data.Abeam  = Abeam;
    data.Y      = Y;

    data.A_L     = A_L;
    data.eta_L   = eta_L;
    data.Y_L     = Y_L;

    data.A_L_g1     = A_L_g1;
    data.A_L_g5     = A_L_g5;

    data.A_T     = A_T;
    data.eta_T   = eta_T;
    data.Y_T     = Y_T;

    data.AWp     = AWp;
    data.AWm     = AWm;

    data.AWp_g1  = AWp_g1;
    data.AWp_g5  = AWp_g5;
    data.AWm_g1  = AWm_g1;
    data.AWm_g5  = AWm_g5;

    data.g1wm    = g1wm;
    data.g1wp    = g1wp;
    data.g5wm    = g1wm;
    data.g5wp    = g1wp;

    data.f1wm    = f1wm;
    data.f1wp    = f1wp;
    data.f3wm    = f3wm;
    data.f3wp    = f3wp;

    data.mass    = 0.000511;
    data.particle_id = 11;
    data.charge = -1;
//     data.pf =  kf.Energy();
    data.pf =  kf.P();    
    data.Z_ion = ion->GetZ();
    data.N_ion = ion->GetN();

    data.vx = vert.X();
    data.vy = vert.Y();
    data.vz = vert.Z();
 

    // Look at pdf data
    double useQ2;
    if( Q2 < 1.0 ){ useQ2 = 1.0;} else {useQ2 = Q2;}
    data.pdf[0] = cteq_pdf_evolvepdf(fPDF, 1, x, sqrt(useQ2));
    data.pdf[1] = cteq_pdf_evolvepdf(fPDF, 2, x, sqrt(useQ2));
    data.pdf[2] = cteq_pdf_evolvepdf(fPDF, 3, x, sqrt(useQ2));
    data.pdf[3] = cteq_pdf_evolvepdf(fPDF,-1, x, sqrt(useQ2));
    data.pdf[4] = cteq_pdf_evolvepdf(fPDF,-2, x, sqrt(useQ2));
    data.pdf[5] = cteq_pdf_evolvepdf(fPDF,-3, x, sqrt(useQ2));

    for( int ipdf = 0; ipdf < 5; ipdf++ ){
	data.Dpdf[ipdf] = getDeltaq(x, useQ2, ipdf);
    }
 
    ev->SetEventData(data);

    return;
}

//hadron generator according to wiser (uniformly distributed with weight)
void eicPhysics::MakeEvent2(eicBeam *beam, eicIon *ion, eicEvent *ev , eicModel *model) {

  double radlen=0., mass=0., weight_v;
  int particle_id, charge;
  int modelsig = model->GetModel();
  radlen = model->GetRadLen();
  double tglx = model->GetLx();
  double tgly = model->GetLy();
  double tglength = model->GetLength();

  TVector3 tgtoff = model->GetTgtOffset();


  // Initialize variables
  particle_id = 1e9;
  charge      = 1e9;  
  
  TLorentzVector ion_lab(ion->GetMom()*sin(ion->GetTheta())*cos(ion->GetPhi()),ion->GetMom()*sin(ion->GetTheta())*sin(ion->GetPhi()),ion->GetMom()*cos(ion->GetTheta()),sqrt(pow(ion->GetMom(),2)+pow(ion->GetMass(),2)));

  TLorentzVector beam_lab(beam->GetMom()*sin(beam->GetTheta())*cos(beam->GetPhi()),beam->GetMom()*sin(beam->GetTheta())*sin(beam->GetPhi()),beam->GetMom()*cos(beam->GetTheta()),sqrt(pow(beam->GetMom(),2)+pow(beam->GetMass(),2)));
     
  TVector3 bv_ion_lab=ion_lab.BoostVector();
  TLorentzVector beam_ionrest=beam_lab;
  beam_ionrest.Boost(-bv_ion_lab);
//   cout << "!!! " << beam_lab.P() << "\t" <<beam_lab.Theta() << "\t" << beam_lab.Phi() << "\t" << beam_lab.E() << endl; 
//   cout << beam_ionrest.P() << "\t" <<beam_ionrest.Theta() << "\t" << beam_ionrest.Phi() << "\t" << beam_ionrest.E() << endl;
  
  double En_beam = beam_ionrest.E(); 
  
//   nucl n;
  double A = ((double) (ion->GetZ()+ion->GetN()));
  double prot_prob = ((double) ion->GetZ())/A;
  // Determine which type of nucleon we hit
//   if( fRandom->Uniform() < prot_prob ){
//     n = kProton;
//   } else {
//     n = kNeutron;
//   }
  
  //  TF2 *func;  
  if (modelsig == 0) { // pi+
    particle_id = 211; 
    charge = +1;
    mass = 0.1396; // mass in GeV
    //    func = new TF2("sigma_pip",Wiser_func_pip,0, En_beam,0,360,2);
  } 
  else if (modelsig == 1) { //pi-
    particle_id = -211;
    charge = -1;
    mass = 0.1396; // mass in GeV
    //   func = new TF2("sigma_pip",Wiser_func_pip,0, En_beam,0,360,2);
  }
  else if (modelsig == 2) { //pi0
    particle_id = 111;
    charge = 0;
    mass = 0.1350; // mass in GeV
    //   func = new TF2("sigma_pip",Wiser_func_pi0,0, En_beam,0,360,2);
  }  
  else if (modelsig == 3) { //K+
    particle_id = 321;
    charge = +1;
    mass = 0.4937; // mass in GeV
    //   func = new TF2("sigma_pip",Wiser_func_pi0,0, En_beam,0,360,2);
  }  
  else if (modelsig == 4) { //K-
    particle_id = -321;
    charge = -1;
    mass = 0.4937; // mass in GeV
    //   func = new TF2("sigma_pip",Wiser_func_pi0,0, En_beam,0,360,2);
  }  
  else if (modelsig == 5) { //Ks
    particle_id = 310;
    charge = 0;
    mass = 0.4976; // mass in GeV
    //   func = new TF2("sigma_pip",Wiser_func_pi0,0, En_beam,0,360,2);
  }  
  else if (modelsig == 6) { //Kl
    particle_id = 130;
    charge = 0;
    mass = 0.4976; // mass in GeV
    //   func = new TF2("sigma_pip",Wiser_func_pi0,0, En_beam,0,360,2);
  }    
  else if (modelsig == 7) { //p
    particle_id = 2212;
    charge = +1;
    mass = 0.9383; // mass in GeV
    //   func = new TF2("sigma_pip",Wiser_func_pi0,0, En_beam,0,360,2);
  }  
  else if (modelsig == 8) { //p-bar
    particle_id = -2212;
    charge = -1;
    mass = 0.9383; // mass in GeV
    //   func = new TF2("sigma_pip",Wiser_func_pi0,0, En_beam,0,360,2);
  }
        
 
  //  func->SetParameters(En_beam,radlen);
  //  double mom_pi= 0, theta_pi= 0;
  //  func->GetRandom2(mom_pi,theta_pi);
  double x_pi=0.,y_pi=0.,z_pi=0.; 
  double mom_pi = fRandom->Uniform(En_beam); 
  fRandom->Sphere(x_pi,y_pi,z_pi,mom_pi); // generate random  vectors, uniformly distributed over the surface
  // of a sphere of radius mom_pi
  TVector3 v3_pi(x_pi,y_pi,z_pi);
  double theta_pi = v3_pi.Theta();
  double phi_pi = v3_pi.Phi();
  
  double ef = sqrt(pow(mom_pi,2) + pow(mass,2));
  
  TLorentzVector pi_ionrest(0,0,0,0),pi_lab(0,0,0,0);    
  pi_ionrest.SetPxPyPzE(mom_pi*sin(theta_pi)*cos(phi_pi),mom_pi*sin(theta_pi)*sin(phi_pi),mom_pi*cos(theta_pi),ef);
//   cout << pi_ionrest.P() << "\t" <<pi_ionrest.Theta() << "\t" << pi_ionrest.Phi() << "\t" << pi_ionrest.E() << endl;   
  pi_ionrest.Rotate(TVector3(0,0,1).Angle(beam_ionrest.Vect()),TVector3(0,0,1).Cross(beam_ionrest.Vect()));   //wiser only takes theta relative to +z as assumed electron beam direction,need to rotate to actual electron beam direction
//   cout << pi_ionrest.P() << "\t" <<pi_ionrest.Theta() << "\t" << pi_ionrest.Phi() << "\t" << pi_ionrest.E() << endl;   
  pi_lab=pi_ionrest;  
  pi_lab.Boost(bv_ion_lab);
//   cout << pi_ionrest.P() << "\t" <<pi_ionrest.Theta() << "\t" << pi_ionrest.Phi() << "\t" << pi_ionrest.E() << endl; 
//   cout << pi_lab.P() << "\t" <<pi_lab.Theta() << "\t" << pi_lab.Phi() << "\t" << pi_lab.E() << endl;
  
  // Generating the vertex randomly in the target
  TVector3 vert;
  double vert_x, vert_y,vert_z;

  vert_x = fRandom->Uniform((-tglx/2),(tglx/2))+tgtoff.X(); ;
  vert_y = fRandom->Uniform((-tgly/2),(tgly/2))+tgtoff.Y(); ;
  double targprop = fRandom->Uniform();
//  vert_z = fRandom->Uniform((-tglength/2),(tglength/2))+tgtoff.Z(); ;
  vert_z = (targprop-0.5)*tglength+tgtoff.Z();

  vert.SetXYZ(vert_x,vert_y,vert_z);

  double intrad = 2.0*log(En_beam/0.000511)/(137.0*3.14159);

  radlen = targprop*radlen*100.*(4.0/3.0) + intrad*100.0; 

  //  The radiation spectrum is proportional to 4/3/k (for small k).
  //  Internal factor for lab system when E >> Mn is (alpha/pi)ln(Elab/me)


  float En_beam2 = float(En_beam);
  float mom_pi2 = float(mom_pi);
  float theta_pi2 = float(theta_pi*180.0/3.14159); // theta pi in deg
  float radlen2 = float(radlen);
  float weight_f;
  int type;
  Gamma1.SetXYZT(0.,0.,0.,0.);
  Gamma2.SetXYZT(0.,0.,0.,0.);
  Gamma1_vt.SetXYZ(0.,0.,0.);
 
  if (particle_id == 211) {
    //   weight_v = WISER_ALL_FIT(mom_pi);
    //  weight_v = WISER_ALL_SIG(En_beam,mom_pi,theta_pi,radlen,1);   
    type = 1;
    wiser_all_sig_(&En_beam2,&mom_pi2,&theta_pi2,&radlen2,&type,&weight_f);
    weight_v = prot_prob*double(weight_f);
    type = 2;    
    wiser_all_sig_(&En_beam2,&mom_pi2,&theta_pi2,&radlen2,&type,&weight_f);
    weight_v = (1-prot_prob)*double(weight_f) + weight_v;        
    //   cout << weight_f << " \t " << weight_v << endl;
    //    printf("%f \n",weight_v);
  }
  else if (particle_id == -211) {
    // weight_v = WISER_ALL_SIG(En_beam,mom_pi,theta_pi,radlen,2);
    type = 2;
    wiser_all_sig_(&En_beam2,&mom_pi2,&theta_pi2,&radlen2,&type,&weight_f);
    weight_v = prot_prob*double(weight_f);
    type = 1;    
    wiser_all_sig_(&En_beam2,&mom_pi2,&theta_pi2,&radlen2,&type,&weight_f);
    weight_v = (1-prot_prob)*double(weight_f) + weight_v;        
  }
  else if (particle_id == 111) {
    //  weight_v = 0.5 * ( WISER_ALL_SIG(En_beam,mom_pi,theta_pi,radlen,1) + WISER_ALL_SIG(En_beam,mom_pi,theta_pi,radlen,2)) ;
    type = 1;
    wiser_all_sig_(&En_beam2,&mom_pi2,&theta_pi2,&radlen2,&type,&weight_f);
    weight_v = 2.0 * double(weight_f);
    type = 2;
    wiser_all_sig_(&En_beam2,&mom_pi2,&theta_pi2,&radlen2,&type,&weight_f);
    weight_v = 2.0 * double(weight_f) + weight_v;
//     TVector3 vp(1.,1.,1.);
//     vp.SetMag(mom_pi);
//     vp.SetTheta(theta_pi);
//     vp.SetPhi(phi_pi);
//     Decay_pi0(vp,vert);    
  }
  else if (particle_id == 321) {
    //   weight_v = WISER_ALL_FIT(mom_pi);
    //  weight_v = WISER_ALL_SIG(En_beam,mom_pi,theta_pi,radlen,1);
    type = 3;
    wiser_all_sig_(&En_beam2,&mom_pi2,&theta_pi2,&radlen2,&type,&weight_f);
    weight_v = prot_prob*double(weight_f);
    type = 4;    
    wiser_all_sig_(&En_beam2,&mom_pi2,&theta_pi2,&radlen2,&type,&weight_f);
    weight_v = (1-prot_prob)*double(weight_f) + weight_v;        
    //   cout << weight_f << " \t " << weight_v << endl;
    //    printf("%f \n",weight_v);
  }
  else if (particle_id == -321) {
    // weight_v = WISER_ALL_SIG(En_beam,mom_pi,theta_pi,radlen,2);
    type = 4;
    wiser_all_sig_(&En_beam2,&mom_pi2,&theta_pi2,&radlen2,&type,&weight_f);
    weight_v = prot_prob*double(weight_f);
    type = 3;    
    wiser_all_sig_(&En_beam2,&mom_pi2,&theta_pi2,&radlen2,&type,&weight_f);
    weight_v = (1-prot_prob)*double(weight_f) + weight_v;        
  }
  else if (particle_id == 310 || particle_id == 130) {
    //  weight_v = 0.5 * ( WISER_ALL_SIG(En_beam,mom_pi,theta_pi,radlen,1) + WISER_ALL_SIG(En_beam,mom_pi,theta_pi,radlen,2)) ;
    type = 3;
    wiser_all_sig_(&En_beam2,&mom_pi2,&theta_pi2,&radlen2,&type,&weight_f);
    weight_v = 0.5 * double(weight_f);
    type = 4;
    wiser_all_sig_(&En_beam2,&mom_pi2,&theta_pi2,&radlen2,&type,&weight_f);
    weight_v = 0.5 * double(weight_f) + weight_v;
    weight_v = 0.5*weight_v;   ///Ks or Kl is only half of K0 which is average of K+ and K-
  }  
  else if (particle_id == 2212) {
    //   weight_v = WISER_ALL_FIT(mom_pi);
    //  weight_v = WISER_ALL_SIG(En_beam,mom_pi,theta_pi,radlen,1);
    type=5;
    wiser_all_sig_(&En_beam2,&mom_pi2,&theta_pi2,&radlen2,&type,&weight_f);
    weight_v = double(weight_f);
    //   cout << weight_f << " \t " << weight_v << endl;
    //    printf("%f \n",weight_v);
  }  
  else if (particle_id == -2212) {
    //   weight_v = WISER_ALL_FIT(mom_pi);
    //  weight_v = WISER_ALL_SIG(En_beam,mom_pi,theta_pi,radlen,1);
    type=6;
    wiser_all_sig_(&En_beam2,&mom_pi2,&theta_pi2,&radlen2,&type,&weight_f);
    weight_v = double(weight_f);
    //   cout << weight_f << " \t " << weight_v << endl;
    //    printf("%f \n",weight_v);
  }  
  else weight_v = 0.;  

  //printf("sigma = %e\n", weight_v*1e-3 );
  
  eventdata data;
  if( 0.0 < weight_v && weight_v < 1e9 ){
    data.weight  = weight_v * En_beam * 4 * TMath::Pi() * A; // nanobars/GeV-str-nuclei * (DeltaE sample generated) * (Full angle generated) * (Number of nucleons)
    //    data.weight *= beam->GetLumin();
    data.weight *= 1e-37 ; // nb to m^2

    data.weight *= beam->GetLumin();
	
  } else {
    // Unphysical for some reason
    data.weight = 0.0;
  }

  // Fixed target kinematics...  This needs to be done properly
  // for particle other than scattering electron, these vlaues are not meaningful
  data.Q2     = 0.0;
  data.x      = 0.0;
  data.y      = 0.0;
  data.W      = 0.0;
  data.Wpweight = 0.0;
  data.Wmweight = 0.0;   

  data.ef = pi_lab.E();
  data.pf = pi_lab.P();
  data.theta = pi_lab.Theta();
  data.phi = pi_lab.Phi();
  
  data.particle_id = particle_id;
  data.charge = charge;
  data.mass = mass;
  data.Z_ion = ion->GetZ();
  data.N_ion = ion->GetN();
  data.vx = vert.X();
  data.vy = vert.Y();
  data.vz = vert.Z();

  Gamma1.Rotate(TVector3(0,0,1).Angle(beam_ionrest.Vect()),TVector3(0,0,1).Cross(beam_ionrest.Vect()));   //wiser only takes theta relative to +z as assumed electron beam direction,need to rotate to actual electron beam direction    
  Gamma1.Boost(bv_ion_lab);
  Gamma2.Rotate(TVector3(0,0,1).Angle(beam_ionrest.Vect()),TVector3(0,0,1).Cross(beam_ionrest.Vect()));   //wiser only takes theta relative to +z as assumed electron beam direction,need to rotate to actual electron beam direction    
  Gamma2.Boost(bv_ion_lab);
  
  data.g1_theta = Gamma1.Theta();
  data.g1_phi   = Gamma1.Phi();
  data.g1_p     = Gamma1.P();

  data.g2_theta = Gamma2.Theta();
  data.g2_phi   = Gamma2.Phi();
  data.g2_p     = Gamma2.P();
  
  ev->SetEventData(data);

  return;
  
}

// e- elastic generator
void eicPhysics::MakeEvent3(eicBeam *beam, eicIon *ion, eicEvent *ev , eicModel *model) {

  double tglx = model->GetLx();
  double tgly = model->GetLy();
  double tglength = model->GetLength();

  TVector3 tgtoff = model->GetTgtOffset();


  // Initialize variables

  double ionp  = sqrt( pow(ion->GetEnergy(),2.0) - pow(ion->GetMass(),2.0) );
  double beta  = ionp/ion->GetEnergy();
  double gamma = 1.0/sqrt(1.0-beta*beta);
  
  double e_lab = beam->GetEnergy()*gamma*(1.0+beta);
  
  TLorentzVector e_rest( 0.0, 0.0, e_lab, e_lab);
  TVector3 blab(0.0, 0.0, beta);

  nucl n;
  double A = ((double) (ion->GetZ()+ion->GetN()));
  double prot_prob = ((double) ion->GetZ())/A;
  // Determine which type of nucleon we hit
  if( fRandom->Uniform() < prot_prob ){
    n = kProton;
  } else {
    n = kNeutron;
  }

 
  

  // this is the part that I got from my code
     
  double d_sig, tau, sin2, cos2, G_E2, G_M2, GEp, GMp, GMn, G_Mn2;                                                 
  double MV2=0.71;                                                                                                
  double mu_p=2.792847; //magnetic moment of the proton                                             
  double mu_n = 1.91; // magnetic moment of the neutron = (-1.91), but here I need the abs             
  double alpha= 7.297352*pow(10,-3); // alpha EM interaction                                               
  //  double m_e= 0.00051099892; // electron mass in GeV                                           
  double m_p= 0.93827203; // proton mass in GeV                                           
  double m_n= 0.93956536; // neutron mass in GeV                                                        
  double dsigdOdE_el_nb;                                                     
  double nbarn = 0.389 * pow(10,6); //  barn: (1 GeV)**-2 = 0.389e-3 barn
  double x_e=0.,y_e=0.,z_e=0.;

 
  // double mom_e = fRandom->Uniform(En_beam); 
  fRandom->Sphere(x_e,y_e,z_e,1.0); // generate random  vectors, uniformly distributed over the surface
  // of a sphere of radius mom_pi
  TVector3 v3_e(x_e,y_e,z_e);
  double theta_e = v3_e.Theta();
  double phi_e = v3_e.Phi();
  double E_2_v = 0;
  double q2_dis = 0;
  int type;
  switch( n ){
    case kProton:
      type = 1;
      break;
    case kNeutron:
      type = 2;
      break;
    default:
      type = 1;
      break;
  }

  // Initialize these
  sin2 = 1e-9;                                                                                
  cos2 = 1e-9;                                                       
     
  if (type==1) {                                                                 
    E_2_v = e_lab/ (1 + e_lab/m_p * (1 - cos(theta_e))); // value applying the delta of Dirac from eq A.13 in Mo Tsai                     
    q2_dis = 4 * e_lab * E_2_v * pow(sin(theta_e/2),2); // determine again q2_dis with the new E_2             
    GEp = 1. / pow( 1. + q2_dis / MV2 , 2.);                                                                     
    GMp = mu_p / pow( 1. + q2_dis / MV2 , 2.);                                                                    
    sin2 = pow(sin(theta_e/2.0),2);                                                                                
    cos2 = pow(cos(theta_e/2.0),2);                                                       
    G_E2 = pow(GEp,2);                                                                                         
    G_M2 = pow(GMp,2);                                               
    tau = q2_dis/4/pow(m_p,2); 

    d_sig = pow(alpha,2)/4./pow(e_lab,2)/pow(sin2,2)/(1.+2.*e_lab/m_p*sin2) *((G_E2+tau*G_M2)/(1.+tau)*cos2 + 2.*tau*G_M2*sin2); // this cross se 
    // ction is already d2sigma/dOdE for elastic scattering, eta and delta of equation of A.13 are already applied and the delta is already expanded 

                                                                        
    dsigdOdE_el_nb = d_sig*nbarn ; // cross section from elas.c is determined in mubarn, so I substitute hbc2 with nbarn for termining it in nanobarn                                                        
  }                                                                
  else {                                  
    E_2_v = e_lab/ (1 + e_lab/m_n * (1 - cos(theta_e))); // value modified for the neutron                    
    q2_dis = 4 * e_lab * E_2_v * pow(sin(theta_e/2),2); // value modified for the neutron     
    tau = q2_dis/4/pow(m_n,2); // value modified for the neutron                    
    GMn = mu_n /  pow( 1. + q2_dis / MV2 , 2.);       
    G_Mn2 = pow(GMn,2);                                          
    d_sig = pow(alpha,2)/4./pow(e_lab,2)/pow(sin2,2)/(1.+2.*e_lab/m_p*sin2) *((tau*G_Mn2)/(1.+tau)*cos2 + 2.*tau*G_Mn2*sin2); // same for proton, but now G_En = 0          
    dsigdOdE_el_nb = d_sig*nbarn ;                        
  }               

 
  // Generating the vertex randomly in the target
  TVector3 vert;
  double vert_x, vert_y,vert_z;

  vert_x = fRandom->Uniform((-tglx/2),(tglx/2))+tgtoff.X(); ;
  vert_y = fRandom->Uniform((-tgly/2),(tgly/2))+tgtoff.Y(); ;
  vert_z = fRandom->Uniform((-tglength/2),(tglength/2))+tgtoff.Z(); ;

  vert.SetXYZ(vert_x,vert_y,vert_z);
  
  eventdata data;
  
 // Volume term - volume we sample in is just a sphere, so 4pi
  double V = 4.0*3.14159;
 
  if(0.0 < dsigdOdE_el_nb && dsigdOdE_el_nb < 1e29 ){
    data.weight  = dsigdOdE_el_nb*(1e-37); // nanobarn to m^2
    data.weight *= beam->GetLumin();
    data.weight *= V*A;  // multiply by sampling volume and number of nucleons
  } else {
    // Unphysical for some reason
    data.weight = 0.0;
  }

  

  data.ef     = E_2_v;
  data.theta  = theta_e;
  data.phi    = phi_e;
  data.Q2     = q2_dis;
  data.mass    = 0.000511;
  data.particle_id = 11;
  data.charge = -1;
  data.pf =  E_2_v;
  data.Z_ion = ion->GetZ();
  data.N_ion = ion->GetN();

  data.x      = q2_dis/(2.0*MASS_P*(e_lab-E_2_v));
  data.y      = (e_lab-E_2_v)/e_lab;
  data.W      = sqrt(MASS_P*MASS_P + 2.0*MASS_P*(e_lab-E_2_v) -q2_dis);


  data.vx = vert.X();
  data.vy = vert.Y();
  data.vz = vert.Z();

  ev->SetEventData(data);

  return;
    


}

//e- Moller generator
void eicPhysics::MakeEvent4(eicBeam *beam, eicIon *ion, eicEvent *ev , eicModel *model) {
  const double me = 0.511e-3; // GeV
  const double alpha = 1.0/137.0;

    // Moller scattering
  double tglx = model->GetLx();
  double tgly = model->GetLy();
  double tglength = model->GetLength();

  TVector3 tgtoff = model->GetTgtOffset();

  double ionp  = sqrt( pow(ion->GetEnergy(),2.0) - pow(ion->GetMass(),2.0) );
  double beta_lab  = ionp/ion->GetEnergy();
  double gamma_lab = 1.0/sqrt(1.0-beta_lab*beta_lab);
  
  double e_lab = beam->GetEnergy()*gamma_lab*(1.0+beta_lab);

  TLorentzVector e_rest( 0.0, 0.0, e_lab, e_lab);
  TVector3 blab(0.0, 0.0, beta_lab);

  // CoM <-> lab boost
  double beta_com  = sqrt( (e_lab - me)/(e_lab + me) );
  double gamma_com = 1.0/sqrt(1.0 - beta_com*beta_com);

  double e_com = me*gamma_com;

  // we have Z electrons per nucleon to scatter off of
  double Z = ion->GetZ()+ion->GetN();

  // Sample in CoM frame for now?

  double thcom = acos(2.0*fRandom->Uniform()-1.0);
  double phcom = 2.0*TMath::Pi()*fRandom->Uniform();

  // Technically, there are two Moller electrons we need to track
  // Because we're just interested in rates, let's randomly choose
  // one and then multiply the weight factor by two

  double V = 4.0*TMath::Pi(); // sampling space for dOmega

  double sigma = alpha*alpha*pow(3.0+cos(thcom)*cos(thcom),2.0)/pow(sin(thcom),4.0)/(2.0*me*e_lab);
  double weight_v = sigma*(0.197*0.197*1e-30)*V; // GeV^-2 -> m^2
                                                //  here there would be a factor of 2 for undercounting pairs, 
						//  but we must also divide by two because we are double covering 
						//  phasespace because of identical particles

  double pperp = e_com*sin(thcom);
  double ppar  = e_com*cos(thcom);
  
  /////////////////////////////////////
  eventdata data;

  data.weight = weight_v*beam->GetLumin()*Z;

  TLorentzVector lef;

  double primary = 0.0;
  if( fRandom->Integer(2) == 0 ){
      // Choose primary electron
      lef = TLorentzVector( pperp*cos(phcom), pperp*sin(phcom), ppar, e_com );
      primary = 1.0;
  } else {
      lef = TLorentzVector( -pperp*cos(phcom), -pperp*sin(phcom),-ppar, e_com );
      primary = -1.0;
  }

  // Put back into lab frame
  lef.Boost(TVector3(0.0, 0.0,  beta_com));

  // Fixed target kinematics...  This needs to be done properly
  data.Q2     = 2.0*e_lab*lef.Energy()*(1.0-cos(lef.Theta()));
  data.y      = 1.0 - lef.Energy()/e_lab;

  // and collision frame
  lef.Boost(TVector3(0.0, 0.0, -beta_lab));

  //  We'll just hijack these variables since they don't mean
  //  much for these processes.  Sorry!  -SPR 8/29/12
  data.x      = thcom;
  data.W      = primary;

  data.ef = lef.Energy();
  data.pf = lef.P();
  data.theta = lef.Theta();
  data.phi = lef.Phi();
  
  data.particle_id = 11;
  data.charge = -1;
  data.mass = me;
  data.Z_ion = ion->GetZ();
  data.N_ion = ion->GetN();

  double vert_x, vert_y, vert_z;

  vert_x = fRandom->Uniform((-tglx/2),(tglx/2)) + tgtoff.X();
  vert_y = fRandom->Uniform((-tgly/2),(tgly/2)) + tgtoff.Y();
  vert_z = fRandom->Uniform((-tglength/2),(tglength/2))+ tgtoff.Z();;

  TVector3 vert(vert_x,vert_y,vert_z);

  data.vx = vert.X();
  data.vy = vert.Y();
  data.vz = vert.Z();

  ev->SetEventData(data);

}

void eicPhysics::MakeEvent5_findmax(eicBeam *beam, eicIon *ion, eicEvent *ev , eicModel *model) {
  double radlen=0., mass=0., weight_v;
  double weight_v_pip, weight_v_pim,weight_v_Kp,weight_v_Km,weight_v_p,weight_v_pbar;
  int particle_id, charge;
  int modelsig = model->GetModel();
  radlen = model->GetRadLen();
  double tglx = model->GetLx();
  double tgly = model->GetLy();
  double tglength = model->GetLength();


  TVector3 tgtoff = model->GetTgtOffset();


  // Initialize variables
  particle_id = 1e9;
  charge      = 1e9;

  TLorentzVector ion_lab(ion->GetMom()*sin(ion->GetTheta())*cos(ion->GetPhi()),ion->GetMom()*sin(ion->GetTheta())*sin(ion->GetPhi()),ion->GetMom()*cos(ion->GetTheta()),sqrt(pow(ion->GetMom(),2)+pow(ion->GetMass(),2)));

  TLorentzVector beam_lab(beam->GetMom()*sin(beam->GetTheta())*cos(beam->GetPhi()),beam->GetMom()*sin(beam->GetTheta())*sin(beam->GetPhi()),beam->GetMom()*cos(beam->GetTheta()),sqrt(pow(beam->GetMom(),2)+pow(beam->GetMass(),2)));
     
  TVector3 bv_ion_lab=ion_lab.BoostVector();
  TLorentzVector beam_ionrest=beam_lab;
  beam_ionrest.Boost(-bv_ion_lab);

  double A = ((double) (ion->GetZ()+ion->GetN()));
  double prot_prob = ((double) ion->GetZ())/A;
  
  double En_beam = beam_ionrest.E();   

  if(  En_beam < 0.3 ){ 
      fprintf( stderr, "ERROR:  Beam energy too low for implemented uniform generator\n");
      exit(1);
  }

  //  TF2 *func;  
  if (modelsig == 10) { // pi+
    particle_id = 211; 
    charge = +1;
    mass = 0.1396; // mass in GeV
    //    func = new TF2("sigma_pip",Wiser_func_pip,0, En_beam,0,360,2);
  } 
  else if (modelsig == 11) { //pi-
    particle_id = -211;
    charge = -1;
    mass = 0.1396; // mass in GeV
    //   func = new TF2("sigma_pip",Wiser_func_pip,0, En_beam,0,360,2);
  }
  else if (modelsig == 12) { //pi0
    particle_id = 111;
    charge = 0;
    mass = 0.1350; // mass in GeV
    //   func = new TF2("sigma_pip",Wiser_func_pi0,0, En_beam,0,360,2);
  }  
  else if (modelsig == 13) { //K+
    particle_id = 321;
    charge = +1;
    mass = 0.4937; // mass in GeV
    //   func = new TF2("sigma_pip",Wiser_func_pi0,0, En_beam,0,360,2);
  }  
  else if (modelsig == 14) { //K-
    particle_id = -321;
    charge = -1;
    mass = 0.4937; // mass in GeV
    //   func = new TF2("sigma_pip",Wiser_func_pi0,0, En_beam,0,360,2);
  }  
  else if (modelsig == 15) { //Ks
    particle_id = 310;
    charge = 0;
    mass = 0.4976; // mass in GeV
    //   func = new TF2("sigma_pip",Wiser_func_pi0,0, En_beam,0,360,2);
  }  
  else if (modelsig == 16) { //Kl
    particle_id = 130;
    charge = 0;
    mass = 0.4976; // mass in GeV
    //   func = new TF2("sigma_pip",Wiser_func_pi0,0, En_beam,0,360,2);
  }    
  else if (modelsig == 17) { //p
    particle_id = 2212;
    charge = +1;
    mass = 0.9383; // mass in GeV
    //   func = new TF2("sigma_pip",Wiser_func_pi0,0, En_beam,0,360,2);
  }  
  else if (modelsig == 18) { //p-bar
    particle_id = -2212;
    charge = -1;
    mass = 0.9383; // mass in GeV
    //   func = new TF2("sigma_pip",Wiser_func_pi0,0, En_beam,0,360,2);
  }

  double intrad = 2.0*log(En_beam/0.000511)/(137.0*3.14159);

//   double radratio = 3.0*intrad/(4.0*radlen);
  
//   targprop = radratio*
//       (sqrt( 1.0 + fRandom->Uniform()*( pow(radratio, -2.0) + 2.0/radratio ) ) - 1.0 );

  // This is fixed
//   double fracradlen = targprop*radlen*(4.0/3.0) + intrad; 

//   float radlen2 = float(fracradlen);
  
  float En_beam2 = float(En_beam);

  float targprop = 1.0;
  float radlen2=targprop*float(radlen)*(4.0/3.0) + intrad; //max always happens at full target length
  
  // Scan for maximum
  max =0.0;
  
  int npidx  = 300;
  int nthidx = 300;
  
  float max_mom_pi2=0,max_theta_pi2=0;
  for(int i = 0; i < npidx; i++ ){
      for(int j = 0; j < nthidx; j++ ){	
	  float mom_pi2, theta_pi2;	
	  // Scan around 2 GeV
// 	  mom_pi2   = 0.1*((double) i)*En_beam2/npidx + mass;  // These are guesses, but they work for E down to 0.3 GeV
// 	  theta_pi2 = (10.0*((double) j)/nthidx/En_beam2)*3.14159/180;
	  mom_pi2   = 0.9*((double) i)*En_beam2/npidx + mass;
	  theta_pi2 = (180.0*((double) j)/nthidx)*3.14159/180;
	  
	  weight_v = 0.0;

	  weight_v_pip = wiser_sigma( En_beam2, mom_pi2, theta_pi2, radlen2, 0);
	  weight_v_pim = wiser_sigma( En_beam2, mom_pi2, theta_pi2, radlen2, 1);
	  weight_v_Kp = wiser_sigma( En_beam2, mom_pi2, theta_pi2, radlen2, 2);
	  weight_v_Km = wiser_sigma( En_beam2, mom_pi2, theta_pi2, radlen2, 3);
	  weight_v_p = wiser_sigma( En_beam2, mom_pi2, theta_pi2, radlen2, 4);
	  weight_v_pbar = wiser_sigma( En_beam2, mom_pi2, theta_pi2, radlen2, 5);//Returns the cross section in units of nanobarns/GeV/str

	  if (particle_id == 211) {
	      weight_v = prot_prob*weight_v_pip + (1.0-prot_prob)*weight_v_pim;
	  }
	  else if (particle_id == -211) {
	      weight_v = prot_prob*weight_v_pim + (1.0-prot_prob)*weight_v_pip;
	  }
	  else if (particle_id == 111) {
	      weight_v = (weight_v_pip + weight_v_pim) *2.0;
	  }
	  else if (particle_id == 321) {
	      weight_v = prot_prob*weight_v_Kp + (1.0-prot_prob)*weight_v_Km;
	  }
	  else if (particle_id == -321) {
	      weight_v = prot_prob*weight_v_Km + (1.0-prot_prob)*weight_v_Kp;
	  }
	  else if (particle_id == 310 || particle_id == 130) {
	      weight_v = 0.5*(weight_v_Kp + weight_v_Km) /2.0;
	  }
	  else if (particle_id == 2212) {
	      weight_v = weight_v_p;
	  }
	  else if (particle_id == -2212) {
	      weight_v = weight_v_pbar;
	  }	  
	  
	  if( weight_v > max ){ 
	      max   = weight_v;	 
	      max_mom_pi2 = mom_pi2;
	      max_theta_pi2 = theta_pi2;
	  }
	  
      }
  }

  double scale = 1.1;  
  max *= scale;

  if( max < 0.0 ){ printf("Kinematics too close to threshold\n"); exit(1);}
  else printf("Found Max crossection %f at p = %f, th = %f\n", max, max_mom_pi2, max_theta_pi2);
  
  if( !fHaveTotalXs ){
      printf("Calculating total wiser cross section\n");
      fTotalXs_pip = wiser_total_sigma( En_beam2, intrad, radlen*4.0/3.0, 0);
      fTotalXs_pim = wiser_total_sigma( En_beam2, intrad, radlen*4.0/3.0, 1);
      fTotalXs_Kp = wiser_total_sigma( En_beam2, intrad, radlen*4.0/3.0, 2);
      fTotalXs_Km = wiser_total_sigma( En_beam2, intrad, radlen*4.0/3.0, 3);
      fTotalXs_p = wiser_total_sigma( En_beam2, intrad, radlen*4.0/3.0, 4);
      fTotalXs_pbar = wiser_total_sigma( En_beam2, intrad, radlen*4.0/3.0, 5);      
      printf("Calculated\n");
      fHaveTotalXs = true;
  } 	    
  
}

//hadron generator according to wiser (distributed according to crossection)
void eicPhysics::MakeEvent5(eicBeam *beam, eicIon *ion, eicEvent *ev , eicModel *model) {

  double radlen=0., mass=0., weight_v;
  double weight_v_pip, weight_v_pim,weight_v_Kp,weight_v_Km,weight_v_p,weight_v_pbar;
  int particle_id, charge;
  int modelsig = model->GetModel();
  radlen = model->GetRadLen();
  double tglx = model->GetLx();
  double tgly = model->GetLy();
  double tglength = model->GetLength();


  TVector3 tgtoff = model->GetTgtOffset();


  // Initialize variables
  particle_id = 1e9;
  charge      = 1e9;

  TLorentzVector ion_lab(ion->GetMom()*sin(ion->GetTheta())*cos(ion->GetPhi()),ion->GetMom()*sin(ion->GetTheta())*sin(ion->GetPhi()),ion->GetMom()*cos(ion->GetTheta()),sqrt(pow(ion->GetMom(),2)+pow(ion->GetMass(),2)));

  TLorentzVector beam_lab(beam->GetMom()*sin(beam->GetTheta())*cos(beam->GetPhi()),beam->GetMom()*sin(beam->GetTheta())*sin(beam->GetPhi()),beam->GetMom()*cos(beam->GetTheta()),sqrt(pow(beam->GetMom(),2)+pow(beam->GetMass(),2)));
     
  TVector3 bv_ion_lab=ion_lab.BoostVector();
  TLorentzVector beam_ionrest=beam_lab;
  beam_ionrest.Boost(-bv_ion_lab);

  double A = ((double) (ion->GetZ()+ion->GetN()));
  double prot_prob = ((double) ion->GetZ())/A;
  
  double En_beam = beam_ionrest.E();   

  if(  En_beam < 0.3 ){ 
      fprintf( stderr, "ERROR:  Beam energy too low for implemented uniform generator\n");
      exit(1);
  }

  //  TF2 *func;  
  if (modelsig == 10) { // pi+
    particle_id = 211; 
    charge = +1;
    mass = 0.1396; // mass in GeV
    //    func = new TF2("sigma_pip",Wiser_func_pip,0, En_beam,0,360,2);
  } 
  else if (modelsig == 11) { //pi-
    particle_id = -211;
    charge = -1;
    mass = 0.1396; // mass in GeV
    //   func = new TF2("sigma_pip",Wiser_func_pip,0, En_beam,0,360,2);
  }
  else if (modelsig == 12) { //pi0
    particle_id = 111;
    charge = 0;
    mass = 0.1350; // mass in GeV
    //   func = new TF2("sigma_pip",Wiser_func_pi0,0, En_beam,0,360,2);
  }  
  else if (modelsig == 13) { //K+
    particle_id = 321;
    charge = +1;
    mass = 0.4937; // mass in GeV
    //   func = new TF2("sigma_pip",Wiser_func_pi0,0, En_beam,0,360,2);
  }  
  else if (modelsig == 14) { //K-
    particle_id = -321;
    charge = -1;
    mass = 0.4937; // mass in GeV
    //   func = new TF2("sigma_pip",Wiser_func_pi0,0, En_beam,0,360,2);
  }  
  else if (modelsig == 15) { //Ks
    particle_id = 310;
    charge = 0;
    mass = 0.4976; // mass in GeV
    //   func = new TF2("sigma_pip",Wiser_func_pi0,0, En_beam,0,360,2);
  }  
  else if (modelsig == 16) { //Kl
    particle_id = 130;
    charge = 0;
    mass = 0.4976; // mass in GeV
    //   func = new TF2("sigma_pip",Wiser_func_pi0,0, En_beam,0,360,2);
  }    
  else if (modelsig == 17) { //p
    particle_id = 2212;
    charge = +1;
    mass = 0.9383; // mass in GeV
    //   func = new TF2("sigma_pip",Wiser_func_pi0,0, En_beam,0,360,2);
  }  
  else if (modelsig == 18) { //p-bar
    particle_id = -2212;
    charge = -1;
    mass = 0.9383; // mass in GeV
    //   func = new TF2("sigma_pip",Wiser_func_pi0,0, En_beam,0,360,2);
  }

  double intrad = 2.0*log(En_beam/0.000511)/(137.0*3.14159);
 
//   double radratio = 3.0*intrad/(4.0*radlen);

//   targprop = radratio*
//       (sqrt( 1.0 + fRandom->Uniform()*( pow(radratio, -2.0) + 2.0/radratio ) ) - 1.0 );

  // This is fixed
//   double fracradlen = targprop*radlen*(4.0/3.0) + intrad; 

//     float radlen2 = float(fracradlen);
  
  float En_beam2 = float(En_beam);

  float targprop;  
  float radlen2;
  
  int cnt = 0;
  double totalxs = 0.0;  
  double testval = 0.0; 
  
  TLorentzVector pi_ionrest(0,0,0,0),pi_lab(0,0,0,0);  
  do { 
      targprop = fRandom->Uniform();
      radlen2 = targprop*float(radlen)*(4.0/3.0) + intrad;
    
      double ef, mom_pi, theta_pi,phi_pi;
      mom_pi = fRandom->Uniform(En_beam2);

      theta_pi = acos( fRandom->Uniform(-1.0, 1.0) );
      phi_pi   = 2.0*TMath::Pi()*fRandom->Uniform();

      ef = sqrt(pow(mom_pi,2) + pow(mass,2));      
      
      pi_ionrest.SetPxPyPzE(mom_pi*sin(theta_pi)*cos(phi_pi),mom_pi*sin(theta_pi)*sin(phi_pi),mom_pi*cos(theta_pi),ef);
      pi_ionrest.Rotate(TVector3(0,0,1).Angle(beam_ionrest.Vect()),TVector3(0,0,1).Cross(beam_ionrest.Vect()));   //wiser only takes theta relative to +z as assumed electron beam direction,need to rotate to actual electron beam direction      
      pi_lab=pi_ionrest;
      pi_lab.Boost(bv_ion_lab);
    //   cout << pi_ionrest.P() << "\t" <<pi_ionrest.Theta() << "\t" << pi_ionrest.Phi() << "\t" << pi_ionrest.E() << endl; 
    //   cout << pi_lab.P() << "\t" <<pi_lab.Theta() << "\t" << pi_lab.Phi() << "\t" << pi_lab.E() << endl;

      double x_pi=0.,y_pi=0.,z_pi=0.;      
      x_pi = mom_pi*sin(theta_pi)*cos(phi_pi);
      y_pi = mom_pi*sin(theta_pi)*sin(phi_pi);
      z_pi = mom_pi*cos(theta_pi);
  
      // of a sphere of radius mom_pi
      TVector3 v3_pi(x_pi,y_pi,z_pi);

      float mom_pi2, theta_pi2;      
      mom_pi2 = float(mom_pi);
      theta_pi2 = float(theta_pi);
      
      weight_v_pip = wiser_sigma( En_beam2, mom_pi2, theta_pi2, radlen2, 0);
      weight_v_pim = wiser_sigma( En_beam2, mom_pi2, theta_pi2, radlen2, 1);
      weight_v_Kp = wiser_sigma( En_beam2, mom_pi2, theta_pi2, radlen2, 2);
      weight_v_Km = wiser_sigma( En_beam2, mom_pi2, theta_pi2, radlen2, 3);
      weight_v_p = wiser_sigma( En_beam2, mom_pi2, theta_pi2, radlen2, 4);
      weight_v_pbar = wiser_sigma( En_beam2, mom_pi2, theta_pi2, radlen2, 5);	  

      if (particle_id == 211) {
	  totalxs = prot_prob*fTotalXs_pip+ (1.0-prot_prob)*fTotalXs_pim;
	  weight_v = prot_prob*weight_v_pip + (1.0-prot_prob)*weight_v_pim;
      }
      else if (particle_id == -211) {
	  totalxs = prot_prob*fTotalXs_pim+ (1.0-prot_prob)*fTotalXs_pip;
	  weight_v = prot_prob*weight_v_pim + (1.0-prot_prob)*weight_v_pip;
      }
      else if (particle_id == 111) {
	      totalxs = (fTotalXs_pip + fTotalXs_pim)*2.0;	
	      weight_v = (weight_v_pip+weight_v_pim) *2.0;
      }
      else if (particle_id == 321) {
	  totalxs = prot_prob*fTotalXs_Kp+ (1.0-prot_prob)*fTotalXs_Km;
	  weight_v = prot_prob*weight_v_Kp + (1.0-prot_prob)*weight_v_Km;
      }
      else if (particle_id == -321) {
	  totalxs = prot_prob*fTotalXs_Km+ (1.0-prot_prob)*fTotalXs_Kp;
	  weight_v = prot_prob*weight_v_Km + (1.0-prot_prob)*weight_v_Kp;
      }
      else if (particle_id == 310 || particle_id == 130) {
	  totalxs  = 0.5*(fTotalXs_Kp + fTotalXs_Km)/2.0;		
	  weight_v = 0.5*(weight_v_Kp + weight_v_Km) /2.0;
      }
      else if (particle_id == 2212) {
	  totalxs  = fTotalXs_p;			
	  weight_v = weight_v_p;
      }
      else if (particle_id == -2212) {
	  totalxs  = fTotalXs_pbar;				
	  weight_v = weight_v_pbar;
      }      
      
      cnt++;

      double v = fRandom->Uniform();

      testval = max*v;
      if( weight_v > max ){
	  printf("out of bounds in uniform part %f > %f (p = %f, th = %f), targprop = %f \n", weight_v, max, mom_pi, theta_pi,targprop);
	  exit(1);
      }

  } while ( weight_v < testval );

//  printf("Number of trials = %d, max = %f, weight_v = %g\n", cnt, max, weight_v );

  // Generating the vertex randomly in the target
  TVector3 vert;
  double vert_x, vert_y,vert_z;

  vert_x = fRandom->Uniform((-tglx/2),(tglx/2))+tgtoff.X(); ;
  vert_y = fRandom->Uniform((-tgly/2),(tgly/2))+tgtoff.Y(); ;
  vert_z = (targprop-0.5)*tglength+tgtoff.Z();

  vert.SetXYZ(vert_x,vert_y,vert_z);
  
  eventdata data;
  if( 0.0 < weight_v && weight_v < 1e9 ){
    data.weight  = totalxs * A; // nanobars/GeV-str-nuclei * (DeltaE sample generated) * (Full angle generated) * (Number of nucleons)
    //    data.weight *= beam->GetLumin();
    data.weight *= 1e-37 ; // nb to m^2

    data.weight *= beam->GetLumin();
	
  } else {
    // Unphysical for some reason

    data.weight = 0.0;
  }

  // Fixed target kinematics...  This needs to be done properly
  // for particle other than scattering electron, these vlaues are not meaningful
  data.Q2     = 0.0;
  data.x      = 0.0;
  data.y      = 0.0;
  data.W      = 0.0;
  data.Wpweight = 0.0;
  data.Wmweight = 0.0;   

  data.ef = pi_lab.E();
  data.pf = pi_lab.P();
  data.theta = pi_lab.Theta();
  data.phi = pi_lab.Phi();
  
  data.particle_id = particle_id;
  data.charge = charge;
  data.mass = mass;
  data.Z_ion = ion->GetZ();
  data.N_ion = ion->GetN();
  data.vx = vert.X();
  data.vy = vert.Y();
  data.vz = vert.Z();

  Gamma1.Rotate(TVector3(0,0,1).Angle(beam_ionrest.Vect()),TVector3(0,0,1).Cross(beam_ionrest.Vect()));   //wiser only takes theta relative to +z as assumed electron beam direction,need to rotate to actual electron beam direction    
  Gamma1.Boost(bv_ion_lab);
  Gamma2.Rotate(TVector3(0,0,1).Angle(beam_ionrest.Vect()),TVector3(0,0,1).Cross(beam_ionrest.Vect()));   //wiser only takes theta relative to +z as assumed electron beam direction,need to rotate to actual electron beam direction    
  Gamma2.Boost(bv_ion_lab);
  
  data.g1_theta = Gamma1.Theta();
  data.g1_phi   = Gamma1.Phi();
  data.g1_p     = Gamma1.P();

  data.g2_theta = Gamma2.Theta();
  data.g2_phi   = Gamma2.Phi();
  data.g2_p     = Gamma2.P();
  
  ev->SetEventData(data);
    

  return;
  
}


//===================Ziheng==============
//for Peter Bosted fit, model Number 23


void eicPhysics::MakeEvent6(eicBeam *beam, eicIon *ion, eicEvent *ev , eicModel *model) {
    
    G2PPhysBase *pPhys =  new G2PPhysPB();
    
    // Radiative correction parameters for P.B. fit: RadLength[2]=(Tb,Ta);
    // Tb: total radiative length before scattering in radiation length;
    // Ta: total radiative length after scattering in radiation length;
    double RadLength[2]={0,0};  // Set Radiative correction. None here.
    const double me = 0.000511; // Set electron mass in GeV
    const double Mp = 0.93825;
    const double Mn = 0.93957;
    
//1. tell pPhy the beam and target
    // Beam: BeamParticle #11(electron), radiative correction {Tb=0,Ta=0};
    // Target: Atom# Z, nucleon# A, Mass M;
    
    // ??we have A electrons per nucleon to scatter off of
    // 11 for electron beam.
    int A = (ion->GetZ()+ion->GetN());
    int Z = ion->GetZ();
    double M = Mp;
    //the same mass used in Bosted.f line 91
    pPhys->SetParticle(11);
    pPhys->SetPars(RadLength,2);
    //"2" for two parameters Ta and Tb; or "0" for nothing
    pPhys->SetTarget(Z, A);
//    pPhys->SetTargetMass(M);
    
//2.get uniform vert_xyz
    //what physics meaning of tglx, vert_x, tftoff?
    double tglx = model->GetLx();
    double tgly = model->GetLy();
    double tglength = model->GetLength();
    TVector3 tgtoff = model->GetTgtOffset();
    
    double vert_x = fRandom->Uniform((-tglx/2),(tglx/2)) + tgtoff.X();
    double vert_y = fRandom->Uniform((-tgly/2),(tgly/2)) + tgtoff.Y();
    double vert_z = fRandom->Uniform((-tglength/2),(tglength/2))+ tgtoff.Z();
    //printf("vx=%f, vy=%f, vz = %f \n", vert_x,vert_y,vert_z);
    
    
//3.calculate crossection in given Ei, Ef, theta
    double Ei = beam->GetEnergy();
    double Ef = Ei*fRandom->Uniform(0,1.0);
    double Pf = sqrt(Ef*Ef-me*me);
    double phi = 2.0*TMath::Pi()*fRandom->Uniform(0,1.0);
    double theta = acos(fRandom->Uniform(-1,1));
    //printf("Ei=%f, Ef=%f, Pf=%f, phi=%f,theta=%f \n", Ei, Ef, Pf, phi, theta);
    
    
    double Q2 = 4.*Ei*Ef*pow((sin(abs(theta)/2.)),2);
    double W2 = M*M + 2*M*( Ei - Ef ) - Q2;
    double x = Q2/(2*M*(Ei-Ef));
    double y = 1.0 - Ef/Ei;
    //printf("x=%f,Q2=%f, w2=%f \n", x,Q2,W2);

    
    double V = 4.0*TMath::Pi()*(Ei-0);
    // V has unit of GeV*sr
    double sigma = 0;
    double weight_v = 0;
 
    
//4. Use P Bosted fit
    
/*-------------------------------------------------------------------------------------
|   Fit to inelastic cross sections for A(e,e')X
|   valid for all W<3 GeV and all Q2<10 GeV2
|
|   Inputs: Z, A (double precision)are Z and A of nucleus
|                                  (use Z=0., A=1. to get free neutron)
|           Qsq (double precision) is 4-vector momentum transfer squared
|                                  (positive in chosen metric)
|           Wsq (double precision) is invarinat mass squared of final
|                                  state calculated assuming electron scattered from
|                                  a free proton
|
|   outputs: F1, F2 (double precision) are structure functions per nucleus
|   
|   Version of 10/20/2006 P. Bosted
--------------------------------------------------------------------------------------*/
    
    
    if ( W2>0 )
    {
        sigma = pPhys->GetXS(Ei,Ef,theta)*1e6;
        //cross section in units of nb (1e-33 cm^2, 1e-37 m^2);
        //the GetXS return in ub/MeV-sr of total nucleus
        //1e6 from ub/MeV-sr to nb/GeV-sr
        // sigma in unit of nb/GeV-sr
        weight_v = sigma*beam->GetLumin()*V*1e-37;
        //1e-37 for nb to m^2;
        //luminosity Hz/m^2
        if (std::isnan(sigma)) {sigma = 0;weight_v = 0;}
    }
    else W2=0;
    // for W2<0 non-physical event, set W2 = 0;
    
    //printf(" sigma=%f, weight_v=%f \n", sigma, weight_v);


 //5. write the date sheet
    eventdata data;
    
        data.weight = weight_v;
        data.crs = sigma;
        data.Q2     = Q2;
        data.x      = x;
        data.y      = y;//should be y; change to debug
        data.W      = sqrt(W2);
        data.ef     = Ef;
        data.pf     = Pf;
        data.theta  = theta;
        data.phi    = phi;
        data.particle_id = 11;
        data.charge = -1;
        data.Z_ion  = ion->GetZ();
        data.N_ion  = ion->GetN();
        data.vx     = vert_x;
        data.vy     = vert_y;
        data.vz     = vert_z;
    
    ev->SetEventData(data);

    }

//=====================END====================



void eicPhysics::SampleWiserPThZ( double En_beam, double &p, double &th, double &z, double *pars ){
  p = fRandom->Uniform(En_beam); 
  th = acos(fRandom->Uniform(-1.0, 1.0)); 
  z = fRandom->Uniform(); 

  return;
}

double eicPhysics::F1( double x, double Q2, nucl n ){
    return F2(x,Q2,n)/(2.0*x);
}

double eicPhysics::F2( double x, double Q2, nucl n ){
    // Sanity checks
    if( x < 0 || 1.0 < x ) return 0.0;
    if( Q2 < 0 ) { printf("Negative Q2?!\n"); exit(1); }
    // To keep things nice, this puts us in the resonance
    // region so there's nothing useful here anyways
    if( Q2 < 1.0 ) { Q2 = 1.0; }

    double e2_u = pow( 2.0/3.0,2.0);
    double e2_d = pow(-1.0/3.0,2.0);
    double f2sum = 0.0;

    int uidx, didx;

    switch( n ){
	case kProton:
	    uidx = 1;
	    didx = 2;
	    break;
	case kNeutron:
    	    uidx = 2;
	    didx = 1;
	    break;
	default:
    	    uidx = 1;
	    didx = 2;
	    break;
    }

    f2sum += e2_u*x*cteq_pdf_evolvepdf(fPDF, uidx, x, sqrt(Q2));
    f2sum += e2_d*x*cteq_pdf_evolvepdf(fPDF, didx, x, sqrt(Q2));

    // Strange quark
    f2sum += e2_d*x*cteq_pdf_evolvepdf(fPDF,  3, x, sqrt(Q2));
    f2sum += e2_d*x*cteq_pdf_evolvepdf(fPDF, -3, x, sqrt(Q2));

    // Charm quark
    f2sum += e2_u*x*cteq_pdf_evolvepdf(fPDF,  4, x, sqrt(Q2));
    f2sum += e2_u*x*cteq_pdf_evolvepdf(fPDF, -4, x, sqrt(Q2));

    // ubar
    f2sum += e2_u*x*cteq_pdf_evolvepdf(fPDF, -uidx, x, sqrt(Q2));
    // dbar
    f2sum += e2_d*x*cteq_pdf_evolvepdf(fPDF, -didx, x, sqrt(Q2));

    return f2sum;
}

double eicPhysics::F1gZ( double x, double Q2, nucl n ){
    // Sanity checks
    if( x < 0 || 1.0 < x ) return 0.0;
    if( Q2 < 0 ) { printf("Negative Q2?!\n"); exit(1); }
    // To keep things nice, this puts us in the resonance
    // region so there's nothing useful here anyways
    if( Q2 < 1.0 ) { Q2 = 1.0; }

    double e_u =  2.0/3.0;
    double e_d = -1.0/3.0;

    double sin2thW = 0.23119;
    double gV_u =  0.5 - 4.0*sin2thW/3.0;
    double gV_d = -0.5 + 2.0*sin2thW/3.0;

    double f1sum = 0.0;

    int uidx, didx;

    switch( n ){
	case kProton:
	    uidx = 1;
	    didx = 2;
	    break;
	case kNeutron:
    	    uidx = 2;
	    didx = 1;
	    break;
	default:
    	    uidx = 1;
	    didx = 2;
	    break;
    }

    f1sum += e_u*gV_u*cteq_pdf_evolvepdf(fPDF, uidx, x, sqrt(Q2));
    f1sum += e_d*gV_d*cteq_pdf_evolvepdf(fPDF, didx, x, sqrt(Q2));

    // Strange quark
    f1sum += e_d*gV_d*cteq_pdf_evolvepdf(fPDF,  3, x, sqrt(Q2));
    f1sum += e_d*gV_d*cteq_pdf_evolvepdf(fPDF, -3, x, sqrt(Q2));

    // Charm quark
    f1sum += e_u*gV_u*cteq_pdf_evolvepdf(fPDF,  4, x, sqrt(Q2));
    f1sum += e_u*gV_u*cteq_pdf_evolvepdf(fPDF, -4, x, sqrt(Q2));

    // ubar
    f1sum += e_u*gV_u*cteq_pdf_evolvepdf(fPDF, -uidx, x, sqrt(Q2));
    // dbar
    f1sum += e_d*gV_d*cteq_pdf_evolvepdf(fPDF, -didx, x, sqrt(Q2));

    return f1sum;
}

double eicPhysics::F3gZ( double x, double Q2, nucl n ){
    // Sanity checks
    if( x < 0 || 1.0 < x ) return 0.0;
    if( Q2 < 0 ) { printf("Negative Q2?!\n"); exit(1); }
    // To keep things nice, this puts us in the resonance
    // region so there's nothing useful here anyways
    if( Q2 < 1.0 ) { Q2 = 1.0; }

    double e_u =  2.0/3.0;
    double e_d = -1.0/3.0;

    double gA_u =  1.0/2.0;
    double gA_d = -1.0/2.0;

    double f3sum = 0.0;


    int uidx, didx;

    switch( n ){
	case kProton:
	    uidx = 1;
	    didx = 2;
	    break;
	case kNeutron:
    	    uidx = 2;
	    didx = 1;
	    break;
	default:
    	    uidx = 1;
	    didx = 2;
	    break;
    }

    f3sum += e_u*gA_u*cteq_pdf_evolvepdf(fPDF, uidx, x, sqrt(Q2));
    f3sum += e_d*gA_d*cteq_pdf_evolvepdf(fPDF, didx, x, sqrt(Q2));

    // Strange quark
    f3sum += e_d*gA_d*cteq_pdf_evolvepdf(fPDF,  3, x, sqrt(Q2));
    f3sum -= e_d*gA_d*cteq_pdf_evolvepdf(fPDF, -3, x, sqrt(Q2));

   // Charm quark
    f3sum += e_u*gA_u*cteq_pdf_evolvepdf(fPDF,  4, x, sqrt(Q2));
    f3sum -= e_u*gA_u*cteq_pdf_evolvepdf(fPDF, -4, x, sqrt(Q2));

    // ubar
    f3sum -= e_u*gA_u*cteq_pdf_evolvepdf(fPDF, -uidx, x, sqrt(Q2));
    // dbar
    f3sum -= e_d*gA_d*cteq_pdf_evolvepdf(fPDF, -didx, x, sqrt(Q2));

    return f3sum*2.0;
}

double eicPhysics::F1Z( double x, double Q2, nucl n ){
    // Sanity checks
    if( x < 0 || 1.0 < x ) return 0.0;
    if( Q2 < 0 ) { printf("Negative Q2?!\n"); exit(1); }
    // To keep things nice, this puts us in the resonance
    // region so there's nothing useful here anyways
    if( Q2 < 1.0 ) { Q2 = 1.0; }

    double sin2thW = 0.23119;
    double gA_u =  0.5;
    double gA_d = -0.5;
    double gV_u =  0.5 - 4.0*sin2thW/3.0;
    double gV_d = -0.5 + 2.0*sin2thW/3.0;

    double f1sum = 0.0;

    int uidx, didx;

    switch( n ){
	case kProton:
	    uidx = 1;
	    didx = 2;
	    break;
	case kNeutron:
    	    uidx = 2;
	    didx = 1;
	    break;
	default:
    	    uidx = 1;
	    didx = 2;
	    break;
    }

    f1sum += 0.5*(gV_u*gV_u+gA_u*gA_u)*cteq_pdf_evolvepdf(fPDF, uidx, x, sqrt(Q2));
    f1sum += 0.5*(gV_d*gV_d+gA_d*gA_d)*cteq_pdf_evolvepdf(fPDF, didx, x, sqrt(Q2));

    // Strange quark
    f1sum += 0.5*(gV_d*gV_d+gA_d*gA_d)*cteq_pdf_evolvepdf(fPDF,  3, x, sqrt(Q2));
    f1sum += 0.5*(gV_d*gV_d+gA_d*gA_d)*cteq_pdf_evolvepdf(fPDF, -3, x, sqrt(Q2));
    // Charm quark
    f1sum += 0.5*(gV_u*gV_u+gA_u*gA_u)*cteq_pdf_evolvepdf(fPDF,  4, x, sqrt(Q2));
    f1sum += 0.5*(gV_u*gV_u+gA_u*gA_u)*cteq_pdf_evolvepdf(fPDF, -4, x, sqrt(Q2));

    // ubar
    f1sum += 0.5*(gV_u*gV_u+gA_u*gA_u)*cteq_pdf_evolvepdf(fPDF, -uidx, x, sqrt(Q2));
    // dbar
    f1sum += 0.5*(gV_d*gV_d+gA_d*gA_d)*cteq_pdf_evolvepdf(fPDF, -didx, x, sqrt(Q2));

    return f1sum;
}

double eicPhysics::F3Z( double x, double Q2, nucl n ){
    // Sanity checks
    if( x < 0 || 1.0 < x ) return 0.0;
    if( Q2 < 0 ) { printf("Negative Q2?!\n"); exit(1); }
    // To keep things nice, this puts us in the resonance
    // region so there's nothing useful here anyways
    if( Q2 < 1.0 ) { Q2 = 1.0; }

    double sin2thW = 0.23119;

    double gA_u =  1.0/2.0;
    double gA_d = -1.0/2.0;
    double gV_u =  0.5 - 4.0*sin2thW/3.0;
    double gV_d = -0.5 + 2.0*sin2thW/3.0;

    double f3sum = 0.0;


    int uidx, didx;

    switch( n ){
	case kProton:
	    uidx = 1;
	    didx = 2;
	    break;
	case kNeutron:
    	    uidx = 2;
	    didx = 1;
	    break;
	default:
    	    uidx = 1;
	    didx = 2;
	    break;
    }

    f3sum += gV_u*gA_u*cteq_pdf_evolvepdf(fPDF, uidx, x, sqrt(Q2));
    f3sum += gV_d*gA_d*cteq_pdf_evolvepdf(fPDF, didx, x, sqrt(Q2));

    // Strange quark
    f3sum += gV_d*gA_d*cteq_pdf_evolvepdf(fPDF,  3, x, sqrt(Q2));
    f3sum -= gV_d*gA_d*cteq_pdf_evolvepdf(fPDF, -3, x, sqrt(Q2));
    
    // Charm quark
    f3sum += gV_u*gA_u*cteq_pdf_evolvepdf(fPDF,  4, x, sqrt(Q2));
    f3sum -= gV_u*gA_u*cteq_pdf_evolvepdf(fPDF, -4, x, sqrt(Q2));

    // ubar
    f3sum -= gV_u*gA_u*cteq_pdf_evolvepdf(fPDF, -uidx, x, sqrt(Q2));
    // dbar
    f3sum -= gV_d*gA_d*cteq_pdf_evolvepdf(fPDF, -didx, x, sqrt(Q2));

    return f3sum*2.0;
}

double eicPhysics::GetTheta( double amin, double amax ){
    double min = cos(amax);
    double max = cos(amin);

    double cth = fRandom->Uniform(min,max);
    return acos(cth);
}

double eicPhysics::GetPhi( double min, double max ){
    return fRandom->Uniform(min,max);
}

double eicPhysics::GetEf( double en, double th ){
    double Emax = (MASS_P*en)/(MASS_P + en*(1.0-cos(th)));
    
    return fRandom->Uniform(0.0, Emax);
}

double eicPhysics::GetEfmax( double en, double th ){
    double Emax = (MASS_P*en)/(MASS_P + en*(1.0-cos(th)));
    return Emax;
}

double eicPhysics::g1gZ( double x, double Q2, nucl n ){
    // Sanity checks
    if( x < 0 || 1.0 < x ) return 0.0;
    if( Q2 < 0 ) { printf("Negative Q2?!\n"); exit(1); }
    // To keep things nice, this puts us in the resonance
    // region so there's nothing useful here anyways
    if( Q2 < 1.0 ) { Q2 = 1.0; }

    double e_u =  2.0/3.0;
    double e_d = -1.0/3.0;

    double sin2thW = 0.23119;
    double gV_u =  0.5 - 4.0*sin2thW/3.0;
    double gV_d = -0.5 + 2.0*sin2thW/3.0;
    //double gA_u =  1.0/2.0;
    //double gA_d = -1.0/2.0;

    double g1sum = 0.0;

    int uidx, didx;

    switch( n ){
	case kProton:
	    uidx = 0;
	    didx = 1;
	    break;
	case kNeutron:
    	    uidx = 1;
	    didx = 0;
	    break;
	default:
    	    uidx = 0;
	    didx = 1;
	    break;
    }

    g1sum += e_u*gV_u*getDeltaq(x, Q2, uidx);
    g1sum += e_d*gV_d*getDeltaq(x, Q2, didx);

    // Strange quark
    g1sum += 2.0*e_d*gV_d*getDeltaq(x, Q2, 4);
 
    // ubar
    g1sum += e_u*gV_u*getDeltaq(x, Q2, uidx+2);
    // dbar
    g1sum += e_d*gV_d*getDeltaq(x, Q2, didx+2);
    
    return g1sum;
}

double eicPhysics::g5gZ( double x, double Q2, nucl n ){
// Sanity checks
    if( x < 0 || 1.0 < x ) return 0.0;
    if( Q2 < 0 ) { printf("Negative Q2?!\n"); exit(1); }
    // To keep things nice, this puts us in the resonance
    // region so there's nothing useful here anyways
    if( Q2 < 1.0 ) { Q2 = 1.0; }

    double e_u =  2.0/3.0;
    double e_d = -1.0/3.0;

    //double sin2thW = 0.23119;
    //double gV_u =  0.5 - 4.0*sin2thW/3.0;
    //double gV_d = -0.5 + 2.0*sin2thW/3.0;
    double gA_u =  1.0/2.0;
    double gA_d = -1.0/2.0;

    double g5sum = 0.0;

    int uidx, didx;

    switch( n ){
	case kProton:
	    uidx = 0;
	    didx = 1;
	    break;
	case kNeutron:
    	    uidx = 1;
	    didx = 0;
	    break;
	default:
    	    uidx = 0;
	    didx = 1;
	    break;
    }

    g5sum += e_u*gA_u*getDeltaq(x, Q2, uidx);
    g5sum += e_d*gA_d*getDeltaq(x, Q2, didx);

    // ubar
    g5sum -= e_u*gA_u*getDeltaq(x, Q2, uidx+2);
    // dbar
    g5sum -= e_d*gA_d*getDeltaq(x, Q2, didx+2);

    return g5sum;
}


double eicPhysics::g1Z( double x, double Q2, nucl n ){
    // Sanity checks
    if( x < 0 || 1.0 < x ) return 0.0;
    if( Q2 < 0 ) { printf("Negative Q2?!\n"); exit(1); }
    // To keep things nice, this puts us in the resonance
    // region so there's nothing useful here anyways
    if( Q2 < 1.0 ) { Q2 = 1.0; }

    //double e_u =  2.0/3.0;
    //double e_d = -1.0/3.0;

    double sin2thW = 0.23119;
    double gV_u =  0.5 - 4.0*sin2thW/3.0;
    double gV_d = -0.5 + 2.0*sin2thW/3.0;
    double gA_u =  1.0/2.0;
    double gA_d = -1.0/2.0;

    double g1sum = 0.0;

    int uidx, didx;

    switch( n ){
	case kProton:
	    uidx = 0;
	    didx = 1;
	    break;
	case kNeutron:
    	    uidx = 1;
	    didx = 0;
	    break;
	default:
    	    uidx = 0;
	    didx = 1;
	    break;
    }

    g1sum += (gV_u*gV_u+gA_u*gA_u)*getDeltaq(x, Q2, uidx);
    g1sum += (gV_d*gV_d+gA_d*gA_d)*getDeltaq(x, Q2, didx);

    // Strange quark
    g1sum += 2.0*(gV_d*gV_d+gA_d*gA_d)*getDeltaq(x, Q2, 4);
 
    // ubar
    g1sum += (gV_u*gV_u+gA_u*gA_u)*getDeltaq(x, Q2, uidx+2);
    // dbar
    g1sum += (gV_d*gV_d+gA_d*gA_d)*getDeltaq(x, Q2, didx+2);
    
    return g1sum/2.0;
}

double eicPhysics::g5Z( double x, double Q2, nucl n ){
// Sanity checks
    if( x < 0 || 1.0 < x ) return 0.0;
    if( Q2 < 0 ) { printf("Negative Q2?!\n"); exit(1); }
    // To keep things nice, this puts us in the resonance
    // region so there's nothing useful here anyways
    if( Q2 < 1.0 ) { Q2 = 1.0; }

    //double e_u =  2.0/3.0;
    //double e_d = -1.0/3.0;

    double sin2thW = 0.23119;
    double gV_u =  0.5 - 4.0*sin2thW/3.0;
    double gV_d = -0.5 + 2.0*sin2thW/3.0;
    double gA_u =  1.0/2.0;
    double gA_d = -1.0/2.0;

    double g5sum = 0.0;

    int uidx, didx;

    switch( n ){
	case kProton:
	    uidx = 0;
	    didx = 1;
	    break;
	case kNeutron:
    	    uidx = 1;
	    didx = 0;
	    break;
	default:
    	    uidx = 0;
	    didx = 1;
	    break;
    }

    g5sum += gV_u*gA_u*getDeltaq(x, Q2, uidx);
    g5sum += gV_d*gA_d*getDeltaq(x, Q2, didx);

    // ubar
    g5sum -= gV_u*gA_u*getDeltaq(x, Q2, uidx+2);
    // dbar
    g5sum -= gV_d*gA_d*getDeltaq(x, Q2, didx+2);

    return g5sum;
}


double eicPhysics::F1W( double x, double Q2, nucl n, int C ){
    return F2W(x,Q2,n,C)/(2.0*x);
}

double eicPhysics::F2W( double x, double Q2, nucl n, int C ){
    // Sanity checks
    if( x < 0 || 1.0 < x ) return 0.0;
    if( Q2 < 0 ) { printf("Negative Q2?!\n"); exit(1); }
    // To keep things nice, this puts us in the resonance
    // region so there's nothing useful here anyways
    if( Q2 < 1.0 ) { Q2 = 1.0; }

    double f2sum = 0.0;

    int uidx, didx, sidx, cidx;

    switch( n ){
	case kProton:
	    uidx = 1;
	    didx = 2;
	    sidx = 3;
	    cidx = 4;
	    break;
	case kNeutron:
    	    uidx = 2;
	    didx = 1;
	    sidx = 3;
	    cidx = 4;
	    break;
	default:
    	    uidx = 1;
	    didx = 2;
	    sidx = 3;
	    cidx = 4;
	    break;
    }

    assert( C==1 || C==-1 );

    // Swap quarks if we have positive charged current
    int tmp;
    if( C == 1 ){
	tmp = didx;
	didx = uidx;
	uidx = tmp;

    	tmp = cidx;
	cidx = sidx;
	sidx = tmp;
    }

    f2sum += 2.0*x*cteq_pdf_evolvepdf(fPDF, uidx, x, sqrt(Q2));
    f2sum += 2.0*x*cteq_pdf_evolvepdf(fPDF,-didx, x, sqrt(Q2));

    // Strange quark
    f2sum += 2.0*x*cteq_pdf_evolvepdf(fPDF, -sidx, x, sqrt(Q2));
    f2sum += 2.0*x*cteq_pdf_evolvepdf(fPDF,  cidx, x, sqrt(Q2));

    return f2sum;
}

double eicPhysics::F3W( double x, double Q2, nucl n, int C ){
    // Sanity checks
    if( x < 0 || 1.0 < x ) return 0.0;
    if( Q2 < 0 ) { printf("Negative Q2?!\n"); exit(1); }
    // To keep things nice, this puts us in the resonance
    // region so there's nothing useful here anyways
    if( Q2 < 1.0 ) { Q2 = 1.0; }

    double f3sum = 0.0;

    int uidx, didx, sidx, cidx;

    switch( n ){
	case kProton:
	    uidx = 1;
	    didx = 2;
	    sidx = 3;
	    cidx = 4;
	    break;
	case kNeutron:
    	    uidx = 2;
	    didx = 1;
	    sidx = 3;
	    cidx = 4;
	    break;
	default:
    	    uidx = 1;
	    didx = 2;
	    sidx = 3;
	    cidx = 4;
	    break;
    }

    assert( C==1 || C==-1 );

    // Swap quarks if we have positive charged current
    int tmp;
    if( C == 1 ){
	tmp = didx;
	didx = uidx;
	uidx = tmp;

    	tmp = cidx;
	cidx = sidx;
	sidx = tmp;
    }


    f3sum += 2.0*cteq_pdf_evolvepdf(fPDF, uidx, x, sqrt(Q2));
    f3sum -= 2.0*cteq_pdf_evolvepdf(fPDF,-didx, x, sqrt(Q2));

    // Strange quark
    f3sum -= 2.0*cteq_pdf_evolvepdf(fPDF, -sidx, x, sqrt(Q2));
    f3sum += 2.0*cteq_pdf_evolvepdf(fPDF,  cidx, x, sqrt(Q2));

    return f3sum;
}

double eicPhysics::g1W( double x, double Q2, nucl n, int C ){
    // Sanity checks
    if( x < 0 || 1.0 < x ) return 0.0;
    if( Q2 < 0 ) { printf("Negative Q2?!\n"); exit(1); }
    // To keep things nice, this puts us in the resonance
    // region so there's nothing useful here anyways
    if( Q2 < 1.0 ) { Q2 = 1.0; }

    double g1sum = 0.0;

    int uidx, didx;

    switch( n ){
	case kProton:
	    uidx = 0;
	    didx = 1;
	    break;
	case kNeutron:
    	    uidx = 1;
	    didx = 0;
	    break;
	default:
    	    uidx = 0;
	    didx = 1;
	    break;
    }

    assert( C==1 || C==-1 );

    int tmp;
    if( C == 1 ){
	tmp = didx;
	didx = uidx;
	uidx = tmp;
    }

    g1sum += getDeltaq(x, Q2, uidx);
    g1sum += getDeltaq(x, Q2, didx+2);
    g1sum += getDeltaq(x, Q2, 4);

    return g1sum;
}

double eicPhysics::g4W( double x, double Q2, nucl n, int C ){
    return 2.0*x*g5W(x,Q2,n,C);
}

double eicPhysics::g5W( double x, double Q2, nucl n, int C ){
// Sanity checks
    if( x < 0 || 1.0 < x ) return 0.0;
    if( Q2 < 0 ) { printf("Negative Q2?!\n"); exit(1); }
    // To keep things nice, this puts us in the resonance
    // region so there's nothing useful here anyways
    if( Q2 < 1.0 ) { Q2 = 1.0; }


    double g5sum = 0.0;

    int uidx, didx;

    switch( n ){
	case kProton:
	    uidx = 0;
	    didx = 1;
	    break;
	case kNeutron:
    	    uidx = 1;
	    didx = 0;
	    break;
	default:
    	    uidx = 0;
	    didx = 1;
	    break;
    }

    int tmp;
    if( C == 1 ){
	tmp = didx;
	didx = uidx;
	uidx = tmp;
    }

    g5sum -= getDeltaq(x, Q2, uidx);
    g5sum += getDeltaq(x, Q2, didx+2);
    g5sum -= C*getDeltaq(x, Q2, 4);

    return g5sum;
}

void eicPhysics::ReadPolTable(){
    fNx = 41;
    fNQ2 = 11;

    FILE *f = fopen("PolQuarkDist_fine.txt", "r");

    if( !f ){ printf("Can't open file PolQuarkDist_fine.txt\n"); exit(1);}

    // First three lines are bad
    
    int i,j,k;
    //char dummy[255];
    char *dummy2 = new char[255];
    double doubledummy;
    size_t dsize;

    //int nscan;
    getline(&dummy2, &dsize, f);
    getline(&dummy2, &dsize, f);
    getline(&dummy2, &dsize, f);

    for( i = 0; i < fNQ2; i++ ){
	for( j = 0; j < fNx; j++ ){
	    fscanf(f, "%lf%lf%lf%lf%lf%lf%lf", &fDeltaqQ2[i][j], &fDeltaqx[i][j],
	       &fDeltaq[0][i][j], &fDeltaq[1][i][j], &fDeltaq[2][i][j], &fDeltaq[3][i][j], &fDeltaq[4][i][j]);
	  //  printf("nscan= %d, x = %g, Q2 = %f  Dubar = %f\n", nscan, fDeltaqx[i][j], fDeltaqQ2[i][j], fDeltaq[0][i][j] );
	}
    }

    fclose(f);

    fNxlin = 20;
    fNQ2lin = 20;

    FILE *f2 = fopen("PolPdf_Linear.txt", "r");

    if( !f2 ){ printf("Can't open file PolPdf_Linear.txt\n"); exit(1);}

    // First three lines are bad
    
    for( i = 0; i < fNQ2lin; i++ ){
	for( j = 0; j < fNxlin-1; j++ ){
	    fscanf(f2, "%lf%lf%lf%lf%lf%lf%lf%lf", &fDeltaqxlin[i][j], &fDeltaqQ2lin[i][j],
	       &fDeltaqlin[0][i][j], &fDeltaqlin[1][i][j], &fDeltaqlin[2][i][j], &fDeltaqlin[3][i][j], &fDeltaqlin[4][i][j], &doubledummy);
	  //  printf("nscan= %d, x = %g, Q2 = %f  Dubar = %f\n", nscan, fDeltaqx[i][j], fDeltaqQ2[i][j], fDeltaq[0][i][j] );
	//    printf("%f %f %f %f %f %f %f\n", fDeltaqQ2lin[i][j], fDeltaqxlin[i][j],
	 //      fDeltaqlin[0][i][j], fDeltaqlin[1][i][j], fDeltaqlin[2][i][j], fDeltaqlin[3][i][j], fDeltaqlin[4][i][j]);
	}
    }

    // I guess constrain these to be 0 at 1
    for( i = 0; i < fNQ2lin; i++ ){
	for( j = 0; j < 5; j++ ){
	    fDeltaqlin[j][i][fNxlin-1] = 0.0;
	    fDeltaqxlin[i][fNxlin-1] = 1.0;
	    fDeltaqQ2lin[i][fNxlin-1] = fDeltaqQ2lin[i][fNxlin-2];
	    for( k = 0; k < fNxlin; k++ ){
		// This distribution is multiplied by x...
		fDeltaqlin[j][i][k] /= fDeltaqxlin[i][k];
	    }
	}
    }
    fclose(f2);

    printf("Read in linear\n");

    delete dummy2;
    return;
}

double eicPhysics::getDeltaq( double x, double Q2, int quark ){
    double xmin = log(1.0e-5);
    double xmax = log(1.0);
    double qmin = log(1.0);
    double qmax = log(1.0e5);

    double xlinmin = 0.05;
    double xlinmax = 1.0;
    double qlinmin = 1.0;
    double qlinmax = 20.0;

    double xscale, qscale;
    double xs, qs;
    int    xidx, qidx;
    double sum = 0.0;

    if( xlinmin < x && x < xlinmax && qlinmin < Q2 && Q2 < qlinmax ){
	// Use linearly spaced distributions
	// mostly for fixed target experiments
	sum = 0.0;

	xscale = (x  - xlinmin)*(fNxlin-1)/(xlinmax-xlinmin);
	qscale = (Q2 - qlinmin)*(fNQ2lin-1)/(qlinmax-qlinmin);

	xidx = (int) floor(xscale);
	qidx = (int) floor(qscale);

	xs = xscale - xidx;
	qs = qscale - qidx;

	// Bilinear interpolation

	//    printf("%g < x = %g < %g;  xidx = %d xs = %f\n", fDeltaqx[0][xidx], x,fDeltaqx[0][xidx+1], xidx, xs);
	//printf("%g < Q2 = %g < %g;  qidx = %d qs = %f\n", fDeltaqQ2[qidx][0], Q2,fDeltaqQ2[qidx+1][0], qidx, qs);

	sum = 0.0;

	sum += fDeltaqlin[quark][qidx][xidx]*(1.0-xs)*(1.0-qs);
	sum += fDeltaqlin[quark][qidx+1][xidx]*(1.0-xs)*qs;
	sum += fDeltaqlin[quark][qidx][xidx+1]*xs*(1.0-qs);
	sum += fDeltaqlin[quark][qidx+1][xidx+1]*xs*qs;
    } else {
	// logarithmically spaced distributions
	// mostly for collider experiments
    

	xscale = (log(x)  - xmin)*(fNx-1)/(xmax-xmin);
	qscale = (log(Q2) - qmin)*(fNQ2-1)/(qmax-qmin);

	xidx = (int) floor(xscale);
	qidx = (int) floor(qscale);

	xs = xscale - xidx;
	qs = qscale - qidx;

	// Bilinear interpolation

	//    printf("%g < x = %g < %g;  xidx = %d xs = %f\n", fDeltaqx[0][xidx], x,fDeltaqx[0][xidx+1], xidx, xs);
	//printf("%g < Q2 = %g < %g;  qidx = %d qs = %f\n", fDeltaqQ2[qidx][0], Q2,fDeltaqQ2[qidx+1][0], qidx, qs);

	sum = 0.0;

	sum += fDeltaq[quark][qidx][xidx]*(1.0-xs)*(1.0-qs);
	sum += fDeltaq[quark][qidx+1][xidx]*(1.0-xs)*qs;
	sum += fDeltaq[quark][qidx][xidx+1]*xs*(1.0-qs);
	sum += fDeltaq[quark][qidx+1][xidx+1]*xs*qs;
    }

    return sum;
}

// double eicPhysics::Wiser( int Z, int N, char PART, double E_IN, double P_IN, double TH_IN, double radlen_IN) {
// // Z, N are target info
// // PART is the particle whose cross section is to be calcualted
// // E_IN is beam energy in GeV
// // P_IN is scattered particle momentum in GeV
// // TH_IN is scattering angle in rad
// // radlen_IN is total radiation length in %
// // XSEC is the output cross section in nbarn/GeV/sr, per nucleon
// //    IMPORTANT: to get the cross section per nuclei, times this by (A)**0.8
// //
// // X. Zheng June 2006
// 
//     double E1, P, THP, rad_len, AN, XSEC, DLF, AL, QDF, PTP, E, TP, AJ, D2QD, D2QF, D2DEL, RADLEN;
//     int IA, IP, ITYPE;
//     
//     double PI=3.1416 ,AMP=938.28 ,AMD=1875.63 ,AMPI=139.6 ;
//     E1=E_IN*1000. ;
//     P=P_IN*1000. ;
//     THP=TH_IN*180./PI ;
//     rad_len=radlen_IN ;
// 
//     IA=Z+N ;
//     
//     if (PART=="P") {
//       AN=double(N)/3.+2.*double(Z)/3. ;
//       IP=1 ;
//       ITYPE=5 ;
//     }
//     else if (PART=="N") {
//       AN=double(Z)/3.+2.*double(N)/3.;
//       IP=-1; 
//     }
//     else if (PART=="PI+") { 
//       AN=double(Z)/3.;
//       IP=2 ; 
//       ITYPE=1;
//     }
//     else if (PART=="PI-") { 
//        AN=double(N)/3. ;
//        IP=2; 
//        ITYPE=2;
//     }
//     else if (PART=="PI0") { 
//        AN=2.*double(N+Z)/3. ; 
//        IP=2 ;
//     }
//     else {
//       printf("NOT PARTICLE INCLUDED IN THIS CROSS SECTION?!\n");
//       exit(1);
//     }
//     
//     if ( IP == 1 || IP == -1 ) {
//       if (IA > 1 && IA < 5) DLF=double(IA) ;
//       else  DLF=7. ;
//       AL=DLF ; // LEVINGER FACTOR 
//       QDF=AL*double(N*Z)/double(IA) ; 
//     }
//     
//     PTP=P;
//     
//     if ( IP == 1 || IP == -1 ) {
//       E=sqrt(pow(PTP,2)+pow(AMP,2)) ; 
//       TP=pow(PTP,2)/(E+AMP) ; 
//       AJ=PTP/E ;
//       D2QD=0. ; 
//       D2QF=0. ;
//     }
//     else {
//       E=sqrt(pow(PTP,2)+pow(AMPI,2)) ;
//       TP=pow(PTP,2)/(E+AMPI) ;
//       AJ=PTP/E ; 
//       D2QD=0. ;
//       D2QF=0. ;
//     }
//     D2DEL=0. ;
//     
//     if (IP==2 || IP==-2 || IP==0) {
//        XSEC = WISER_ALL_SIG(E1,PTP,THP,RAD_LEN,ITYPE) ;
//        // XSEC in nanobarn/GeV*str, for proton target
//        // should be timed by "effective number of nucleons per nuclei"???
//        // or (N+Z)???
//        //      XSEC=XSEC*((N+Z)**0.8) ! now is for per nuclei
//        
//     }
//     else XSEC = 0;
//     
//     
//     
//     return XSEC;
// }


// double eicPhysics::QUADMO(double &PLOWER,double &PUPPER,double &EPSLON, int &NLVL) { // In this code the fortran version
// // uses QUADMO just as a QUADMO(WISER_ALL_FIT,E_GAMMA_MIN,E08,EPSILON,N), in C WISER_ALL_FIT can be called inside 
// // the function, without declaring it in the argouments                                                   
//    int   LEVEL, MINLVL=3 ,MAXLVL=24, IRETRN[50],I, NLVL2 ;                 
//    double VALINT[50][2], VMX[50], RX[50], FMX[50], FRX[50], FMRX[50], ESTRX[50], EPSX[50] ;                                 
//    double  FL, FML, FM, FMR, FR, EST, ESTL, ESTR, ESTINT, DIFFEST, AREA, ABAREA, VM, COEF, ROMBRG ;                        
//    double VL, R, EPS;
//     LEVEL = 0  ;                                                   
//     NLVL2 = 0 ;                                                     
//     ABAREA = 0.0 ;                                                 
//     VL = PLOWER ;                                                     
//     R = PUPPER ;                                                    
//     FL = WISER_ALL_FIT(VL) ;                                                 
//     FM = WISER_ALL_FIT(0.5*(VL+R)) ;                                         
//     FR = WISER_ALL_FIT(R) ;                                                
//     EST = 0.0 ;                                                    
//     EPS = EPSLON ; 
    
// L100:
//     ++LEVEL;
//     VM = 0.5*(VL+R) ;
//     COEF = R-VL;
//     if (COEF != 0.) {
//       goto L150;
//     }
//     ROMBRG = EST ;
//     goto L300;
// L150:
//     FML = WISER_ALL_FIT(0.5*(VL+VM));                                           
//     FMR = WISER_ALL_FIT(0.5*(VM+R));                                           
//     ESTL = (FL+4.0*FML+FM)*COEF ;                                      
//     ESTR = (FM+4.0*FMR+FR)*COEF ;                                     
//     ESTINT = ESTL+ESTR ;
//     AREA=abs(ESTL)+abs(ESTR) ;                                       
//     ABAREA=AREA+ABAREA-abs(EST) ;
//     if (LEVEL != MAXLVL) {
// 	goto L200;
//     }
//     ++(NLVL2);
//     ROMBRG = ESTINT ;
//     goto L300;
// L200:
//     if ((DIFFEST = EST - ESTINT, abs(DIFFEST)) > EPS * ABAREA || LEVEL < MINLVL) {
// 	goto L400;
//     }
//     ROMBRG = (16.*ESTINT-EST)/15.0 ;
// L300:
//     --LEVEL;
//     I = IRETRN[LEVEL -1] ;
//     VALINT[LEVEL -1][I-1] = ROMBRG; 
//     switch (I) {
//       case 1: goto L500;
//       case 2: goto L600;
//     }
// L400:
//     IRETRN[LEVEL - 1] = 1 ;                                              
//     VMX[LEVEL - 1] = VM ;                                                
//     RX[LEVEL - 1] = R  ;                                               
//     FMX[LEVEL - 1] = FM ;                                              
//     FMRX[LEVEL - 1] = FMR ;                                             
//     FRX[LEVEL - 1] = FR ;                                               
//     ESTRX[LEVEL - 1] = ESTR ;                                          
//     EPSX[LEVEL - 1] = EPS ;                                            
//     EPS = EPS/1.4 ;
//     R = VM;                                                         
//     FR = FM ;                                                      
//     FM = FML ;                                                     
//     EST = ESTL ;
//     goto L100;
// L500:
//     IRETRN[LEVEL - 1] = 2 ;                                             
//     VL = VMX[LEVEL - 1] ;                                                 
//     R = RX[LEVEL - 1] ;                                                 
//     FL = FMX[LEVEL - 1] ;                                              
//     FM = FMRX[LEVEL - 1] ;                                              
//     FR = FRX[LEVEL - 1]  ;                                              
//     EST = ESTRX[LEVEL - 1] ;                                            
//     EPS = EPSX[LEVEL - 1] ;
//     goto L100;
// L600:
//     ROMBRG = VALINT[LEVEL-1][0]+VALINT[LEVEL][1] ;
//     if (LEVEL > 1) {
//       goto L300;
//     }
//     double ret_val = ROMBRG /12. ;
//     return ret_val;
    
// }


// double eicPhysics::WISER_ALL_FIT(double E_GAMMA) {
//     //---------------------------------------------------------
// // Calculates  pi, k, p  cross section for gamma + p -> k
// //  It is already divided by E_GAMMA, the bremstrulung spectra
// // David Wiser's fit from Thesis, eq. IV-A-2 and Table III.
// // Can be called from WISER_SIG using integration routine QUADMO
// // E,P are KAON energy and momentum
// // P_t is KAON transverse momentum
// // P_CM is KAON center of mass momentum
// // P_CM_L is KAON center of mass longitudinal momentum
// // TYPE:     1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar
// // E_GAMMA is photon energy.
// //             Steve Rock 2/21/96
// //---------------------------------------------------------
  
//   double  E,P,COST,P_T,M_X,U_MAN, ret_val ;
//   int  TYPE ;  //  1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar
//   int PARTICLE ;   // 1= pi, 2= K, 3 =P
// //  Wiser's fit    pi+     pi-    k+     k-     p+       p- 
//    double A1[6] = {566.,  486.,   368., 18.2,  1.33E5,  1.63E3 } ; 
//    double A2[6] = {829.,  115.,   1.91, 307.,  5.69E4, -4.30E3 } ; 
//    double A3[6] = {1.79,  1.77,   1.91, 0.98,  1.41,    1.79 }; 
//    double A4[6] = {2.10,  2.18,   1.15, 1.83,   .72,    2.24 } ;
//    double A6 = 1.90 ,A7 =-.0117;  //proton only
//    double MASS2[3] = {.019488, .2437, .8804} ;
//    double MASS[3] = {.1396, .4973, .9383} ; 
//    double MP2=.8804 ,MP=.9383 , RADDEG=.0174533 ;
//    double X_R, S,B_CM, GAM_CM,  P_CM ;
//    double P_CM_MAX, P_CM_L ;
//  //Mandlestam variables                                                
//    S = MP2 + 2.* E_GAMMA * MP ;   

// //Go to Center of Mass to get X_R
//    B_CM = E_GAMMA/(E_GAMMA+MP) ;
//    GAM_CM = 1./sqrt(1.-pow(B_CM,2)) ;
//    P_CM_L = -GAM_CM *B_CM *E + GAM_CM * P * COST ;
//    P_CM = sqrt(pow(P_CM_L,2) + pow(P_T,2)) ;  


//    P_CM_MAX =sqrt(S +pow((pow(M_X,2)-MASS2[PARTICLE-1]),2)/S -2.*(pow(M_X,2) +MASS2[PARTICLE-1]))/2. ;
//    X_R =  P_CM/P_CM_MAX ;  
//    if(X_R>1.) { // Out of kinematic range
//      ret_val = 0. ;
//    }
//    else if(TYPE != 5) { // not the proton
//      ret_val = (A1[TYPE-1] + A2[TYPE-1]/sqrt(S)) * pow((1. -X_R + pow(A3[TYPE-1],2)/S),A4[TYPE-1])/E_GAMMA ;
//    }
//    else  { // special formula for proton
//      ret_val = ( (A1[TYPE-1] + A2[TYPE-1]/sqrt(S)) * pow((1. -X_R + pow(A3[TYPE-1],2)/S),A4[TYPE-1]) / pow((1.+U_MAN),(A6+A7*S)) )/E_GAMMA  ;
//    }
      
//    return ret_val;
       
// }

// double eicPhysics::WISER_ALL_SIG(double E0,double P,double THETA_DEG,double RAD_LEN,double TYPE)  {
  
// //------------------------------------------------------------------------------
// // Calculate pi,K,p  cross section for electron beam on a proton target
// // IntegrateQs over function WISER_FIT using integration routine QUADMO
// // E0         is electron beam energy, OR max of Brem spectra
// // P,E       is scattered particle  momentum,energy
// // THETA_DEG  is kaon angle in degrees
// // RAD_LEN (%)is the radiation length of target, including internal
// //                (typically 5%)
// //               = .5 *(target radiation length in %) +5.
// //       ***  =100. IF BREMSTRULUNG PHOTON BEAM OF 1 EQUIVIVENT QUANTA ***
// // TYPE:     1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar
// // SIGMA      is output cross section in nanobars/GeV-str
// //------------------------------------------------------------------------------
    
//     double SIGMA;
//     double E,P_COM,COST,P_T,M_X,U_MAN ;
//     int TYPE_COM,PARTICLE ;
// //  Wiser's fit    pi+     pi-    k+     k-     p+      p-   
//     double A5[6] = {-5.49,  -5.23, -5.91, -4.45, -6.77,  -6.53} ;
//     double A6[6] = {-1.73,  -1.82, -1.74, -3.23,  1.90,  -2.45} ;
//     double MASS2[3] = {.019488, .2437, .8804} ;
//     double MASS[3] = {.1396, .4973, .9383} ; 
//     double MP=.9383 ,  MP2=.8804 , RADDEG=.0174533 ;
//     double  M_L,SIG_E ;
//     double E_GAMMA_MIN,E08,EPSILON=.003 ;                      
//     int N,CHARGE ;
    
//     P_COM = P  ;
//     TYPE_COM = int(TYPE) ;
//     PARTICLE = (TYPE_COM+1)/2 ; // 1= pi, 2= K, 3 =P
//     CHARGE = TYPE_COM -2*PARTICLE +2 ; // 1 for + charge, 2 for - charge
//     E08 =E0 ;
    
//     E =sqrt(MASS2[PARTICLE-1] + pow(P,2)) ;

//     COST = cos(RADDEG * THETA_DEG) ;
//     P_T = P * sin(RADDEG * THETA_DEG) ;
    
//     if(TYPE_COM<=4) {  //mesons
//        if(CHARGE==1) {   // K+ n final state
//         M_X = MP ;
//        }
//        else {   // K- K+ P final state
//         M_X = MP+ MASS[PARTICLE-1] ;
//        }
//     }
//     else {  // baryons 
//        if(CHARGE==1) {   // pi p  final state
//         M_X = MASS[0] ; // pion mass
//        }
//        else { // P P-bar  P final state
//         M_X = 2.*MP ;
//        }  
//     }
    
//     E_GAMMA_MIN = (pow(M_X,2) -MASS2[PARTICLE -1] - MP2+2.*MP*E)/ (2.*(MP -E +P*COST)) ;
    
//     if (E_GAMMA_MIN>.1){ //Kinematically allowed?
//        M_L = sqrt(pow(P_T,2) + MASS2[PARTICLE-1]) ;    

//        if (TYPE_COM != 5) {  // everything but proton
//         SIG_E = QUADMO(E_GAMMA_MIN,E08,EPSILON,N)  *  exp(A5[TYPE_COM-1] *M_L) * exp(A6[TYPE_COM-1] *pow(P_T,2)/E) ;
//        }
//        else { // proton 
//         U_MAN = abs(MP2 + MASS2[PARTICLE-1] -2.*MP*E) ;
//         SIG_E = QUADMO(E_GAMMA_MIN,E08,EPSILON,N)  *  exp(A5[TYPE_COM-1] *M_L) ; 
//        }
//        SIGMA = pow(P,2)/E * 1000. * RAD_LEN/100. *SIG_E ;
//     }
//     else {   // Kinematically forbidden
//        SIGMA = 0. ;
//     }    
  
//     return SIGMA;
// }



// double eicPhysics::Wiser_func_pip(double *x, double *par) {
//   // WISER_ALL_SIG(double E0,double P,double THETA_DEG,double RAD_LEN,int TYPE)
//   // E0 = par[0]
//   // RAD_LEN = par[1] 
//   // TYPE:     1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar
//   // P  = x[0]
//   // THETA_DEG = x[1]
//   double value = 0;
//   // need to modify TYPE
//   value = WISER_ALL_SIG(par[0],x[0],x[1],par[1],1);

//   return value;

// }

// double eicPhysics::Wiser_func_pim(double *x, double *par) {
//   // WISER_ALL_SIG(double E0,double P,double THETA_DEG,double RAD_LEN,int TYPE)
//   // E0 = par[0]
//   // RAD_LEN = par[1] 
//   // TYPE:     1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar
//   // P  = x[0]
//   // THETA_DEG = x[1]
//   double value = 0;
//   // need to modify TYPE
//   value = WISER_ALL_SIG(par[0],x[0],x[1],par[1],2);

//   return value;

// }

// double eicPhysics::Wiser_func_pi0(double *x, double *par) {
//   // WISER_ALL_SIG(double E0,double P,double THETA_DEG,double RAD_LEN,int TYPE)
//   // E0 = par[0]
//   // RAD_LEN = par[1] 
//   // TYPE:     1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar
//   // P  = x[0]
//   // THETA_DEG = x[1]
//   double value = 0;
//   // need to modify TYPE
//   value = 0.5 * (  WISER_ALL_SIG(par[0],x[0],x[1],par[1],1) + WISER_ALL_SIG(par[0],x[0],x[1],par[1],2)  );

//   return value;

// }


void eicPhysics::Decay_pi0(TVector3 vp, TVector3 vert) {
  double mp0 = 0.1349766 ; // mass pi0 in GeV
  TLorentzVector Vp_4(vp,sqrt(vp.Mag2() + pow(mp0,2))); // pi0 quadrimomentum
  TVector3 b_3 ; // beta to boost the LAB frame for going in the pi0 rest frame 
  b_3 = Vp_4.BoostVector(); // return (px/E,py/E,pz/E) (is all in GeV)

  // decay of the pi0 in two photons

  double x_pi=0.,y_pi=0.,z_pi=0.;
  fRandom->Sphere(x_pi,y_pi,z_pi,mp0/2); // In the pi0 rest frame the two photons are generated isotropically with half the energy of the pion (the mass)
  TVector3 g1(x_pi,y_pi,z_pi);
  TVector3 g2(1.,1.,1.); // two photons 3-momentum, set to 1 , so that they can be stretched and turned, setting the magnitude and the angles


  g2 = -g1; // still rest frame of the pion
  Gamma1.SetVect(g1); // defining the 4-vectors
  Gamma1.SetE(mp0/2);
  Gamma2.SetVect(g2);
  Gamma2.SetE(mp0/2);  
  Gamma1.Boost(b_3); // boosting in the LAB frame
  Gamma2.Boost(b_3);

  //define vertex of the two gammas .... guess one was fine too...
  TF1 *fr = new TF1("fr","exp(-x/(8.4e-17))",0,1.0e-14) ; // 8.4e-17 is the mean lifetime of the pi0
  double time = fr->GetRandom(0.,1.0e-14);
  TLorentzVector move(0.,0.,0.,time); // displacement for the creation of the two gammas in the pi0 rest frame
  move.Boost(b_3); // displacement for the creation of the two gammas in the LAB frame
  Gamma1_vt.SetX( vert.X() + move.X() );
  Gamma1_vt.SetY( vert.Y() + move.Y() );
  Gamma1_vt.SetZ( vert.Z() + move.Z() );

  delete fr;
   
  return;

}
