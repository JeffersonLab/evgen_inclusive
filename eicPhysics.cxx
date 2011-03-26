#include "eicPhysics.h"

#include "cteqpdf.h"

#include <stdio.h>
#include <math.h>

#include "TRandom.h"
#include "TLorentzVector.h"


eicPhysics::eicPhysics(){
    int mode = 400;
    fPDF = cteq_pdf_alloc_id(mode);
    
    fRandom = new TRandom2(0);
    printf("Seed number %d\n",fRandom->GetSeed());
    ReadPolTable();

    return;
}

eicPhysics::~eicPhysics(){
    cteq_pdf_free(fPDF);
    return;
}

void eicPhysics::MakeEvent(eicBeam *beam, eicIon *ion, eicEvent *ev ){
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
    double s     = (ki+pi).Mag2();

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


    // From Halzen and Marten
//    double ds_dxdy 
//	= 2.0*3.14159*(x*y*y*f1 + ((1.0-y)-M*x*y/(2.0*numax))*f2)/(x*x*y*y*M*numax*137.0*137.0);

    //  From PDG
    double ds_dxdy 
	= 4.0*3.14159*((1.0-y-pow(x*y*MASS_P,2.0)/Q2)*f2+y*y*x*f1)
	  /(x*y*Q2*137.0*137.0);


//    double ds_dOmegadE = ds_dxdy*ef*J/(2.0*M*e_rest.Energy()*3.14159*y);

    ////////////////////////////////////////////////////
    //  Parity violation
    double G_fermi = 1.16637e-5; // GeV^-2
    double MZ = 91.1876; // GeV
    double sin2thW = 0.23119;

    double Ynum = 2.0*x*(1.0-y/2.0);
    double Yden = 2.0*x*y + 4.0*x*(1.0-y-x*y*M/(2.0*numax))/y;

    double Y = Ynum/Yden;

    double eta_gZ = G_fermi*MZ*MZ*Q2*137.0/(2.0*sqrt(2.0)*3.14159)/(Q2+MZ*MZ);

    double gA = -1.0/2.0;
    double gV = -1.0/2.0 + 2.0*sin2thW;

    double Abeam = -eta_gZ*(gA*f1gz/f1 + Y*gV*f3gz/f1);

    double g1gz = A*g1gZ(x,Q2,n);
    double g5gz = A*g5gZ(x,Q2,n);

    double eta_L = eta_gZ*(x*y*y + x*(1.0-y)*(2.0-x*y*M/numax))
	           /(x*y*y + 2.0*x*(1.0-y-x*y*M/(2.0*numax)));
    double Y_L   = y*(2.0-y-x*y*M/numax)/(y*y + (1.0-y)*(2.0-x*y*M/numax));

    double eta_T = 0.0;
    double Y_T   = 0.0;

    double A_L = eta_L*(gV*g5gz/f1 + Y_L*gA*g1gz/f1);
    double A_T = 0.0;

    ////////////////////////////////////////////////////

    // Boost kf back to EIC frame
    kf.Boost(-blab);

    double th_eic = kf.Theta();

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

    if( 0.0 < ds_dxdy && ds_dxdy < 1e9 ){
	data.weight  = ds_dxdydphi*(0.197*0.197*1e-30)*V*J; // GeV^-2 to m^2
	data.weight *= beam->GetLumin();
    } else {
	// Unphysical for some reason
	data.weight = 0.0;
    }

    data.rate   = 0.0;

//    data.crs    = ds_dOmegadE;
    data.crs    = ds_dxdy;

    data.ef     = kf.Energy();
    data.theta  = kf.Theta();

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

    data.A_T     = A_T;
    data.eta_T   = eta_T;
    data.Y_T     = Y_T;

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

double eicPhysics::F1( double x, double Q2, nucl n ){
    return F2(x,Q2,n)/(2.0*x);
}

double eicPhysics::Dp( double x) {
  if (x<0.3 || x>1.0)   return 0.0;
  if (x>=0.3 && x<0.6)  return 1.0-6.083333*x+9.16666666*pow(x,2);
  else return 9.149226-2.991434e+01*x+2.610939e+01*pow(x,2);
}

double eicPhysics::Dn( double x) {
  if (x<0.3 || x>1.0)   return 0.0;
  else return 2*(1.771210e-01+2.659544*x-1.596289e+01*pow(x,2)+2.055688e+01*pow(x,3))-Dp(x);
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
    double DHT = 0.0;
    int uidx, didx;

    switch( n ){
	case kProton:
	    uidx = 1;
	    didx = 2;
	    DHT = Dp(x);
	    break;
	case kNeutron:
    	    uidx = 2;
	    didx = 1;
	    DHT = Dn(x);
	    break;
	default:
    	    uidx = 1;
	    didx = 2;
	    DHT = Dp(x);
	    break;
    }

    f2sum += e2_u*x*cteq_pdf_evolvepdf(fPDF, uidx, x, sqrt(Q2));
    f2sum += e2_d*x*cteq_pdf_evolvepdf(fPDF, didx, x, sqrt(Q2));

    // Strange quark
    f2sum += e2_d*x*cteq_pdf_evolvepdf(fPDF,  3, x, sqrt(Q2));
    f2sum += e2_d*x*cteq_pdf_evolvepdf(fPDF, -3, x, sqrt(Q2));

    // ubar
    f2sum += e2_u*x*cteq_pdf_evolvepdf(fPDF, -uidx, x, sqrt(Q2));
    // dbar
    f2sum += e2_d*x*cteq_pdf_evolvepdf(fPDF, -didx, x, sqrt(Q2));

    f2sum = (1+DHT/Q2) * f2sum;

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

    // ubar
    f3sum -= e_u*gA_u*cteq_pdf_evolvepdf(fPDF, -uidx, x, sqrt(Q2));
    // dbar
    f3sum -= e_d*gA_d*cteq_pdf_evolvepdf(fPDF, -didx, x, sqrt(Q2));

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

    g5sum += e_u*gA_u*getDeltaq(x, Q2, uidx);
    g5sum += e_d*gA_d*getDeltaq(x, Q2, didx);

    // ubar
    g5sum -= e_u*gA_u*getDeltaq(x, Q2, uidx+2);
    // dbar
    g5sum -= e_d*gA_d*getDeltaq(x, Q2, didx+2);

    return g5sum;
}


void eicPhysics::ReadPolTable(){
    fNx = 41;
    fNQ2 = 11;

    FILE *f = fopen("PolQuarkDist_fine.txt", "r");

    if( !f ){ printf("Can't open file PolQuarkDist_fine.txt\n"); exit(1);}

    // First three lines are bad
    
    int i,j;
    char dummy[255];
    char *dummy2 = new char[255];
    size_t dsize;

    int nscan;
    getline(&dummy2, &dsize, f);
    getline(&dummy2, &dsize, f);
    getline(&dummy2, &dsize, f);

    for( i = 0; i < fNQ2; i++ ){
	for( j = 0; j < fNx; j++ ){
	    nscan = fscanf(f, "%lf%lf%lf%lf%lf%lf%lf", &fDeltaqQ2[i][j], &fDeltaqx[i][j],
	       &fDeltaq[0][i][j], &fDeltaq[1][i][j], &fDeltaq[2][i][j], &fDeltaq[3][i][j], &fDeltaq[4][i][j]);
//	    printf("nscan= %d, x = %g, Q2 = %f  Dubar = %f\n", nscan, fDeltaqx[i][j], fDeltaqQ2[i][j], fDeltaq[0][i][j] );
	}
    }

    delete dummy2;

    return;
}

double eicPhysics::getDeltaq( double x, double Q2, int quark ){
    double xmin = log(1.0e-5);
    double xmax = log(1.0);
    double qmin = log(1.0);
    double qmax = log(1.0e5);

    double xscale, qscale;
    double xs, qs;
    int    xidx, qidx;

    xscale = (log(x)  - xmin)*(fNx-1)/(xmax-xmin);
    qscale = (log(Q2) - qmin)*(fNQ2-1)/(qmax-qmin);

    xidx = (int) floor(xscale);
    qidx = (int) floor(qscale);

    xs = xscale - xidx;
    qs = qscale - qidx;

    // Bilinear interpolation

//    printf("%g < x = %g < %g;  xidx = %d xs = %f\n", fDeltaqx[0][xidx], x,fDeltaqx[0][xidx+1], xidx, xs);
    //printf("%g < Q2 = %g < %g;  qidx = %d qs = %f\n", fDeltaqQ2[qidx][0], Q2,fDeltaqQ2[qidx+1][0], qidx, qs);

    double sum = 0.0;

    sum += fDeltaq[quark][qidx][xidx]*(1.0-xs)*(1.0-qs);
    sum += fDeltaq[quark][qidx+1][xidx]*(1.0-xs)*qs;
    sum += fDeltaq[quark][qidx][xidx+1]*xs*(1.0-qs);
    sum += fDeltaq[quark][qidx+1][xidx+1]*xs*qs;

    return sum;
}
