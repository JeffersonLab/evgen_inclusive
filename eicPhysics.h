#ifndef __eicPhysics_h
#define __eicPhysics_h

#include "cteqpdf.h"

#include "eicBeam.h"
#include "eicIon.h"
#include "eicEvent.h"
#include "eicModel.h"
#include <iostream>


#include "TRandom2.h"

#define MASS_P 0.938


class eicPhysics {
    public:
	eicPhysics();
	~eicPhysics();

	void MakeEvent( eicBeam *, eicIon *, eicEvent * );
	void MakeEvent2( eicBeam *, eicIon *, eicEvent *, eicModel *);

	enum nucl {kProton, kNeutron};
    private:
	 cteq_pdf_t *fPDF;

	 TRandom2 *fRandom;

	 double GetTheta( double min = 0, double max = 3.14159 );
	 double GetPhi(   double min = 0, double max = 3.14159*2 );
	 double GetEf( double, double );
	 double GetEfmax( double, double );

	 double F2(double x, double Q2, nucl n );
	 double F1(double x, double Q2, nucl n );

	 double F1gZ(double x, double Q2, nucl n );
	 double F3gZ(double x, double Q2, nucl n );
	 double g1gZ(double x, double Q2, nucl n );
	 double g5gZ(double x, double Q2, nucl n );

	 double F1Z(double x, double Q2, nucl n );
	 double F3Z(double x, double Q2, nucl n );
	 double g1Z(double x, double Q2, nucl n );
	 double g5Z(double x, double Q2, nucl n );


	 double F1W(double x, double Q2, nucl n, int C );
	 double F2W(double x, double Q2, nucl n, int C );
	 double F3W(double x, double Q2, nucl n, int C );
	 double g1W(double x, double Q2, nucl n, int C );
	 double g4W(double x, double Q2, nucl n, int C );
	 double g5W(double x, double Q2, nucl n, int C );

	 void ReadPolTable();
	 double getDeltaq( double x, double Q2, int quark );

	 int fNx;
	 int fNQ2;
	 double fDeltaq[5][50][50];
	 double fDeltaqx[50][50];
	 double fDeltaqQ2[50][50];

	 int fNxlin;
	 int fNQ2lin;
	 // Linearly spaced pol pdfs
	 double fDeltaqlin[5][50][50];
	 double fDeltaqxlin[50][50];
	 double fDeltaqQ2lin[50][50];
	 
	 // Pion part
//	 double Wiser( int Z, int N, char PART, double E_IN, double P_IN, double TH_IN, double radlen_IN);

	 double QUADMO(double &PLOWER,double &PUPPER,double &EPSLON, int &NLVL);
	 double WISER_ALL_FIT(double E_GAMMA);
	 double WISER_ALL_SIG(double E0,double P,double THETA_DEG,double RAD_LEN,double TYPE);
	 double Wiser_func_pip(double *x, double *par);
	 double Wiser_func_pim(double *x, double *par);
	 double Wiser_func_pi0(double *x, double *par);
	 
};
#endif//__eicPhysics_h
